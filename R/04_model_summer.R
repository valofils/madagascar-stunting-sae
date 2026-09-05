# ---------------------------------------------------------------------------
# 04_model_summer.R
#
# Area-level Bayesian small-area estimation with SUMMER: the smoothed-direct
# (Fay-Herriot-with-spatial-structure) estimator of Mercer et al. (2015) and
# Wakefield et al. (2020), fitted at district (adm2) level with a BYM2 random
# effect on the district adjacency graph.
#
# Model, on the logit scale:
#   logit(p_i)_direct ~ Normal(mu_i, V_i)          V_i = design variance (known)
#   mu_i = alpha + x_i' beta + b_i
#   b_i  = BYM2(phi, tau): phi * spatial (ICAR, scaled) + (1 - phi) * iid
#
# BYM2 (Riebler et al. 2016) is used rather than plain BYM because its single
# mixing parameter phi is interpretable - it reports what share of the residual
# variance is spatially structured, which is itself evidence about whether the
# highland paradox is a geographic phenomenon or a compositional one.
#
# Implementation note: SUMMER 2.0's smoothSurvey() is the right entry point
# here, not smoothDirect(). smoothDirect() is built around the mortality
# time-series case and requires year.label / year.range; smoothSurvey() takes
# the child-level survey data plus a design specification for a single
# cross-section, computes the weighted direct estimates internally, and fits
# the BYM2 smoothing model to them. It is also the function that accepts areal
# covariates through its X argument.
#
# Inputs : data/interim/dhs_child_haz.rds          (from 03)
#          data/processed/adm2_districts.gpkg      (from 01)
#          data/processed/commune_covariates.csv   (from 02, aggregated up)
# Outputs: data/processed/summer_adm2.csv, outputs/figures/04_*.png
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("dplyr", "sf", "ggplot2", "scales", "spdep", "SUMMER", "INLA")

set.seed(20210)

# ===========================================================================
# 1. Inputs
# ===========================================================================
stopifnot(file.exists(OUT$dhs_child), file.exists(OUT$adm2))
child <- readRDS(OUT$dhs_child)
adm2 <- sf::st_read(OUT$adm2, quiet = TRUE)

child <- child[!is.na(child$ADM2_PCODE) & !is.na(child$stunted), ]
msg("districts in the boundary file: ", nrow(adm2),
    " | districts with DHS children: ", dplyr::n_distinct(child$ADM2_PCODE))
msg("children in the model: ", nrow(child))

# ===========================================================================
# 2. District adjacency
# ===========================================================================
# The row/column names of Amat define the region list, so every district ends
# up in the output - including those the DHS never sampled, which are exactly
# the ones the model exists to predict.
nb2 <- spdep::poly2nb(adm2, row.names = adm2$ADM2_PCODE, queen = TRUE)

cent <- sf::st_coordinates(sf::st_point_on_surface(
  sf::st_geometry(sf::st_transform(adm2, CRS_EQA))))
add_edge <- function(nb, i, j) {
  nb[[i]] <- sort(unique(c(as.integer(nb[[i]][nb[[i]] > 0L]), j)))
  nb[[j]] <- sort(unique(c(as.integer(nb[[j]][nb[[j]] > 0L]), i)))
  nb
}
# As in 01: an ICAR prior on a disconnected graph is improper, so connect the
# island districts to the mainland before building the adjacency matrix.
repeat {
  comp <- spdep::n.comp.nb(nb2)
  if (comp$nc <= 1) break
  sizes <- table(comp$comp.id)
  main <- as.integer(names(sizes)[which.max(sizes)])
  inm <- which(comp$comp.id == main); out <- which(comp$comp.id != main)
  d <- outer(seq_along(out), seq_along(inm), function(a, b) {
    sqrt((cent[out[a], 1] - cent[inm[b], 1])^2 + (cent[out[a], 2] - cent[inm[b], 2])^2)
  })
  k <- which(d == min(d), arr.ind = TRUE)[1, ]
  nb2 <- add_edge(nb2, out[k[1]], inm[k[2]])
}
Amat <- spdep::nb2mat(nb2, style = "B", zero.policy = TRUE)
dimnames(Amat) <- list(adm2$ADM2_PCODE, adm2$ADM2_PCODE)
msg("adjacency: ", nrow(Amat), " districts, mean neighbours ",
    round(mean(rowSums(Amat)), 2))

# ===========================================================================
# 3. District-level covariates (population-weighted from the commune stack)
# ===========================================================================
# Covariates enter on the logit scale as fixed effects. They are aggregated
# from communes with under-5 population weights, so a district's covariate is
# the value experienced by its average under-5 child rather than by its average
# square kilometre.
COVARS <- c("elevation", "ruggedness", "temp_min_cold", "temp_seasonality",
            "precip_annual", "precip_seasonality", "frac_crop", "frac_built",
            "frac_water_perm", "cattle_density", "travel_time",
            "travel_time_healthcare", "nightlights")

X <- NULL
if (file.exists(OUT$cov_commune)) {
  cov_com <- utils::read.csv(OUT$cov_commune)
  have <- intersect(COVARS, names(cov_com))

  X <- cov_com |>
    dplyr::mutate(w = pmax(pop_u5, 1e-6)) |>
    dplyr::group_by(ADM2_PCODE) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(have),
                                   ~ stats::weighted.mean(.x, w, na.rm = TRUE)),
                     .groups = "drop")

  # Drop covariates that are missing for any district: SUMMER passes X straight
  # to INLA, which cannot use a fixed effect with an NA in the prediction frame.
  keep <- have[vapply(have, function(v) !anyNA(X[[v]]), logical(1))]
  if (length(keep) < length(have))
    msg("dropping covariates with missing districts: ",
        paste(setdiff(have, keep), collapse = ", "))

  # Standardise so the coefficients are comparable and the default priors sane.
  X <- X |> dplyr::mutate(dplyr::across(dplyr::all_of(keep),
                                        ~ as.numeric(scale(.x))))
  X <- as.data.frame(X[, c("ADM2_PCODE", keep)])
  # SUMMER links X to the data by matching a column of X against regionVar,
  # so the key column has to be named "region", not ADM2_PCODE.
  names(X)[1] <- "region"
  have <- keep
  msg("covariates in the model (", length(have), "): ",
      paste(have, collapse = ", "))
} else {
  have <- character(0)
  msg("NOTE: no covariate file yet - fitting the spatial-only model.")
}

# ===========================================================================
# 4. Smoothed-direct model
# ===========================================================================
# smoothSurvey computes the design-based direct estimates internally from
# strata / cluster / weights, then smooths them with the BYM2 prior. Passing
# the child-level data rather than precomputed direct estimates lets SUMMER
# keep the design and the smoothing consistent with each other.
fit_data <- child |>
  dplyr::transmute(region = ADM2_PCODE,
                   stunted = stunted,
                   strata = strata,
                   weights = wt,
                   cluster = cluster,
                   household = household)

msg("fitting SUMMER smoothSurvey (BYM2, ", length(have), " covariates)")
fit <- SUMMER::smoothSurvey(
  data = fit_data,
  Amat = Amat,
  X = X,
  response.type = "binary",
  responseVar = "stunted",
  strataVar = "strata",
  weightVar = "weights",
  regionVar = "region",
  clusterVar = "~cluster+household",
  CI = 0.95,
  save.draws = TRUE)

saveRDS(fit, file.path(DIR$interim, "04_summer_fit.rds"))

# ===========================================================================
# 5. Collect results
# ===========================================================================
sm <- fit$smooth
ht <- fit$HT

pick <- function(df, cands, default = NA_real_) {
  hit <- cands[cands %in% names(df)]
  if (length(hit) == 0) return(rep(default, nrow(df))) else df[[hit[1]]]
}

out <- data.frame(
  ADM2_PCODE   = as.character(pick(sm, c("region"))),
  summer_est   = pick(sm, c("mean")),
  summer_median = pick(sm, c("median")),
  summer_lower = pick(sm, c("lower")),
  summer_upper = pick(sm, c("upper")),
  summer_var   = pick(sm, c("var")))
# SUMMER reports var on the requested scale; fall back to the interval width.
bad <- is.na(out$summer_var)
out$summer_var[bad] <- ((out$summer_upper - out$summer_lower)[bad] / (2 * 1.96))^2

direct <- data.frame(
  ADM2_PCODE = as.character(pick(ht, c("region"))),
  direct     = pick(ht, c("direct.est", "HT.est")),
  direct_var = pick(ht, c("direct.var", "HT.var")))
direct$se <- sqrt(direct$direct_var)

counts <- child |>
  dplyr::group_by(ADM2_PCODE) |>
  dplyr::summarise(n_children = dplyr::n(),
                   n_clusters = dplyr::n_distinct(cluster), .groups = "drop")

out <- out |>
  dplyr::left_join(direct, by = "ADM2_PCODE") |>
  dplyr::left_join(counts, by = "ADM2_PCODE") |>
  dplyr::left_join(sf::st_drop_geometry(adm2)[, c("ADM2_PCODE", "ADM2_EN", "ADM1_EN")],
                   by = "ADM2_PCODE") |>
  dplyr::mutate(se_reduction = 1 - sqrt(summer_var) / se)

utils::write.csv(out, file.path(DIR$processed, "summer_adm2.csv"), row.names = FALSE)

msg("districts returned: ", nrow(out),
    " | with a direct estimate: ", sum(!is.na(out$direct)))
msg("median SE reduction vs the direct estimator: ",
    round(100 * stats::median(out$se_reduction, na.rm = TRUE), 1), "%")

# ===========================================================================
# 6. Model summary
# ===========================================================================
utils::capture.output(summary(fit),
                      file = file.path(DIR$tables, "04_summer_summary.txt"))

# phi = share of the residual variance that is spatially structured. A phi near
# 1 says the district-level residual is almost entirely spatially smooth, which
# is direct evidence that the paradox is geographic rather than compositional.
if (!is.null(fit$fit$summary.hyperpar)) {
  hp <- fit$fit$summary.hyperpar
  utils::write.csv(hp, file.path(DIR$tables, "04_summer_hyperpar.csv"))
  ph <- grep("Phi", rownames(hp))
  if (length(ph) > 0)
    msg("BYM2 phi = ", round(hp[ph[1], "0.5quant"], 3),
        " [", round(hp[ph[1], "0.025quant"], 3), ", ",
        round(hp[ph[1], "0.975quant"], 3),
        "] - share of residual variation that is spatially structured")
}
if (!is.null(fit$fit$summary.fixed)) {
  fx <- fit$fit$summary.fixed
  utils::write.csv(fx, file.path(DIR$tables, "04_summer_fixed_effects.csv"))
  msg("fixed effects written (", nrow(fx), " terms)")
}

# ===========================================================================
# 7. Figures
# ===========================================================================
map_df <- dplyr::left_join(adm2, out, by = c("ADM2_PCODE", "ADM2_EN", "ADM1_EN"))
theme_map <- ggplot2::theme_void(base_size = 9) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11))

p1 <- ggplot2::ggplot(map_df) +
  ggplot2::geom_sf(ggplot2::aes(fill = summer_est), colour = "grey30",
                   linewidth = 0.1) +
  ggplot2::scale_fill_viridis_c(option = "rocket", direction = -1,
                                labels = scales::percent, name = "stunting") +
  ggplot2::labs(title = "Smoothed-direct stunting prevalence, district level",
                subtitle = "SUMMER BYM2, DHS 2021") + theme_map
save_fig(p1, "04_summer_adm2_estimate.png")

p2 <- ggplot2::ggplot(out[!is.na(out$direct), ],
                      ggplot2::aes(direct, summer_est)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = summer_lower, ymax = summer_upper),
                         width = 0, colour = "grey70") +
  ggplot2::geom_point(ggplot2::aes(size = n_clusters), alpha = 0.6) +
  ggplot2::scale_x_continuous(labels = scales::percent) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::labs(x = "Direct estimate", y = "Smoothed estimate", size = "clusters",
                title = "Shrinkage towards the spatial mean",
                subtitle = "Districts with fewest clusters are pulled furthest") +
  ggplot2::theme_minimal(base_size = 10)
save_fig(p2, "04_summer_shrinkage.png", width = 6.5, height = 5.5)

msg("04_model_summer.R complete")
