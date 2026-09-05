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
# mixing parameter phi is interpretable - it reports what share of the
# residual variance is spatially structured, which is itself evidence about
# whether the highland paradox is a geographic phenomenon or a compositional
# one.
#
# Inputs : data/processed/direct_estimates_adm2.csv  (from 03)
#          data/processed/adm2_districts.gpkg        (from 01)
#          data/processed/commune_covariates.csv     (from 02, aggregated up)
# Outputs: data/processed/summer_adm2.csv, outputs/figures/04_*.png
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("dplyr", "sf", "ggplot2", "scales", "spdep", "SUMMER", "INLA")

set.seed(20210)

# ===========================================================================
# 1. Inputs
# ===========================================================================
stopifnot(file.exists(OUT$direct_adm2), file.exists(OUT$adm2))
direct <- utils::read.csv(OUT$direct_adm2)
adm2 <- sf::st_read(OUT$adm2, quiet = TRUE)

msg("districts in the boundary file: ", nrow(adm2),
    " | with a direct estimate: ", nrow(direct))

# Every district must appear in the model frame, including those the DHS never
# sampled: those are exactly the ones the model exists to predict.
frame <- adm2 |>
  sf::st_drop_geometry() |>
  dplyr::select(ADM2_PCODE, ADM2_EN, ADM1_PCODE, ADM1_EN) |>
  dplyr::left_join(direct, by = "ADM2_PCODE")

msg("districts with no DHS data (predicted from the spatial prior alone): ",
    sum(is.na(frame$direct)))

# ===========================================================================
# 2. District adjacency
# ===========================================================================
nb2 <- spdep::poly2nb(adm2, row.names = adm2$ADM2_PCODE, queen = TRUE)
iso <- which(spdep::card(nb2) == 0)
if (length(iso) > 0) {
  cent <- sf::st_coordinates(sf::st_point_on_surface(
    sf::st_geometry(sf::st_transform(adm2, CRS_EQA))))
  for (i in iso) {
    d <- sqrt((cent[, 1] - cent[i, 1])^2 + (cent[, 2] - cent[i, 2])^2); d[i] <- Inf
    j <- which.min(d)
    nb2[[i]] <- sort(unique(c(as.integer(nb2[[i]][nb2[[i]] > 0L]), j)))
    nb2[[j]] <- sort(unique(c(as.integer(nb2[[j]][nb2[[j]] > 0L]), i)))
  }
  msg("connected ", length(iso), " isolated districts")
}
Amat <- spdep::nb2mat(nb2, style = "B", zero.policy = TRUE)
dimnames(Amat) <- list(adm2$ADM2_PCODE, adm2$ADM2_PCODE)

# ===========================================================================
# 3. District-level covariates (population-weighted from the commune stack)
# ===========================================================================
# Covariates enter on the logit scale as fixed effects. They are aggregated
# from communes with under-5 population weights, so a district's covariate is
# the value experienced by its average under-5 child rather than by its
# average square kilometre.
COVARS <- c("elevation", "ruggedness", "temp_min_cold", "temp_seasonality",
            "precip_annual", "precip_seasonality", "frac_crop", "frac_built",
            "frac_water_perm", "cattle_density", "travel_time", "nightlights")

X <- NULL
if (file.exists(OUT$cov_commune)) {
  cov_com <- utils::read.csv(OUT$cov_commune)
  have <- intersect(COVARS, names(cov_com))
  msg("covariates available: ", paste(have, collapse = ", "))

  X <- cov_com |>
    dplyr::mutate(w = pmax(pop_u5, 1e-6)) |>
    dplyr::group_by(ADM2_PCODE) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(have),
                                   ~ stats::weighted.mean(.x, w, na.rm = TRUE)),
                     pop_u5 = sum(pop_u5, na.rm = TRUE), .groups = "drop")

  # Standardise so the priors on beta are on a common scale and the
  # coefficients are directly comparable in the effect table.
  X <- X |> dplyr::mutate(dplyr::across(dplyr::all_of(have),
                                        ~ as.numeric(scale(.x))))
  frame <- dplyr::left_join(frame, X, by = "ADM2_PCODE")
} else {
  msg("NOTE: no covariate file yet - fitting the spatial-only model.")
}

# ===========================================================================
# 4. Smoothed-direct model
# ===========================================================================
fit_data <- frame |>
  dplyr::transmute(region = ADM2_PCODE,
                   direct.est = direct,
                   direct.var = direct_var,
                   direct.logit.est = logit_direct,
                   direct.logit.var = logit_var)

have_cov <- if (is.null(X)) character(0) else intersect(COVARS, names(frame))
if (length(have_cov) > 0) {
  # Drop covariates that are missing anywhere: INLA cannot use a fixed effect
  # with an NA in the prediction frame.
  keep <- have_cov[vapply(have_cov, function(v) !anyNA(frame[[v]]), logical(1))]
  if (length(keep) < length(have_cov))
    msg("dropping covariates with missing values: ",
        paste(setdiff(have_cov, keep), collapse = ", "))
  have_cov <- keep
  fit_data <- cbind(fit_data, frame[, have_cov, drop = FALSE])
}

msg("fitting SUMMER smoothDirect (BYM2, ", length(have_cov), " covariates)")

fml <- if (length(have_cov) > 0)
  stats::as.formula(paste("~", paste(have_cov, collapse = " + "))) else NULL

fit <- SUMMER::smoothDirect(
  data = fit_data,
  Amat = Amat,
  X = if (length(have_cov) > 0) cbind(region = fit_data$region,
                                      fit_data[, have_cov, drop = FALSE]) else NULL,
  formula = fml,
  year_label = NULL,
  time.model = NULL,          # single cross-section: no temporal component
  spatial.model = "bym2",
  responseType = "binary",    # models the logit of the direct estimate
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

res <- SUMMER::getSmoothed(fit)
if (is.list(res) && !is.data.frame(res)) res <- res[[1]]

out <- res |>
  dplyr::transmute(ADM2_PCODE = as.character(region),
                   summer_est = median,
                   summer_lower = lower,
                   summer_upper = upper,
                   summer_var = (upper - lower)^2 / (2 * 1.96)^2) |>
  dplyr::left_join(frame[, c("ADM2_PCODE", "ADM2_EN", "ADM1_EN", "direct",
                             "direct_var", "se", "n_children", "n_clusters")],
                   by = "ADM2_PCODE")

out$se_reduction <- 1 - sqrt(out$summer_var) / out$se

utils::write.csv(out, file.path(DIR$processed, "summer_adm2.csv"), row.names = FALSE)

msg("median SE reduction vs the direct estimator: ",
    round(100 * stats::median(out$se_reduction, na.rm = TRUE), 1), "%")

# ===========================================================================
# 5. Model summary
# ===========================================================================
sm <- summary(fit)
utils::capture.output(sm, file = file.path(DIR$tables, "04_summer_summary.txt"))

# phi = share of the residual variance that is spatially structured.
if (!is.null(fit$fit$summary.hyperpar)) {
  hp <- fit$fit$summary.hyperpar
  utils::write.csv(hp, file.path(DIR$tables, "04_summer_hyperpar.csv"))
  phi_row <- grep("Phi", rownames(hp))
  if (length(phi_row) > 0)
    msg("BYM2 phi = ", round(hp[phi_row[1], "0.5quant"], 3),
        " (share of residual variation that is spatially structured)")
}
if (!is.null(fit$fit$summary.fixed)) {
  utils::write.csv(fit$fit$summary.fixed,
                   file.path(DIR$tables, "04_summer_fixed_effects.csv"))
}

# ===========================================================================
# 6. Figures
# ===========================================================================
map_df <- dplyr::left_join(adm2, out, by = "ADM2_PCODE")
theme_map <- ggplot2::theme_void(base_size = 9) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11))

p1 <- ggplot2::ggplot(map_df) +
  ggplot2::geom_sf(ggplot2::aes(fill = summer_est), colour = "grey30", linewidth = 0.1) +
  ggplot2::scale_fill_viridis_c(option = "rocket", direction = -1, labels = scales::percent,
                                name = "stunting", limits = c(0, NA)) +
  ggplot2::labs(title = "Smoothed-direct stunting prevalence, district level",
                subtitle = "SUMMER BYM2, DHS 2021") + theme_map
save_fig(p1, "04_summer_adm2_estimate.png")

p2 <- ggplot2::ggplot(out, ggplot2::aes(direct, summer_est)) +
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
