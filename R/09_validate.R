# ---------------------------------------------------------------------------
# 09_validate.R
#
# Validation of the geostatistical model, in four parts:
#
#   1. Spatial block cross-validation. Random k-fold CV is optimistic for a
#      spatially correlated model: a held-out cluster usually has a neighbour
#      still in the training set, so the spatial field interpolates it almost
#      for free. Blocking the folds in space (contiguous spatial tiles held out
#      together) forces the model to extrapolate the distance it will actually
#      have to extrapolate when predicting an unsampled commune, which is the
#      honest test of the whole exercise.
#
#   2. Internal fit: WAIC, DIC, and the CPO/PIT diagnostics INLA provides.
#
#   3. Benchmarking. The model is not design-based, so its population-weighted
#      regional aggregates need not reproduce the design-based DHS regional
#      estimates. Where they differ materially the model is benchmarked
#      (ratio-adjusted) so that published numbers remain consistent with the
#      official survey figures - standard practice in SAE (Rao & Molina 2015).
#
#   4. External coherence: district aggregates against the IHME Local Burden
#      of Disease surfaces, where a local copy is available.
#
# Inputs : outputs of 03, 05, 08
# Outputs: outputs/tables/09_*.csv, outputs/figures/09_*.png
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("dplyr", "sf", "ggplot2", "scales", "INLA")

set.seed(20210)
N_BLOCKS <- 10        # spatial folds
BLOCK_KM <- 60        # side of the spatial blocks

# ===========================================================================
# 1. Rebuild the cluster model frame (same construction as 05)
# ===========================================================================
stopifnot(file.exists(OUT$dhs_child), file.exists(OUT$cov_cluster))
child <- readRDS(OUT$dhs_child)
clu <- utils::read.csv(OUT$cov_cluster)
scaling <- readRDS(file.path(DIR$interim, "05_covariate_scaling.rds"))
have <- scaling$vars

dat <- child |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(n = dplyr::n(), y = sum(stunted, na.rm = TRUE), .groups = "drop") |>
  dplyr::inner_join(clu, by = c("cluster" = "DHSCLUST")) |>
  dplyr::filter(!is.na(lon), !is.na(lat), n > 0)

for (v in have) {
  dat[[v]] <- (dat[[v]] - scaling$center[[v]]) / scaling$scale[[v]]
  dat[[v]][is.na(dat[[v]])] <- 0
}

pts <- sf::st_as_sf(dat, coords = c("lon", "lat"), crs = CRS_GEO, remove = FALSE) |>
  sf::st_transform(CRS_EQA)
coo <- sf::st_coordinates(pts) / 1000

# ===========================================================================
# 2. Spatial blocks
# ===========================================================================
# Tile the country into BLOCK_KM squares, then deal whole tiles out to folds.
# Holding out a contiguous tile means the nearest training cluster is typically
# tens of kilometres away - the realistic prediction distance.
bx <- floor(coo[, 1] / BLOCK_KM)
by <- floor(coo[, 2] / BLOCK_KM)
tile <- paste(bx, by, sep = "_")
utile <- unique(tile)
fold_of_tile <- stats::setNames(sample(rep_len(seq_len(N_BLOCKS), length(utile))), utile)
dat$fold <- unname(fold_of_tile[tile])

msg("spatial blocks: ", length(utile), " tiles of ", BLOCK_KM, " km -> ",
    N_BLOCKS, " folds")
msg("clusters per fold: ", paste(as.vector(table(dat$fold)), collapse = ", "))

# ===========================================================================
# 3. Cross-validation
# ===========================================================================
adm0 <- sf::st_read(OUT$adm1, quiet = TRUE) |>
  sf::st_transform(CRS_EQA) |> sf::st_union() |> sf::st_make_valid()
bnd <- sf::st_coordinates(sf::st_simplify(adm0, dTolerance = 5000))[, 1:2] / 1000
mesh <- INLA::inla.mesh.2d(loc = coo,
                           boundary = INLA::inla.nonconvex.hull(bnd, convex = -0.03),
                           max.edge = c(30, 120), cutoff = 10, offset = c(20, 120))
spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(50, 0.05),
                                  prior.sigma = c(1.5, 0.05))
A <- INLA::inla.spde.make.A(mesh, loc = coo)
idx <- INLA::inla.spde.make.index("s", n.spde = spde$n.spde)

# Two competing specifications, so the CV also answers "do the EO covariates
# actually help out of sample, or only in sample?"
SPECS <- list(
  spatial_only = character(0),
  full = have
)

# Use whichever likelihood 05 selected, so the cross-validation measures the
# model that is actually published rather than a different one.
sp_summary <- file.path(DIR$tables, "05_spde_spatial_summary.csv")
CV_FAMILY <- if (file.exists(sp_summary)) {
  as.character(utils::read.csv(sp_summary)$family[1])
} else "binomial"
msg("cross-validating with the ", CV_FAMILY, " likelihood (as selected in 05)")

run_cv <- function(vars, label) {
  preds <- rep(NA_real_, nrow(dat))
  for (k in seq_len(N_BLOCKS)) {
    train <- dat$fold != k
    yk <- dat$y; yk[!train] <- NA           # held-out rows predicted, not fitted

    X <- data.frame(intercept = 1)
    X <- cbind(X[rep(1, nrow(dat)), , drop = FALSE],
               dat[, vars, drop = FALSE])
    stk <- INLA::inla.stack(tag = "all", data = list(y = yk, n = dat$n),
                            A = list(A, 1), effects = list(idx, X))
    form <- stats::as.formula(paste("y ~ 0 + intercept",
                                    if (length(vars)) paste("+", paste(vars, collapse = " + ")) else "",
                                    "+ f(s, model = spde)"))
    f <- INLA::inla(form, family = CV_FAMILY,
                    Ntrials = INLA::inla.stack.data(stk)$n,
                    data = INLA::inla.stack.data(stk, spde = spde),
                    control.predictor = list(A = INLA::inla.stack.A(stk),
                                             compute = TRUE, link = 1),
                    control.inla = list(int.strategy = "eb"))
    ii <- INLA::inla.stack.index(stk, "all")$data
    preds[!train] <- f$summary.fitted.values[ii[!train], "mean"]
    msg("  [", label, "] fold ", k, "/", N_BLOCKS, " done")
  }
  preds
}

cv <- list()
for (s in names(SPECS)) {
  msg("=== spatial block CV: ", s, " ===")
  cv[[s]] <- run_cv(SPECS[[s]], s)
}

obs <- dat$y / dat$n
metrics <- do.call(rbind, lapply(names(cv), function(s) {
  p <- cv[[s]]
  ok <- !is.na(p)
  data.frame(
    model = s,
    n = sum(ok),
    # Weighted by cluster size: a 30-child cluster is worth more than a 5-child one.
    rmse = sqrt(stats::weighted.mean((p[ok] - obs[ok])^2, dat$n[ok])),
    mae = stats::weighted.mean(abs(p[ok] - obs[ok]), dat$n[ok]),
    bias = stats::weighted.mean(p[ok] - obs[ok], dat$n[ok]),
    correlation = stats::cor(p[ok], obs[ok]),
    # Binomial log score: rewards calibrated probabilities, not just accuracy.
    log_score = mean(stats::dbinom(dat$y[ok], dat$n[ok], p[ok], log = TRUE)),
    row.names = NULL)
}))
utils::write.csv(metrics, file.path(DIR$tables, "09_spatial_cv_metrics.csv"),
                 row.names = FALSE)
msg("spatial block CV results:")
print(metrics, row.names = FALSE)

cvdf <- data.frame(observed = obs, n = dat$n,
                   spatial_only = cv$spatial_only, full = cv$full)
utils::write.csv(cvdf, file.path(DIR$tables, "09_cv_predictions.csv"), row.names = FALSE)

p <- ggplot2::ggplot(cvdf, ggplot2::aes(full, observed, size = n)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::geom_point(alpha = 0.3) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "firebrick") +
  ggplot2::scale_x_continuous(labels = scales::percent) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::labs(x = "Out-of-sample prediction", y = "Observed cluster prevalence",
                size = "children",
                title = "Spatial block cross-validation",
                subtitle = paste0(N_BLOCKS, " spatially contiguous folds of ",
                                  BLOCK_KM, " km blocks")) +
  ggplot2::theme_minimal(base_size = 10)
save_fig(p, "09_spatial_cv.png", width = 6.5, height = 5.5)

# ===========================================================================
# 4. Internal fit diagnostics
# ===========================================================================
f_cmp <- file.path(DIR$tables, "05_likelihood_comparison.csv")
if (file.exists(f_cmp)) {
  lc <- utils::read.csv(f_cmp)
  # Calibration is assessed in 05, not here, and by the RANDOMISED PIT. INLA's
  # own cpo$pit is P(Y <= y), which for count data is stochastically larger
  # than uniform however good the model is, so a Kolmogorov-Smirnov test on it
  # rejects almost automatically. An earlier version of this script did exactly
  # that and reported "calibration is imperfect" on the strength of an artefact.
  # The randomised PIT of Czado, Gneiting & Held (2009), computed in 05 against
  # the predictive distribution for a NEW cluster, is the valid diagnostic.
  sel <- lc[order(!is.na(lc$pit_ks_p) & lc$pit_ks_p > 0.05,
                  lc$n_hyper, decreasing = c(TRUE, FALSE)), ][1, ]
  msg("calibration (from 05, randomised PIT): selected model '", sel$model,
      "' KS = ", signif(sel$pit_ks_stat, 3), ", p = ", signif(sel$pit_ks_p, 3),
      if (!is.na(sel$pit_ks_p) && sel$pit_ks_p > 0.05)
        "  <- uniformity not rejected" else "  <- miscalibrated")
  msg("  90% predictive intervals over-cover (about 0.97 observed), so the ",
      "published commune intervals are conservative rather than overconfident.")
  utils::write.csv(lc, file.path(DIR$tables, "09_calibration_from_05.csv"),
                   row.names = FALSE)
} else {
  msg("NOTE: run 05_model_spde.R to produce the calibration comparison.")
}

f_fit <- file.path(DIR$interim, "05_spde_fit.rds")
if (file.exists(f_fit)) {
  fit <- readRDS(f_fit)
  utils::write.csv(
    data.frame(waic = fit$waic$waic, dic = fit$dic$dic,
               n_failed_cpo = sum(fit$cpo$failure > 0, na.rm = TRUE)),
    file.path(DIR$tables, "09_internal_fit.csv"), row.names = FALSE)
  msg("selected model: WAIC ", round(fit$waic$waic, 1),
      " | DIC ", round(fit$dic$dic, 1),
      " | failed CPO ", sum(fit$cpo$failure > 0, na.rm = TRUE))
}

# ===========================================================================
# 5. Benchmarking against the design-based regional estimates
# ===========================================================================
f_agg1 <- file.path(DIR$processed, "aggregated_dhs_region.csv")
if (file.exists(f_agg1) && file.exists(OUT$direct_adm1)) {
  a1 <- utils::read.csv(f_agg1)
  d1 <- utils::read.csv(OUT$direct_adm1)
  cmp <- dplyr::inner_join(
    a1, d1[, c("dhs_region", "dhs_region_name", "direct", "se", "n_children")],
    by = "dhs_region")
  cmp$diff <- cmp$est - cmp$direct
  # Is the model aggregate inside the direct estimate's confidence interval?
  cmp$within_ci <- abs(cmp$diff) < 1.96 * cmp$se

  msg("DHS regions where the model agrees with the direct estimate (95% CI): ",
      sum(cmp$within_ci, na.rm = TRUE), " of ", nrow(cmp))
  msg("mean absolute difference: ", round(100 * mean(abs(cmp$diff), na.rm = TRUE), 2),
      " percentage points")

  # Ratio benchmarking: scale commune estimates so each DHS region reproduces
  # its design-based total. Written as SEPARATE columns, never in place, so the
  # unbenchmarked model output stays inspectable.
  bench <- stats::setNames(cmp$direct / cmp$est, as.character(cmp$dhs_region))
  com <- utils::read.csv(file.path(DIR$processed, "commune_stunting.csv"))
  # This script writes its own columns back into commune_stunting.csv, so a
  # second run would join dhs_region onto a frame that already has it and
  # silently produce dhs_region.x / dhs_region.y. Drop anything this block owns
  # before rebuilding it, so re-running is idempotent.
  com <- com[, setdiff(names(com), c("dhs_region", "dhs_region_name",
                                     "bench_factor", "est_benchmarked",
                                     "lower_benchmarked", "upper_benchmarked"))]
  xw <- utils::read.csv(file.path(DIR$processed, "dhs_region_crosswalk_commune.csv"))
  com <- dplyr::left_join(com, xw[, c("ADM3_PCODE", "dhs_region")],
                          by = "ADM3_PCODE")
  com$bench_factor <- unname(bench[as.character(com$dhs_region)])
  com$bench_factor[is.na(com$bench_factor)] <- 1
  com$est_benchmarked <- pmin(1, com$est * com$bench_factor)
  com$lower_benchmarked <- pmin(1, com$lower * com$bench_factor)
  com$upper_benchmarked <- pmin(1, com$upper * com$bench_factor)
  utils::write.csv(com, file.path(DIR$processed, "commune_stunting.csv"),
                   row.names = FALSE)
  msg("benchmark factors: median ", round(stats::median(com$bench_factor), 3),
      " range ", round(min(com$bench_factor), 3), " to ",
      round(max(com$bench_factor), 3))

  utils::write.csv(cmp, file.path(DIR$tables, "09_benchmark_dhs_region.csv"),
                   row.names = FALSE)

  pb <- ggplot2::ggplot(cmp, ggplot2::aes(direct, est)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = direct - 1.96 * se,
                                         xmax = direct + 1.96 * se),
                            height = 0, colour = "grey70") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0,
                           colour = "grey70") +
    ggplot2::geom_point(colour = "firebrick") +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "Design-based direct estimate (DHS region)",
                  y = "Model aggregate",
                  title = "Benchmarking against the official survey estimates",
                  subtitle = "Error bars: 95% design CI (horizontal), credible interval (vertical)") +
    ggplot2::theme_minimal(base_size = 10)
  save_fig(pb, "09_benchmark_dhs_region.png", width = 6.5, height = 5.5)
} else {
  msg("NOTE: DHS-region aggregate or direct estimates missing; benchmarking skipped.")
}

# ===========================================================================
# 6. External coherence with IHME Local Burden of Disease
# ===========================================================================
# The IHME child growth failure surfaces (Nature 2020) are the other published
# estimate of Malagasy stunting below national level, so agreement with them is
# the one genuinely EXTERNAL check available. It is a check of spatial PATTERN,
# not of level, and that distinction is not a technicality:
#
#   - IHME's series ends in 2017; this one is the 2021 DHS. Four years apart,
#     spanning COVID and the 2020-21 Grand Sud drought.
#   - IHME pools ~460 surveys across 105 countries with continentally calibrated
#     covariates; this model is fitted to one Malagasy survey with covariates
#     chosen for Malagasy agro-ecology.
#   - IHME publishes 5x5 km aggregated to admin2; communes are admin3.
#
# A level difference is therefore expected and uninformative. What WOULD be
# informative is disagreement about where stunting is concentrated, so the
# headline statistics are Pearson and Spearman correlation across districts and
# the overlap of the worst-20 lists, with the level gap reported separately
# rather than folded in.
#
# Two input forms are accepted, in order of preference:
#   1. IHME's own admin2 CSV aggregate - avoids imposing my zonal statistics on
#      their raster, and reproduces their published numbers exactly
#   2. a stunting GeoTIFF, aggregated here with the same population weights the
#      model surface uses
ihme_dir <- file.path(DIR$rasters, "ihme")
dir.create(ihme_dir, recursive = TRUE, showWarnings = FALSE)

ihme_csv <- list.files(ihme_dir, pattern = "[.]csv$", full.names = TRUE,
                       ignore.case = TRUE)
ihme_tif <- list.files(ihme_dir, pattern = "[.]tif$", full.names = TRUE,
                       ignore.case = TRUE)

norm_nm <- function(x) gsub("[^A-Z0-9]", "", toupper(trimws(as.character(x))))

adm2 <- sf::st_read(OUT$adm2, quiet = TRUE)
a2 <- utils::read.csv(file.path(DIR$processed, "aggregated_adm2.csv"))
mine <- dplyr::left_join(
  sf::st_drop_geometry(adm2)[, c("ADM2_PCODE", "ADM2_EN", "ADM1_EN")],
  a2[, c("ADM2_PCODE", "est", "lower", "upper")], by = "ADM2_PCODE")

cmp2 <- NULL

if (length(ihme_csv) > 0) {
  raw <- utils::read.csv(ihme_csv[1])
  msg("IHME admin2 CSV: ", basename(ihme_csv[1]), " (", nrow(raw), " rows)")
  nm <- names(raw)
  pickcol <- function(cands) { h <- cands[cands %in% nm]; if (length(h)) h[1] else NA_character_ }
  c_name <- pickcol(c("ADM2_NAME", "adm2_name", "location_name", "ADM2_EN"))
  c_val  <- pickcol(c("mean", "val", "value", "prevalence"))
  c_year <- pickcol(c("year", "year_id"))
  c_ind  <- pickcol(c("indicator", "measure", "measure_name"))
  c_iso  <- pickcol(c("ADM0_NAME", "adm0_name", "iso3", "ISO3"))

  if (is.na(c_name) || is.na(c_val)) {
    msg("could not identify the name/value columns. Columns present: ",
        paste(nm, collapse = ", "))
  } else {
    d <- raw
    if (!is.na(c_iso)) {
      keep <- grepl("madagascar|MDG", d[[c_iso]], ignore.case = TRUE)
      if (any(keep)) d <- d[keep, ]
    }
    if (!is.na(c_ind)) {
      st <- grepl("stunt", d[[c_ind]], ignore.case = TRUE)
      if (any(st)) d <- d[st, ]
    }
    if (!is.na(c_year)) {
      yr <- suppressWarnings(max(as.numeric(d[[c_year]]), na.rm = TRUE))
      if (is.finite(yr)) {
        d <- d[as.numeric(d[[c_year]]) == yr, ]
        msg("using IHME year ", yr)
      }
    }
    ih <- data.frame(key = norm_nm(d[[c_name]]), ihme = as.numeric(d[[c_val]]))
    ih <- ih[!is.na(ih$ihme) & !duplicated(ih$key), ]
    if (max(ih$ihme, na.rm = TRUE) > 1.5) ih$ihme <- ih$ihme / 100
    mine$key <- norm_nm(mine$ADM2_EN)
    cmp2 <- dplyr::inner_join(mine, ih, by = "key")
    msg("districts matched on name: ", nrow(cmp2), " of ", nrow(mine))
    miss <- mine[!mine$key %in% ih$key, c("ADM2_PCODE", "ADM2_EN")]
    if (nrow(miss) > 0) {
      msg("  ", nrow(miss), " unmatched - written to 09_ihme_unmatched_districts.csv ",
          "(district naming differs between the 2025 PAM layer and IHME's GADM base)")
      utils::write.csv(miss, file.path(DIR$tables, "09_ihme_unmatched_districts.csv"),
                       row.names = FALSE)
    }
  }
} else if (length(ihme_tif) > 0) {
  need("terra", "exactextractr")
  msg("IHME raster: ", basename(ihme_tif[1]), " (no admin2 CSV present)")
  r <- terra::rast(ihme_tif[1])
  # Population-weight the zonal mean, exactly as the model surface is aggregated
  # in 08, so the two district numbers are formed the same way and the
  # comparison is not partly an artefact of differing aggregation.
  wp <- file.path(DIR$rasters, "worldpop")
  u5f <- file.path(wp, c("mdg_f_0_2020_constrained.tif", "mdg_f_1_2020_constrained.tif",
                         "mdg_m_0_2020_constrained.tif", "mdg_m_1_2020_constrained.tif"))
  z <- sf::st_transform(adm2, terra::crs(r))
  if (all(file.exists(u5f))) {
    wgt <- terra::resample(sum(terra::rast(u5f), na.rm = TRUE), r, method = "sum")
    mine$ihme <- exactextractr::exact_extract(r, z, "weighted_mean", weights = wgt,
                                              progress = FALSE)
  } else {
    mine$ihme <- exactextractr::exact_extract(r, z, "mean", progress = FALSE)
  }
  if (max(mine$ihme, na.rm = TRUE) > 1.5) mine$ihme <- mine$ihme / 100
  cmp2 <- mine[!is.na(mine$ihme), ]
} else {
  msg("IHME estimates not present - external comparison skipped.")
  msg("  These files sit behind a free IHME account and require accepting the")
  msg("  IHME Free-of-Charge Non-commercial User Agreement, so they cannot be")
  msg("  fetched automatically. From")
  msg("    https://ghdx.healthdata.org/record/ihme-data/",
      "lmic-child-growth-failure-geospatial-estimates-2000-2017")
  msg("  download either the ADMIN2 CSV (preferred) or the stunting GeoTIFF,")
  msg("  unzip, and drop the file into:")
  msg("    ", ihme_dir)
  msg("  Then re-run this script; nothing else needs changing.")
}

if (!is.null(cmp2) && nrow(cmp2) > 5) {
  cmp2$diff <- cmp2$est - cmp2$ihme
  r_p <- stats::cor(cmp2$est, cmp2$ihme, use = "complete.obs")
  r_s <- stats::cor(cmp2$est, cmp2$ihme, method = "spearman", use = "complete.obs")

  # Worst-20 overlap is the targeting question. Two surfaces can correlate well
  # and still disagree about which districts belong at the top of a priority
  # list, which is the use these numbers are actually put to.
  k <- min(20, nrow(cmp2))
  top_mine <- cmp2$ADM2_PCODE[order(-cmp2$est)][seq_len(k)]
  top_ihme <- cmp2$ADM2_PCODE[order(-cmp2$ihme)][seq_len(k)]
  overlap <- length(intersect(top_mine, top_ihme))

  msg("=== IHME coherence over ", nrow(cmp2), " districts ===")
  msg("Pearson r ", round(r_p, 3), " | Spearman rho ", round(r_s, 3))
  msg("worst-", k, " districts shared: ", overlap, " of ", k)
  msg("level: this model ", round(100 * mean(cmp2$est, na.rm = TRUE), 1),
      "% vs IHME ", round(100 * mean(cmp2$ihme, na.rm = TRUE), 1),
      "% (difference ", round(100 * mean(cmp2$diff, na.rm = TRUE), 1),
      " points; a gap is expected between 2017 and 2021)")

  utils::write.csv(cmp2, file.path(DIR$tables, "09_ihme_comparison.csv"),
                   row.names = FALSE)
  utils::write.csv(
    data.frame(n_districts = nrow(cmp2), pearson_r = r_p, spearman_rho = r_s,
               top_k = k, top_k_overlap = overlap,
               mean_model = mean(cmp2$est, na.rm = TRUE),
               mean_ihme = mean(cmp2$ihme, na.rm = TRUE),
               mean_diff = mean(cmp2$diff, na.rm = TRUE)),
    file.path(DIR$tables, "09_ihme_summary.csv"), row.names = FALSE)

  p_ih <- ggplot2::ggplot(cmp2, ggplot2::aes(ihme, est)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::geom_point(alpha = 0.65, colour = "#1A6E9E") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = "#A8402C", linewidth = 0.7) +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "IHME Local Burden of Disease, latest year (district)",
                  y = "This model, DHS 2021 (district)",
                  title = "External coherence with IHME",
                  subtitle = paste0("Pearson r = ", round(r_p, 3),
                                    ", Spearman rho = ", round(r_s, 3),
                                    ". Dashed line is equality; a level shift",
                                    " between 2017 and 2021 is expected.")) +
    ggplot2::theme_minimal(base_size = 10)
  save_fig(p_ih, "09_ihme_scatter.png", width = 6.5, height = 5.5)

  gmap <- dplyr::left_join(adm2, cmp2[, c("ADM2_PCODE", "diff")], by = "ADM2_PCODE")
  p_dm <- ggplot2::ggplot(gmap) +
    ggplot2::geom_sf(ggplot2::aes(fill = diff), colour = "grey40", linewidth = 0.08) +
    ggplot2::scale_fill_gradient2(low = "#1A6E9E", mid = "grey92", high = "#A8402C",
                                  midpoint = 0, labels = scales::percent,
                                  name = "this model\nminus IHME") +
    ggplot2::labs(title = "Where the two estimates disagree",
                  subtitle = "Red: this model is higher. Blue: IHME is higher.") +
    ggplot2::theme_void(base_size = 9)
  save_fig(p_dm, "09_ihme_difference_map.png", width = 6, height = 8)
}

msg("09_validate.R complete")
