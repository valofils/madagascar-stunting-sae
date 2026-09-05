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
    f <- INLA::inla(form, family = "binomial",
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
f_fit <- file.path(DIR$interim, "05_spde_fit.rds")
if (file.exists(f_fit)) {
  fit <- readRDS(f_fit)
  # PIT should be uniform if the predictive distribution is calibrated.
  pit <- fit$cpo$pit[seq_len(nrow(dat))]
  pit <- pit[is.finite(pit)]
  ks <- stats::ks.test(pit, "punif")
  msg("PIT uniformity (Kolmogorov-Smirnov) p = ", signif(ks$p.value, 3),
      if (ks$p.value < 0.05) "  <- calibration is imperfect" else "  <- calibrated")
  utils::write.csv(
    data.frame(waic = fit$waic$waic, dic = fit$dic$dic,
               n_failed_cpo = sum(fit$cpo$failure > 0, na.rm = TRUE),
               pit_ks_p = ks$p.value),
    file.path(DIR$tables, "09_internal_fit.csv"), row.names = FALSE)

  ph <- ggplot2::ggplot(data.frame(pit = pit), ggplot2::aes(pit)) +
    ggplot2::geom_histogram(bins = 20, fill = "steelblue", colour = "white") +
    ggplot2::geom_hline(yintercept = length(pit) / 20, linetype = "dashed") +
    ggplot2::labs(x = "PIT", y = "clusters",
                  title = "Probability integral transform",
                  subtitle = "A calibrated model gives a flat histogram") +
    ggplot2::theme_minimal(base_size = 10)
  save_fig(ph, "09_pit_histogram.png", width = 6, height = 4)
}

# ===========================================================================
# 5. Benchmarking against the design-based regional estimates
# ===========================================================================
f_agg1 <- file.path(DIR$processed, "aggregated_adm1.csv")
if (file.exists(f_agg1) && file.exists(OUT$direct_adm1)) {
  a1 <- utils::read.csv(f_agg1)
  d1 <- utils::read.csv(OUT$direct_adm1)
  cmp <- dplyr::inner_join(a1, d1[, c("ADM1_PCODE", "direct", "se", "n_children")],
                           by = "ADM1_PCODE")
  cmp$diff <- cmp$est - cmp$direct
  # Is the model aggregate inside the direct estimate's confidence interval?
  cmp$within_ci <- abs(cmp$diff) < 1.96 * cmp$se

  msg("regions where the model agrees with the direct estimate (95% CI): ",
      sum(cmp$within_ci, na.rm = TRUE), " of ", nrow(cmp))
  msg("mean absolute difference: ", round(100 * mean(abs(cmp$diff), na.rm = TRUE), 2),
      " percentage points")

  # Ratio benchmarking: scale commune estimates so each region reproduces its
  # design-based total. Applied as a separate output column, never in place,
  # so the unbenchmarked model results stay inspectable.
  bench_factor <- stats::setNames(cmp$direct / cmp$est, cmp$ADM1_PCODE)
  com <- utils::read.csv(file.path(DIR$processed, "commune_stunting.csv"))
  com$bench_factor <- unname(bench_factor[com$ADM1_PCODE])
  com$bench_factor[is.na(com$bench_factor)] <- 1
  com$est_benchmarked <- pmin(1, com$est * com$bench_factor)
  com$lower_benchmarked <- pmin(1, com$lower * com$bench_factor)
  com$upper_benchmarked <- pmin(1, com$upper * com$bench_factor)
  utils::write.csv(com, file.path(DIR$processed, "commune_stunting.csv"),
                   row.names = FALSE)

  utils::write.csv(cmp, file.path(DIR$tables, "09_benchmark_adm1.csv"),
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
                  title = "Benchmarking against the official survey estimates") +
    ggplot2::theme_minimal(base_size = 10)
  save_fig(pb, "09_benchmark_adm1.png", width = 6.5, height = 5.5)
}

# ===========================================================================
# 6. External coherence with IHME
# ===========================================================================
# Optional: the IHME Local Burden of Disease stunting surface (2000-2017) is
# not redistributable here, so this section runs only if a local copy exists.
ihme <- file.path(DIR$rasters, "ihme", "IHME_stunting_prev_2017.tif")
if (file.exists(ihme)) {
  need("terra", "exactextractr")
  adm2 <- sf::st_read(OUT$adm2, quiet = TRUE)
  r <- terra::rast(ihme)
  adm2$ihme <- exactextractr::exact_extract(r, adm2, "mean", progress = FALSE)
  a2 <- utils::read.csv(file.path(DIR$processed, "aggregated_adm2.csv"))
  cmp2 <- dplyr::inner_join(sf::st_drop_geometry(adm2)[, c("ADM2_PCODE", "ADM2_EN", "ihme")],
                            a2, by = "ADM2_PCODE")
  msg("correlation with IHME at district level: ",
      round(stats::cor(cmp2$est, cmp2$ihme, use = "complete.obs"), 3),
      "  (IHME is a 2017 surface; a difference in level is expected)")
  utils::write.csv(cmp2, file.path(DIR$tables, "09_ihme_comparison.csv"),
                   row.names = FALSE)
} else {
  msg("IHME raster not present - external comparison skipped.")
  msg("  Place IHME_stunting_prev_2017.tif in ", file.path(DIR$rasters, "ihme"),
      " to enable it.")
}

msg("09_validate.R complete")
