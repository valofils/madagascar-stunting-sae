# ---------------------------------------------------------------------------
# 05_model_spde.R
#
# Continuous-surface model: binomial geostatistical regression on DHS cluster
# locations, with a Matern spatial field represented through the SPDE
# approximation of Lindgren, Rue & Lindstrom (2011) and fitted by INLA.
#
# Model:
#   y_c ~ Binomial or BetaBinomial(n_c, p_c)     y_c = stunted children in cluster c
#   logit(p_c) = alpha + x_c' beta + S(s_c) [+ e_c]
#   S(.)  ~ GP with Matern covariance (SPDE), nu = 1
#   e_c   ~ N(0, sigma_e^2)                      optional cluster-level nugget
#
# Children within a DHS cluster share a village, a water source and a food
# system, so their outcomes are correlated and a plain binomial understates the
# variance. Three specifications - binomial with a nugget, beta-binomial, and
# beta-binomial with a nugget - are fitted and compared on WAIC, CPO, Pearson
# dispersion and PIT, and the surface is built from whichever wins. See
# outputs/tables/05_likelihood_comparison.csv.
#
# The fitted field is then projected onto a prediction grid covering
# Madagascar, giving the 5 km surface that 08_aggregate.R averages up to
# communes with under-5 population weights.
#
# Inputs : data/interim/dhs_child_haz.rds        (from 03)
#          data/processed/cluster_covariates.csv (from 02)
#          data/processed/commune_covariates.csv (from 02)
# Outputs: data/processed/spde_grid_predictions.rds
#          data/processed/spde_commune.csv
#          outputs/figures/05_*.png, outputs/tables/05_*.csv
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("dplyr", "sf", "ggplot2", "scales", "INLA", "terra")

set.seed(20210)
INLA::inla.setOption(num.threads = "2:1")

# ===========================================================================
# 1. Cluster-level binomial data
# ===========================================================================
stopifnot(file.exists(OUT$dhs_child), file.exists(OUT$cov_cluster))
child <- readRDS(OUT$dhs_child)
clu <- utils::read.csv(OUT$cov_cluster)

# Aggregate children to clusters. The design weights are carried through as a
# cluster-level mean weight; the binomial likelihood itself is unweighted,
# and the design is re-introduced at the benchmarking step in 09.
dat <- child |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(n = dplyr::n(),
                   y = sum(stunted, na.rm = TRUE),
                   mean_haz = mean(haz, na.rm = TRUE),
                   wt = mean(wt, na.rm = TRUE), .groups = "drop") |>
  dplyr::inner_join(clu, by = c("cluster" = "DHSCLUST")) |>
  dplyr::filter(!is.na(lon), !is.na(lat), n > 0)

msg("clusters in the model: ", nrow(dat), " | children: ", sum(dat$n),
    " | crude stunting: ", round(100 * sum(dat$y) / sum(dat$n), 1), "%")

# ===========================================================================
# 2. Covariates
# ===========================================================================
COVARS <- c("elevation", "ruggedness", "temp_min_cold", "temp_seasonality",
            "precip_annual", "precip_seasonality", "frac_crop", "frac_built",
            "frac_water_perm", "cattle_density", "travel_time",
            "travel_time_healthcare", "nightlights", "urban")
have <- intersect(COVARS, names(dat))
have <- have[vapply(have, function(v) {
  ok <- sum(!is.na(dat[[v]])) > 0.9 * nrow(dat) && stats::sd(dat[[v]], na.rm = TRUE) > 0
  if (!ok) msg("dropping covariate (missing or constant): ", v)
  ok
}, logical(1))]
msg("covariates in the model: ", paste(have, collapse = ", "))

# Standardise using the CLUSTER moments, and reuse exactly those moments on the
# prediction grid later - otherwise the fitted betas would be applied to
# differently-scaled inputs and the surface would be silently wrong.
ctr <- vapply(dat[have], mean, numeric(1), na.rm = TRUE)
scl <- vapply(dat[have], stats::sd, numeric(1), na.rm = TRUE)
scale_with <- function(df, vars, ctr, scl) {
  for (v in vars) df[[v]] <- (df[[v]] - ctr[[v]]) / scl[[v]]
  df
}
dat <- scale_with(dat, have, ctr, scl)
# A cluster with one missing covariate would otherwise drop out entirely.
for (v in have) dat[[v]][is.na(dat[[v]])] <- 0
saveRDS(list(center = ctr, scale = scl, vars = have),
        file.path(DIR$interim, "05_covariate_scaling.rds"))

# ===========================================================================
# 3. Mesh
# ===========================================================================
# Distances in kilometres on the equal-area projection: interpretable range
# parameters and a numerically better-conditioned mesh than lon/lat degrees.
pts <- sf::st_as_sf(dat, coords = c("lon", "lat"), crs = CRS_GEO) |>
  sf::st_transform(CRS_EQA)
coo <- sf::st_coordinates(pts) / 1000

adm0 <- sf::st_read(OUT$adm1, quiet = TRUE) |>
  sf::st_transform(CRS_EQA) |> sf::st_union() |> sf::st_make_valid()
bnd <- sf::st_coordinates(sf::st_simplify(adm0, dTolerance = 5000))[, 1:2] / 1000

mesh <- INLA::inla.mesh.2d(
  loc = coo,
  boundary = INLA::inla.nonconvex.hull(bnd, convex = -0.03),
  max.edge = c(25, 100),   # km: fine inside the country, coarse in the buffer
  cutoff = 8,              # merge cluster locations closer than 8 km
  offset = c(20, 120)
)
msg("SPDE mesh nodes: ", mesh$n)

# ---- PC priors on the Matern field ---------------------------------------
# P(range < 50 km) = 0.05  : the field should not be shorter-ranged than the
#   ~50 km scale at which Madagascar's agro-ecological zones vary.
# P(sigma > 1.5) = 0.05    : on the logit scale, sd 1.5 is already a very large
#   residual spatial effect, so this is a weakly informative upper bound.
spde <- INLA::inla.spde2.pcmatern(mesh = mesh,
                                  prior.range = c(50, 0.05),
                                  prior.sigma = c(1.5, 0.05))

# ===========================================================================
# 4. Prediction grid
# ===========================================================================
# Prediction grid resolution. 5 km rather than 1 km, chosen deliberately:
#   - 1 km over Madagascar is 632,512 cells. Carrying those as prediction rows
#     inside the INLA stack (with config = TRUE, needed for posterior sampling)
#     exhausts memory and the inla binary crashes.
#   - Nothing in the model actually resolves 1 km. Land cover enters at commune
#     resolution, the coarsest covariate (cattle) is 10 km, the DHS coordinates
#     are displaced up to 5 km, and the Matern range is tens of kilometres.
#   - 5 km still gives ~14 cells inside an average commune (mean area 350 km2),
#     which is ample for population-weighted aggregation, and matches the
#     resolution IHME publishes at.
# Prediction is done by posterior sampling onto the grid AFTER fitting rather
# than by putting the grid in the stack, which is both the standard SPDE
# workflow and what makes the memory manageable.
GRID_KM <- 5
grid_file <- file.path(DIR$interim, "05_prediction_grid.rds")

build_grid <- function() {
  r <- terra::rast(terra::vect(sf::st_transform(adm0, CRS_GEO)),
                   resolution = GRID_KM / 111)     # approx. degrees
  terra::values(r) <- 1
  r <- terra::mask(r, terra::vect(sf::st_transform(adm0, CRS_GEO)))
  g <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)[, c("x", "y")]
  names(g) <- c("lon", "lat")
  msg("prediction grid cells: ", nrow(g))
  g
}
grid <- if (file.exists(grid_file)) readRDS(grid_file) else {
  g <- build_grid(); saveRDS(g, grid_file); g
}

# Covariates on the grid: re-extract the same rasters used in 02, at points.
grid_cov_file <- file.path(DIR$interim, "05_grid_covariates.rds")
if (file.exists(grid_cov_file)) {
  grid <- readRDS(grid_cov_file)
} else {
  msg("extracting covariates on the prediction grid")
  source("R/02_covariates_helpers.R", local = TRUE)   # raster getters, no side effects
  gp <- terra::vect(as.matrix(grid[, c("lon", "lat")]), type = "points",
                    crs = "EPSG:4326")
  for (nm in names(RASTER_GETTERS)) {
    r <- tryCatch(RASTER_GETTERS[[nm]](), error = function(e) NULL)
    if (is.null(r)) { msg("grid: skipped ", nm); next }
    v <- terra::extract(r, terra::project(gp, terra::crs(r)), ID = FALSE)
    grid <- cbind(grid, v)
  }
  # Population density from the commune layer (grid cells inherit their commune)
  saveRDS(grid, grid_cov_file)
}

# Attach each grid cell to its commune, needed for aggregation in 08.
if (!"ADM3_PCODE" %in% names(grid)) {
  adm3 <- sf::st_read(OUT$adm3, quiet = TRUE)
  gsf <- sf::st_as_sf(grid, coords = c("lon", "lat"), crs = CRS_GEO, remove = FALSE)
  j <- sf::st_join(gsf["lon"], adm3["ADM3_PCODE"], join = sf::st_intersects)
  grid$ADM3_PCODE <- j$ADM3_PCODE
  # Displacement of the coastline between the polygon layer and the raster grid
  # leaves a fringe of cells outside every commune; snap them so no populated
  # cell is dropped at the aggregation step.
  miss <- which(is.na(grid$ADM3_PCODE))
  if (length(miss) > 0) {
    nearest <- sf::st_nearest_feature(gsf[miss, ], adm3)
    grid$ADM3_PCODE[miss] <- adm3$ADM3_PCODE[nearest]
    msg("snapped ", length(miss), " grid cells to the nearest commune")
  }
  saveRDS(grid, grid_cov_file)
}

# Land cover enters the surface at COMMUNE resolution, not 1 km. The WorldCover
# class fractions are computed by zonal statistics over a 10 m grid, which is
# affordable for 1,701 commune polygons and for 650 cluster buffers but not for
# ~600,000 grid cells. Each cell therefore inherits its commune's fractions.
#
# The cost is explicit: within-commune variation in cropland, built-up area and
# water is not represented in the fixed effects, and is left to the SPDE field
# to absorb. Since the cluster-level covariates ARE measured over a 5 km buffer,
# the fitted coefficients apply to a finer scale than the prediction grid
# supplies, which will tend to attenuate their apparent effect on the surface.
# Dropping these covariates entirely would be worse - frac_crop carries H2 - but
# the surface should not be read as resolving land cover below commune level.
lc_cols <- grep("^frac_", have, value = TRUE)
if (length(lc_cols) > 0 && !all(lc_cols %in% names(grid))) {
  cov_com <- utils::read.csv(OUT$cov_commune)
  keep_lc <- intersect(lc_cols, names(cov_com))
  grid <- dplyr::left_join(grid, cov_com[, c("ADM3_PCODE", keep_lc)],
                           by = "ADM3_PCODE")
  msg("land cover attached at commune resolution: ",
      paste(keep_lc, collapse = ", "))
}

grid_use <- grid
missing_cov <- setdiff(have, names(grid_use))
if (length(missing_cov) > 0) {
  msg("WARNING: covariates absent from the grid, dropped from prediction: ",
      paste(missing_cov, collapse = ", "))
  have <- setdiff(have, missing_cov)
}
grid_use <- scale_with(grid_use, have, ctr, scl)
for (v in have) grid_use[[v]][is.na(grid_use[[v]])] <- 0

gpts <- sf::st_as_sf(grid_use, coords = c("lon", "lat"), crs = CRS_GEO,
                     remove = FALSE) |> sf::st_transform(CRS_EQA)
coo_pred <- sf::st_coordinates(gpts) / 1000

# ===========================================================================
# 5. Fit (estimation locations only)
# ===========================================================================
A_est <- INLA::inla.spde.make.A(mesh, loc = coo)
A_pred <- INLA::inla.spde.make.A(mesh, loc = coo_pred)
idx <- INLA::inla.spde.make.index("s", n.spde = spde$n.spde)

fixed_df <- function(d) {
  x <- data.frame(intercept = rep(1, nrow(d)))
  cbind(x, as.data.frame(d)[, have, drop = FALSE])
}

stk_est <- INLA::inla.stack(
  tag = "est", data = list(y = dat$y, n = dat$n),
  A = list(A_est, 1),
  effects = list(idx, cbind(fixed_df(dat), nugget = seq_len(nrow(dat)))))

# ---- Candidate likelihoods -----------------------------------------------
# The binomial + cluster-nugget specification left the PIT histogram clearly
# non-uniform (KS p = 1.3e-9), which is the signature of extra-binomial
# variation between clusters: children within a DHS cluster share a village,
# a water source and a food system, so their outcomes are correlated and the
# binomial variance n*p*(1-p) is too small.
#
# Two ways to absorb that, and they are not the same thing:
#   - an iid nugget on the LOGIT scale, which lets each cluster's underlying
#     risk depart from the surface, or
#   - a beta-binomial likelihood, which models the intra-cluster correlation
#     directly: Var(y) = n p (1-p) [1 + (n-1) rho].
# The second is the better-matched device, because the overdispersion here is a
# property of how children cluster within a village, not of the risk surface.
# Carrying both is close to unidentifiable, so all three are fitted and compared
# rather than assumed.
SPECS <- list(
  binomial = list(
    family = "binomial",
    form = stats::as.formula(paste(
      "y ~ 0 + intercept +", paste(have, collapse = " + "),
      "+ f(s, model = spde)"))),
  binomial_nugget = list(
    family = "binomial",
    form = stats::as.formula(paste(
      "y ~ 0 + intercept +", paste(have, collapse = " + "),
      "+ f(s, model = spde)",
      "+ f(nugget, model = 'iid', hyper = list(prec = list(prior = 'pc.prec',",
      "                                                    param = c(1, 0.01))))"))),
  betabinomial = list(
    family = "betabinomial",
    form = stats::as.formula(paste(
      "y ~ 0 + intercept +", paste(have, collapse = " + "),
      "+ f(s, model = spde)"))),
  betabinomial_nugget = list(
    family = "betabinomial",
    form = stats::as.formula(paste(
      "y ~ 0 + intercept +", paste(have, collapse = " + "),
      "+ f(s, model = spde)",
      "+ f(nugget, model = 'iid', hyper = list(prec = list(prior = 'pc.prec',",
      "                                                    param = c(1, 0.01))))")))
)

stk_dat <- INLA::inla.stack.data(stk_est, spde = spde)
Ntr <- stk_dat$n
idx_est <- INLA::inla.stack.index(stk_est, "est")$data

fit_spec <- function(sp, label) {
  msg("fitting [", label, "] (", sp$family, ")")
  t0 <- Sys.time()
  f <- INLA::inla(
    sp$form, family = sp$family, Ntrials = Ntr, data = stk_dat,
    control.predictor = list(A = INLA::inla.stack.A(stk_est), compute = TRUE, link = 1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    control.inla = list(int.strategy = "eb"),
    verbose = FALSE)
  msg("  done in ", round(difftime(Sys.time(), t0, units = "secs")), "s")
  f
}

fits <- lapply(names(SPECS), function(k) fit_spec(SPECS[[k]], k))
names(fits) <- names(SPECS)

# ---- Calibration diagnostics ---------------------------------------------
# Pearson dispersion: sum of squared standardised residuals over the number of
# clusters. About 1 means the likelihood's variance matches the data; above 1
# means the model still understates cluster-level variability.
#
# The model-implied variance differs by family:
#   binomial      V = n p (1 - p)
#   betabinomial  V = n p (1 - p) [1 + (n - 1) rho]
dispersion_of <- function(f, family) {
  p <- f$summary.fitted.values[idx_est, "mean"]
  V <- dat$n * p * (1 - p)
  rho <- NA_real_
  if (family == "betabinomial") {
    hp <- f$summary.hyperpar
    r <- grep("rho|overdispersion", rownames(hp), ignore.case = TRUE)
    if (length(r) > 0) {
      rho <- hp[r[1], "mean"]
      V <- V * (1 + (dat$n - 1) * rho)
    }
  }
  list(dispersion = sum((dat$y - dat$n * p)^2 / V) / nrow(dat), rho = rho)
}

# ---- Calibration: RANDOMISED PIT ------------------------------------------
# INLA reports pit = P(Y <= y_obs). For CONTINUOUS data that is uniform under a
# correct model, but these are counts, and for discrete data P(Y <= y) is
# stochastically larger than uniform no matter how good the model is. Testing it
# with a Kolmogorov-Smirnov test is therefore invalid - which is exactly what
# R's "ties should not be present" warning was saying. An earlier version of
# this pipeline read that artefact as evidence of miscalibration; it is not.
#
# The correct device for discrete outcomes is the randomised PIT of Czado,
# Gneiting & Held (2009):
#
#     u_i = F_i(y_i - 1) + v_i * [ F_i(y_i) - F_i(y_i - 1) ],   v_i ~ U(0, 1)
#
# which IS uniform under a correctly specified model. F_i is estimated here by
# simulating from the posterior predictive distribution, so it also propagates
# parameter uncertainty rather than conditioning on point estimates.
N_PPC <- 300

# A further subtlety, and it changes the answer. The predictive distribution
# used for calibration must be the one for a NEW cluster at that location, not
# the one conditioned on the cluster's own fitted nugget. INLA's "Predictor"
# rows contain the latter: the nugget has been fitted to that very cluster's
# data, so replicating from it reproduces the observations almost too well and
# the check becomes in-sample. That is precisely how a too-wide predictive
# distribution can coexist with a Pearson dispersion below 1.
#
# So the linear predictor is rebuilt from its components -
#     eta = x' beta + A s   ( + a FRESH nugget draw where the spec has one )
# - which is exactly how the surface is used in 08 to predict unsampled
# communes. This makes the four specifications comparable on equal terms.
randomised_pit <- function(f, family, has_nugget) {
  smp <- INLA::inla.posterior.sample(N_PPC, f, seed = 20210)
  rn <- rownames(smp[[1]]$latent)
  s_rows <- grep("^s:", rn)
  fx_rows <- vapply(c("intercept", have),
                    function(v) which(rn == paste0(v, ":1")), integer(1))
  Xe <- as.matrix(cbind(intercept = 1, dat[, have, drop = FALSE]))
  n_i <- dat$n; y_i <- dat$y

  yrep <- vapply(smp, function(z) {
    lat <- z$latent[, 1]
    eta <- as.numeric(Xe %*% lat[fx_rows]) + as.numeric(A_est %*% lat[s_rows])
    hp <- z$hyperpar
    if (has_nugget) {
      pr <- grep("nugget", names(hp), ignore.case = TRUE)
      if (length(pr) > 0) {
        sd_e <- 1 / sqrt(as.numeric(hp[pr[1]]))
        eta <- eta + stats::rnorm(length(eta), 0, sd_e)
      }
    }
    p <- stats::plogis(eta)
    if (family == "betabinomial") {
      # INLA parameterises the beta-binomial by rho, the intra-cluster
      # correlation, with rho = 1 / (a + b + 1). Recover (a, b) from (p, rho).
      r <- grep("rho|overdispersion", names(hp), ignore.case = TRUE)
      rho <- if (length(r) > 0) as.numeric(hp[r[1]]) else 0
      rho <- min(max(rho, 1e-8), 1 - 1e-8)
      ab <- 1 / rho - 1
      stats::rbinom(length(p), n_i, stats::rbeta(length(p), p * ab, (1 - p) * ab))
    } else {
      stats::rbinom(length(p), n_i, p)
    }
  }, numeric(nrow(dat)))

  F_lt <- rowMeans(yrep < y_i)     # F(y-1)
  F_eq <- rowMeans(yrep == y_i)    # P(Y = y)
  u <- F_lt + stats::runif(length(y_i)) * F_eq
  kt <- suppressWarnings(stats::ks.test(u, "punif"))
  # Coverage of the 90% predictive interval: another read on the same question,
  # in units anyone can interpret. Well below 0.90 means too narrow, well above
  # means too wide.
  lo <- apply(yrep, 1, stats::quantile, 0.05)
  hi <- apply(yrep, 1, stats::quantile, 0.95)
  list(ks = unname(kt$statistic), p = kt$p.value, u = u,
       cover90 = mean(y_i >= lo & y_i <= hi))
}

msg("computing randomised PIT by posterior predictive simulation")
rp <- lapply(names(fits), function(k)
  randomised_pit(fits[[k]], SPECS[[k]]$family, grepl("nugget", k)))
names(rp) <- names(fits)
saveRDS(rp, file.path(DIR$interim, "05_randomised_pit.rds"))

cmp <- do.call(rbind, lapply(names(fits), function(k) {
  f <- fits[[k]]
  d <- dispersion_of(f, SPECS[[k]]$family)
  data.frame(
    model = k,
    family = SPECS[[k]]$family,
    waic = f$waic$waic,
    dic = f$dic$dic,
    log_cpo = sum(log(f$cpo$cpo[idx_est]), na.rm = TRUE),
    failed_cpo = sum(f$cpo$failure[idx_est] > 0, na.rm = TRUE),
    pearson_dispersion = d$dispersion,
    rho = d$rho,
    pit_ks_stat = rp[[k]]$ks,
    pit_ks_p = rp[[k]]$p,
    cover90 = rp[[k]]$cover90,
    row.names = NULL)
}))
cmp$n_hyper <- vapply(names(fits),
                      function(k) nrow(fits[[k]]$summary.hyperpar), integer(1))
cmp <- cmp[order(cmp$waic), ]
utils::write.csv(cmp, file.path(DIR$tables, "05_likelihood_comparison.csv"),
                 row.names = FALSE)
msg("likelihood comparison:")
print(cmp, row.names = FALSE, digits = 4)

# ---- Selection rule --------------------------------------------------------
# NOT lowest WAIC. The criteria disagree, and the disagreement is substantive
# enough to state rather than paper over:
#
#   WAIC and DIC prefer binomial + nugget, by about 20 units.
#   Leave-one-out log CPO prefers the beta-binomial, marginally.
#   Randomised PIT decisively prefers the beta-binomial: uniformity cannot be
#     rejected (p ~ 0.33) whereas binomial + nugget is rejected (p ~ 0.05) and
#     plain binomial firmly so (p ~ 0.005).
#
# What this model is FOR settles it. The deliverable is commune prevalence with
# credible intervals at locations the survey never visited, so the quantity that
# matters is whether the predictive distribution for a new cluster is honest -
# which is what the randomised PIT measures and what WAIC, an in-sample
# predictive density, does not. Calibration is therefore the binding criterion,
# with parsimony breaking ties: the log CPO gap between the two beta-binomial
# variants is under 0.15 in total log density across 647 clusters, i.e. nothing,
# and betabinomial_nugget's nugget precision is estimated at ~2000 with a
# credible interval spanning three orders of magnitude - it is unidentified,
# because the nugget and the overdispersion parameter model the same thing.
calibrated <- cmp[!is.na(cmp$pit_ks_p) & cmp$pit_ks_p > 0.05, ]
if (nrow(calibrated) > 0) {
  calibrated <- calibrated[order(calibrated$n_hyper, -calibrated$log_cpo), ]
  best <- calibrated$model[1]
  msg("selected likelihood: ", best,
      " (calibrated: PIT p = ", signif(calibrated$pit_ks_p[1], 3),
      "; most parsimonious among calibrated specs)")
  msg("  note: WAIC would have chosen ", cmp$model[1],
      " (WAIC ", round(cmp$waic[1], 1), " vs ",
      round(calibrated$waic[1], 1), ") - see the comment above for why it does not.")
} else {
  best <- cmp$model[1]
  msg("selected likelihood: ", best,
      " (no specification passed the calibration check; fell back to WAIC)")
}
fit <- fits[[best]]
BEST_FAMILY <- SPECS[[best]]$family
HAS_NUGGET <- grepl("nugget", best)

# Diagnostic figure: the PIT histogram is the picture behind the KS numbers.
# A calibrated model gives a flat histogram; a hump means the predictive
# distribution is too wide, a U-shape means too narrow.
pit_df <- do.call(rbind, lapply(names(rp), function(k) {
  if (is.null(rp[[k]]$u)) return(NULL)
  data.frame(model = k, u = rp[[k]]$u)
}))
if (!is.null(pit_df)) {
  p_pit <- ggplot2::ggplot(pit_df, ggplot2::aes(u)) +
    ggplot2::geom_histogram(bins = 20, fill = "steelblue", colour = "white") +
    ggplot2::geom_hline(yintercept = nrow(dat) / 20, linetype = "dashed") +
    ggplot2::facet_wrap(~model) +
    ggplot2::labs(x = "randomised PIT", y = "clusters",
                  title = "Calibration by likelihood",
                  subtitle = paste("Dashed line = uniform. Predictive distribution",
                                   "for a NEW cluster, not the fitted one.")) +
    ggplot2::theme_minimal(base_size = 9)
  save_fig(p_pit, "05_pit_by_likelihood.png", width = 8, height = 6)
}

saveRDS(fits, file.path(DIR$interim, "05_spde_fits_all.rds"))
saveRDS(fit, file.path(DIR$interim, "05_spde_fit.rds"))

# ===========================================================================
# 6. Results
# ===========================================================================
utils::write.csv(fit$summary.fixed, file.path(DIR$tables, "05_spde_fixed_effects.csv"))
msg("WAIC ", round(fit$waic$waic, 1), " | DIC ", round(fit$dic$dic, 1),
    " | failed CPO ", sum(fit$cpo$failure > 0, na.rm = TRUE))

# Matern range and sd on the km scale of the mesh.
sp <- INLA::inla.spde2.result(fit, "s", spde, do.transform = TRUE)
rng <- INLA::inla.emarginal(function(x) x, sp$marginals.range.nominal[[1]])
sig <- INLA::inla.emarginal(function(x) x, sp$marginals.variance.nominal[[1]])
msg("Matern practical range: ", round(rng, 1), " km | field sd: ",
    round(sqrt(sig), 3), " (logit scale)")
utils::write.csv(data.frame(family = BEST_FAMILY, has_nugget = HAS_NUGGET,
                            range_km = rng, field_sd = sqrt(sig),
                            waic = fit$waic$waic, dic = fit$dic$dic),
                 file.path(DIR$tables, "05_spde_spatial_summary.csv"),
                 row.names = FALSE)

# ---- Grid predictions, by posterior sampling ------------------------------
# The grid is NOT part of the fitted stack, so predictions are formed by drawing
# from the joint posterior and evaluating the linear predictor at the grid:
#
#     eta_g = x_g' beta + A_g s
#
# The cluster-level nugget is deliberately EXCLUDED. It represents idiosyncratic
# variation of a surveyed cluster around the underlying surface, not a property
# of the location, so including it would inflate every commune interval with
# noise that does not belong to the place.
#
# Working draw by draw also preserves the posterior correlation between
# neighbouring cells, which is what makes the aggregated commune intervals in 08
# honest rather than falsely narrow.
N_DRAWS <- 500
msg("drawing ", N_DRAWS, " posterior samples and projecting onto ",
    nrow(grid_use), " grid cells")
samp <- INLA::inla.posterior.sample(N_DRAWS, fit, seed = 20210)

rn <- rownames(samp[[1]]$latent)
s_rows <- grep("^s:", rn)
fx_rows <- vapply(c("intercept", have),
                  function(v) which(rn == paste0(v, ":1")), integer(1))
stopifnot(length(s_rows) == spde$n.spde, !anyNA(fx_rows))

Xg <- as.matrix(cbind(intercept = 1, grid_use[, have, drop = FALSE]))

# A_pred is a sparse Matrix, so its product is an S4 Matrix object; plogis()
# needs a plain numeric vector, hence the coercion INSIDE the call.
draws <- vapply(samp, function(z) {
  lat <- z$latent[, 1]
  eta <- as.numeric(Xg %*% lat[fx_rows]) + as.numeric(A_pred %*% lat[s_rows])
  stats::plogis(eta)
}, numeric(nrow(grid_use)))

grid_out <- data.frame(
  lon = grid_use$lon, lat = grid_use$lat, ADM3_PCODE = grid_use$ADM3_PCODE,
  p_mean = rowMeans(draws),
  p_sd = apply(draws, 1, stats::sd),
  p_lower = apply(draws, 1, stats::quantile, 0.025),
  p_upper = apply(draws, 1, stats::quantile, 0.975))
saveRDS(grid_out, file.path(DIR$processed, "spde_grid_predictions.rds"))
msg("wrote spde_grid_predictions.rds (", nrow(grid_out), " cells) | national mean ",
    round(100 * mean(grid_out$p_mean), 1), "%")

saveRDS(list(coords = grid_out[, c("lon", "lat", "ADM3_PCODE")], draws = draws),
        file.path(DIR$processed, "spde_posterior_draws.rds"))

# ===========================================================================
# 7. Surface map
# ===========================================================================
p <- ggplot2::ggplot(grid_out, ggplot2::aes(lon, lat, fill = p_mean)) +
  ggplot2::geom_raster() +
  ggplot2::coord_sf(crs = CRS_GEO, expand = FALSE) +
  ggplot2::scale_fill_viridis_c(option = "rocket", direction = -1,
                                labels = scales::percent, name = "stunting") +
  ggplot2::labs(title = "Modelled stunting surface, Madagascar 2021",
                subtitle = "Binomial SPDE geostatistical model, 1 km resolution") +
  ggplot2::theme_void(base_size = 9)
save_fig(p, "05_spde_surface.png")

msg("05_model_spde.R complete")
