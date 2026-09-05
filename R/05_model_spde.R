# ---------------------------------------------------------------------------
# 05_model_spde.R
#
# Continuous-surface model: binomial geostatistical regression on DHS cluster
# locations, with a Matern spatial field represented through the SPDE
# approximation of Lindgren, Rue & Lindstrom (2011) and fitted by INLA.
#
# Model:
#   y_c ~ Binomial(n_c, p_c)                     y_c = stunted children in cluster c
#   logit(p_c) = alpha + x_c' beta + S(s_c) + e_c
#   S(.)  ~ GP with Matern covariance (SPDE), nu = 1
#   e_c   ~ N(0, sigma_e^2)                      cluster-level nugget
#
# The nugget matters: without it the spatial field absorbs cluster-level
# idiosyncrasy and the effective range collapses, which in turn overstates the
# precision of the predicted surface.
#
# The fitted field is then projected onto a prediction grid covering
# Madagascar, giving the 1 km surface that 08_aggregate.R averages up to
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
            "frac_water_perm", "cattle_density", "travel_time", "nightlights",
            "urban")
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
# 1 km grid over the land mass; covariates come from the same commune-level
# raster stack, resampled to the grid.
GRID_KM <- 1
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
  saveRDS(grid, grid_cov_file)
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
# 5. Stacks and fit
# ===========================================================================
A_est <- INLA::inla.spde.make.A(mesh, loc = coo)
A_pred <- INLA::inla.spde.make.A(mesh, loc = coo_pred)
idx <- INLA::inla.spde.make.index("s", n.spde = spde$n.spde)

fixed_df <- function(d) {
  x <- data.frame(intercept = 1)
  x <- cbind(x[rep(1, nrow(d)), , drop = FALSE], d[, have, drop = FALSE])
  rownames(x) <- NULL
  x
}

stk_est <- INLA::inla.stack(
  tag = "est", data = list(y = dat$y, n = dat$n),
  A = list(A_est, 1),
  effects = list(idx, cbind(fixed_df(dat), nugget = seq_len(nrow(dat)))))

stk_pred <- INLA::inla.stack(
  tag = "pred", data = list(y = NA, n = NA),
  A = list(A_pred, 1),
  effects = list(idx, cbind(fixed_df(grid_use), nugget = NA)))

stk <- INLA::inla.stack(stk_est, stk_pred)

form <- stats::as.formula(paste(
  "y ~ 0 + intercept +", paste(have, collapse = " + "),
  "+ f(s, model = spde)",
  "+ f(nugget, model = 'iid', hyper = list(prec = list(prior = 'pc.prec',",
  "                                                    param = c(1, 0.01))))"))

msg("fitting INLA SPDE model (this is the slow step)")
t0 <- Sys.time()
fit <- INLA::inla(
  form, family = "binomial", Ntrials = INLA::inla.stack.data(stk)$n,
  data = INLA::inla.stack.data(stk, spde = spde),
  control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE, link = 1),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  control.inla = list(int.strategy = "eb"),
  verbose = FALSE)
msg("INLA finished in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

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
utils::write.csv(data.frame(range_km = rng, field_sd = sqrt(sig),
                            waic = fit$waic$waic, dic = fit$dic$dic),
                 file.path(DIR$tables, "05_spde_spatial_summary.csv"),
                 row.names = FALSE)

# ---- Grid predictions ------------------------------------------------------
ip <- INLA::inla.stack.index(stk, "pred")$data
grid_out <- data.frame(
  lon = grid_use$lon, lat = grid_use$lat, ADM3_PCODE = grid_use$ADM3_PCODE,
  p_mean = fit$summary.fitted.values[ip, "mean"],
  p_sd = fit$summary.fitted.values[ip, "sd"],
  p_lower = fit$summary.fitted.values[ip, "0.025quant"],
  p_upper = fit$summary.fitted.values[ip, "0.975quant"])
saveRDS(grid_out, file.path(DIR$processed, "spde_grid_predictions.rds"))
msg("wrote spde_grid_predictions.rds (", nrow(grid_out), " cells)")

# ---- Posterior draws, for correct uncertainty when aggregating ------------
# Aggregating the posterior MEAN to communes would give the right point
# estimate but the wrong interval, because it ignores the spatial correlation
# between neighbouring cells. 08_aggregate.R therefore needs draws, not means.
msg("drawing 500 posterior samples for aggregation")
samp <- INLA::inla.posterior.sample(500, fit, seed = 20210)
pred_rows <- ip
draws <- vapply(samp, function(s) stats::plogis(s$latent[pred_rows, 1]),
                numeric(length(pred_rows)))
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
