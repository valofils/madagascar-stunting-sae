# ---------------------------------------------------------------------------
# 07_explain_paradox.R
#
# The explanatory core of the project: why is stunting highest in the
# agro-ecologically richest part of the country?
#
# Strategy. Fit a sequence of nested models for continuous HAZ at cluster
# level, each adding one hypothesis block, and watch what happens to a single
# quantity: the raw highland penalty. Define
#
#     HIGHLAND = cluster elevation > 800 m (the Hauts Plateaux threshold)
#
# and let delta_k be the coefficient on HIGHLAND in model k. Model 0 contains
# HIGHLAND alone, so delta_0 is the paradox as observed. Each subsequent model
# adds a hypothesis block, and the drop in delta is the share of the paradox
# that block accounts for:
#
#     explained_k = (delta_0 - delta_k) / delta_0
#
# This is a mediation-style decomposition, not a causal identification. It
# answers "what covariates make the highland penalty disappear", which is
# exactly the question the World Bank paper left open, but a covariate that
# absorbs the penalty is a candidate mechanism, not a proven cause. Blocks are
# also entered one at a time (not just cumulatively) because they are
# correlated with each other, and the order of entry would otherwise decide
# the answer.
#
# Every model carries the SPDE spatial field, so delta measures the highland
# effect over and above smooth spatial structure - a highland penalty that
# survives the spatial field is a genuine altitude signal rather than a
# rebadged "central Madagascar" indicator.
#
# Inputs : data/interim/dhs_child_haz.rds, data/processed/cluster_covariates.csv
# Outputs: outputs/tables/07_paradox_decomposition.csv
#          outputs/tables/07_covariate_effects.csv
#          outputs/figures/07_*.png
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("dplyr", "sf", "ggplot2", "INLA")

set.seed(20210)
HIGHLAND_M <- 800   # conventional lower bound of the Hauts Plateaux

# ===========================================================================
# 1. Cluster-level HAZ
# ===========================================================================
stopifnot(file.exists(OUT$dhs_child), file.exists(OUT$cov_cluster))
child <- readRDS(OUT$dhs_child)
clu <- utils::read.csv(OUT$cov_cluster)

# Mean HAZ per cluster, with its sampling variance. Modelling the cluster mean
# rather than the child keeps the spatial model tractable; child-level
# compositional variables enter as cluster means too (see block D).
dat <- child |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_haz = mean(haz, na.rm = TRUE),
    sd_haz = stats::sd(haz, na.rm = TRUE),
    stunting = mean(stunted, na.rm = TRUE),
    # Block D, care and maternal condition, averaged within cluster
    mother_edu = mean(mother_edu, na.rm = TRUE),
    wealth_q = mean(wealth_q, na.rm = TRUE),
    mother_bmi = mean(mother_bmi, na.rm = TRUE),
    mother_age1b = mean(mother_age1b, na.rm = TRUE),
    birth_order = mean(birth_order, na.rm = TRUE),
    # H2, from what children actually ate (see 03, section 2b)
    diet_diversity = mean(diet_diversity, na.rm = TRUE),
    mdd = mean(mdd, na.rm = TRUE),
    # H3, direct infection signal rather than settlement density alone
    diarrhea_2w = mean(diarrhea_2w == 1, na.rm = TRUE),
    fever_2w = mean(fever_2w == 1, na.rm = TRUE),
    improved_water = mean(improved_water, na.rm = TRUE),
    basic_sanitation = mean(basic_sanitation, na.rm = TRUE),
    open_defecation = mean(open_defecation, na.rm = TRUE),
    pct_male = mean(sex == 1, na.rm = TRUE),
    mean_age_m = mean(age_month, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::inner_join(clu, by = c("cluster" = "DHSCLUST")) |>
  dplyr::filter(n >= 3, !is.na(mean_haz), !is.na(elevation))

msg("clusters in the explanatory model: ", nrow(dat))
msg("highland clusters (>", HIGHLAND_M, " m): ", sum(dat$elevation > HIGHLAND_M),
    " | lowland: ", sum(dat$elevation <= HIGHLAND_M))

dat$highland <- as.integer(dat$elevation > HIGHLAND_M)

# The raw contrast, before any adjustment - the paradox in one number.
raw_gap <- mean(dat$mean_haz[dat$highland == 1], na.rm = TRUE) -
  mean(dat$mean_haz[dat$highland == 0], na.rm = TRUE)
msg("unadjusted highland-lowland HAZ gap: ", round(raw_gap, 3), " SD")

# ===========================================================================
# 2. Hypothesis blocks
# ===========================================================================
# Each block names the covariates that operationalise one hypothesis. Only
# those actually present in the covariate file are used, so the script still
# runs when a raster could not be obtained.
# Note on what each block can and cannot carry, given 02b:
#   A  altitude and cold are collinear at r = -0.90; this block is a single
#      altitude-temperature construct and must be reported as such.
#   B  cattle_density behaves as a WEALTH proxy here (highlands have MORE
#      cattle), so diet_diversity / mdd from the KR recode is what actually
#      tests H2. Cattle is retained but must not be read as food access.
#   C  frac_built is kept and nightlights deliberately excluded: they correlate
#      0.95, and nightlights sits conceptually in H4, so including both would
#      inflate whichever block enters first.
#   D  the highlands are BETTER served than the lowlands, so this block is
#      expected to WIDEN the highland penalty, not absorb it. A negative
#      "explained" share here is the correct result, not a bug.
BLOCKS <- list(
  A_cold_stress = c("temp_min_cold", "temp_seasonality", "ruggedness"),
  B_diet        = c("diet_diversity", "mdd", "frac_crop", "cattle_density",
                    "precip_annual", "precip_seasonality"),
  C_infection   = c("pop_count", "frac_built", "frac_water_perm",
                    "diarrhea_2w", "fever_2w", "improved_water",
                    "basic_sanitation", "open_defecation"),
  D_care        = c("mother_edu", "wealth_q", "mother_bmi", "mother_age1b",
                    "birth_order", "travel_time", "travel_time_healthcare")
)
BLOCK_LABEL <- c(
  A_cold_stress = "H1 Altitude / cold stress",
  B_diet        = "H2 Rice monoculture / diet quality",
  C_infection   = "H3 Infection load / enteropathy",
  D_care        = "H4 Care practices, maternal condition, access")

BLOCKS <- lapply(BLOCKS, function(v) {
  v <- intersect(v, names(dat))
  v[vapply(v, function(x) stats::sd(dat[[x]], na.rm = TRUE) > 0 &&
             mean(is.na(dat[[x]])) < 0.1, logical(1))]
})
for (b in names(BLOCKS))
  msg(BLOCK_LABEL[[b]], ": ", paste(BLOCKS[[b]], collapse = ", "))

# Always adjusted for: child age and sex composition drive HAZ mechanically
# and have nothing to do with the paradox.
CONTROLS <- intersect(c("mean_age_m", "pct_male", "urban"), names(dat))

all_vars <- unique(c(unlist(BLOCKS), CONTROLS))
ctr <- vapply(dat[all_vars], mean, numeric(1), na.rm = TRUE)
scl <- vapply(dat[all_vars], stats::sd, numeric(1), na.rm = TRUE)
for (v in all_vars) {
  dat[[v]] <- (dat[[v]] - ctr[[v]]) / scl[[v]]
  dat[[v]][is.na(dat[[v]])] <- 0    # standardised: 0 is the mean
}

# ===========================================================================
# 3. Shared spatial structure
# ===========================================================================
pts <- sf::st_as_sf(dat, coords = c("lon", "lat"), crs = CRS_GEO) |>
  sf::st_transform(CRS_EQA)
coo <- sf::st_coordinates(pts) / 1000

adm0 <- sf::st_read(OUT$adm1, quiet = TRUE) |>
  sf::st_transform(CRS_EQA) |> sf::st_union() |> sf::st_make_valid()
bnd <- sf::st_coordinates(sf::st_simplify(adm0, dTolerance = 5000))[, 1:2] / 1000

mesh <- INLA::inla.mesh.2d(loc = coo,
                           boundary = INLA::inla.nonconvex.hull(bnd, convex = -0.03),
                           max.edge = c(30, 120), cutoff = 10, offset = c(20, 120))
spde <- INLA::inla.spde2.pcmatern(mesh, prior.range = c(50, 0.05),
                                  prior.sigma = c(1, 0.05))
A <- INLA::inla.spde.make.A(mesh, loc = coo)
idx <- INLA::inla.spde.make.index("s", n.spde = spde$n.spde)
msg("mesh nodes: ", mesh$n)

# The cluster mean HAZ has sampling variance sigma^2 / n, so clusters with more
# children are more informative. Passing that as the Gaussian scale gives each
# cluster its correct weight instead of treating all clusters as equal.
fit_model <- function(vars, label) {
  vars <- unique(c("highland", CONTROLS, vars))
  X <- data.frame(intercept = 1, dat[, vars, drop = FALSE])
  stk <- INLA::inla.stack(
    tag = "est", data = list(y = dat$mean_haz),
    A = list(A, 1), effects = list(idx, X))
  form <- stats::as.formula(paste("y ~ 0 + intercept +",
                                  paste(vars, collapse = " + "),
                                  "+ f(s, model = spde)"))
  f <- INLA::inla(form, family = "gaussian",
                  data = INLA::inla.stack.data(stk, spde = spde),
                  scale = dat$n,
                  control.predictor = list(A = INLA::inla.stack.A(stk), compute = FALSE),
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                  control.inla = list(int.strategy = "eb"))
  msg("  fitted [", label, "]  WAIC ", round(f$waic$waic, 1),
      "  highland ", round(f$summary.fixed["highland", "mean"], 3))
  f
}

# ===========================================================================
# 4. The decomposition
# ===========================================================================
msg("=== fitting the nested sequence ===")
fits <- list()
fits[["M0_paradox"]] <- fit_model(character(0), "M0: highland only")

# (a) one block at a time: how much of the gap can each hypothesis carry alone?
for (b in names(BLOCKS))
  if (length(BLOCKS[[b]]) > 0)
    fits[[paste0("M1_", b)]] <- fit_model(BLOCKS[[b]], paste("M1:", b))

# (b) cumulative, in hypothesis order
cum <- character(0)
for (b in names(BLOCKS)) {
  if (length(BLOCKS[[b]]) == 0) next
  cum <- c(cum, BLOCKS[[b]])
  fits[[paste0("M2_cum_", b)]] <- fit_model(cum, paste("M2: cumulative through", b))
}

delta0 <- fits[["M0_paradox"]]$summary.fixed["highland", "mean"]

decomp <- do.call(rbind, lapply(names(fits), function(k) {
  s <- fits[[k]]$summary.fixed["highland", ]
  data.frame(
    model = k,
    highland_coef = s[["mean"]],
    highland_sd = s[["sd"]],
    ci_low = s[["0.025quant"]],
    ci_high = s[["0.975quant"]],
    pct_of_gap_explained = 100 * (delta0 - s[["mean"]]) / delta0,
    waic = fits[[k]]$waic$waic,
    dic = fits[[k]]$dic$dic,
    row.names = NULL)
}))
decomp$block_label <- BLOCK_LABEL[sub("^M[12]_(cum_)?", "", decomp$model)]

utils::write.csv(decomp, file.path(DIR$tables, "07_paradox_decomposition.csv"),
                 row.names = FALSE)
msg("decomposition written; unadjusted highland coefficient = ", round(delta0, 3))
print(decomp[, c("model", "highland_coef", "pct_of_gap_explained", "waic")],
      row.names = FALSE)

# ===========================================================================
# 5. Full covariate effect table (the most complete model)
# ===========================================================================
full_key <- utils::tail(grep("^M2_cum_", names(fits), value = TRUE), 1)
full <- fits[[full_key]]

eff <- full$summary.fixed
eff$variable <- rownames(eff)
eff$block <- vapply(eff$variable, function(v) {
  hit <- names(BLOCKS)[vapply(BLOCKS, function(b) v %in% b, logical(1))]
  if (length(hit) == 0) {
    if (v == "highland") "paradox term" else if (v %in% CONTROLS) "control" else "intercept"
  } else BLOCK_LABEL[[hit[1]]]
}, character(1))
# Coefficients are per standard deviation of the covariate, in HAZ SD units.
eff$significant <- sign(eff[["0.025quant"]]) == sign(eff[["0.975quant"]])
eff <- eff[order(-abs(eff$mean)), ]
utils::write.csv(eff, file.path(DIR$tables, "07_covariate_effects.csv"),
                 row.names = FALSE)

sp <- INLA::inla.spde2.result(full, "s", spde, do.transform = TRUE)
msg("residual spatial range in the full model: ",
    round(INLA::inla.emarginal(function(x) x, sp$marginals.range.nominal[[1]]), 1),
    " km")

saveRDS(fits, file.path(DIR$interim, "07_paradox_fits.rds"))

# ===========================================================================
# 6. Figures
# ===========================================================================
## (a) how the highland penalty shrinks as blocks are added
dd <- decomp[grepl("^M0|^M2_cum_", decomp$model), ]
dd$step <- factor(dd$model, levels = dd$model,
                  labels = c("Unadjusted", sub("^M2_cum_", "+ ", dd$model[-1])))
p1 <- ggplot2::ggplot(dd, ggplot2::aes(step, highland_coef)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  ggplot2::geom_point(size = 3, colour = "firebrick") +
  ggplot2::labs(x = NULL, y = "Highland effect on HAZ (SD)",
                title = "How much of the highland penalty each hypothesis absorbs",
                subtitle = paste0("Cumulative adjustment; unadjusted gap = ",
                                  round(delta0, 2), " HAZ SD")) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
save_fig(p1, "07_paradox_decomposition.png", width = 7.5, height = 5)

## (b) covariate effects in the full model
ef <- eff[eff$variable != "intercept", ]
p2 <- ggplot2::ggplot(ef, ggplot2::aes(stats::reorder(variable, mean), mean,
                                       colour = block)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = `0.025quant`, ymax = `0.975quant`),
                         width = 0) +
  ggplot2::geom_point(size = 2) +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Effect on HAZ (SD per SD of covariate)",
                colour = NULL,
                title = "Covariate effects, fully adjusted spatial model") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(legend.position = "bottom")
save_fig(p2, "07_covariate_effects.png", width = 8, height = 7)

## (c) the paradox itself: HAZ against elevation, and against a wealth proxy
p3 <- ggplot2::ggplot(dat, ggplot2::aes(elevation * scl[["elevation"]] + ctr[["elevation"]],
                                        mean_haz)) +
  ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.3) +
  ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                       colour = "firebrick") +
  ggplot2::geom_vline(xintercept = HIGHLAND_M, linetype = "dashed") +
  ggplot2::labs(x = "Elevation (m)", y = "Cluster mean HAZ", size = "children",
                title = "The fertile-highland paradox",
                subtitle = "Dashed line: 800 m, the conventional edge of the Hauts Plateaux") +
  ggplot2::theme_minimal(base_size = 10)
save_fig(p3, "07_haz_vs_elevation.png", width = 7, height = 5)

msg("07_explain_paradox.R complete")
