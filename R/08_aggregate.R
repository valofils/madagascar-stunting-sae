# ---------------------------------------------------------------------------
# 08_aggregate.R
#
# Turn the 1 km stunting surface from 05 into commune-level estimates.
#
# The aggregation is population-weighted, over posterior DRAWS rather than over
# the posterior mean. This matters. A commune's prevalence is
#
#     P_c = sum_g w_g p_g   with   w_g = under-5 population in grid cell g,
#                                        normalised within the commune
#
# Averaging the posterior mean surface gives the right point estimate but the
# wrong uncertainty: neighbouring cells share the same spatial field, so their
# errors are strongly positively correlated and do NOT average away. Working
# draw by draw preserves that correlation, so the credible interval is honest.
#
# The same routine then aggregates to district and region, where the results
# are benchmarked against the design-based direct estimates in 09.
#
# Inputs : data/processed/spde_posterior_draws.rds  (from 05)
#          WorldPop under-5 raster (cached by 01)
#          data/processed/adm3_communes.gpkg
# Outputs: data/processed/commune_stunting.csv / .gpkg
#          data/processed/aggregated_adm{1,2}.csv
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("dplyr", "sf", "terra", "ggplot2")

f_draws <- file.path(DIR$processed, "spde_posterior_draws.rds")
stopifnot(file.exists(f_draws))
post <- readRDS(f_draws)
draws <- post$draws              # cells x draws
coords <- post$coords            # lon, lat, ADM3_PCODE
msg("posterior draws: ", ncol(draws), " over ", nrow(draws), " grid cells")

adm3 <- sf::st_read(OUT$adm3, quiet = TRUE)

# ===========================================================================
# 1. Population weight for every grid cell
# ===========================================================================
# The surface is a per-child risk; a commune's prevalence is the risk averaged
# over its CHILDREN, not over its area. Weighting by under-5 counts is what
# makes the estimate comparable to a DHS prevalence.
wp <- file.path(DIR$rasters, "worldpop")
u5_files <- file.path(wp, c("mdg_f_0_2020_constrained.tif", "mdg_f_1_2020_constrained.tif",
                            "mdg_m_0_2020_constrained.tif", "mdg_m_1_2020_constrained.tif"))
stopifnot(all(file.exists(u5_files)))
u5 <- sum(terra::rast(u5_files), na.rm = TRUE)

cellpts <- terra::vect(as.matrix(coords[, c("lon", "lat")]), type = "points",
                       crs = "EPSG:4326")
# The grid is ~1 km and WorldPop is 100 m, so a point lookup would sample one
# hundredth of the cell. Sum the population within each grid cell instead.
res_deg <- stats::median(diff(sort(unique(coords$lon))))
sq <- terra::rast(terra::ext(min(coords$lon) - res_deg, max(coords$lon) + res_deg,
                             min(coords$lat) - res_deg, max(coords$lat) + res_deg),
                  resolution = res_deg, crs = "EPSG:4326")
u5_grid <- terra::resample(u5, sq, method = "sum")
w <- terra::extract(u5_grid, cellpts, ID = FALSE)[, 1]
w[is.na(w) | w < 0] <- 0
msg("grid cells with zero under-5 population: ", sum(w == 0), " of ", length(w))

# A commune whose every cell has zero modelled population would otherwise
# produce NaN. Fall back to equal weights inside such communes.
wdf <- data.frame(ADM3_PCODE = coords$ADM3_PCODE, w = w)
zero_com <- wdf |>
  dplyr::group_by(ADM3_PCODE) |>
  dplyr::summarise(tot = sum(w), .groups = "drop") |>
  dplyr::filter(tot <= 0)
if (nrow(zero_com) > 0) {
  msg("communes with no population weight, using equal weights: ", nrow(zero_com))
  w[coords$ADM3_PCODE %in% zero_com$ADM3_PCODE] <- 1
}

# ===========================================================================
# 2. Weighted aggregation, draw by draw
# ===========================================================================
aggregate_draws <- function(area_id, weights, draws) {
  keep <- !is.na(area_id) & weights > 0
  area_id <- area_id[keep]; weights <- weights[keep]; draws <- draws[keep, , drop = FALSE]

  f <- factor(area_id)
  # rowsum() gives, for each area, the population-weighted sum over its cells,
  # in one pass per draw - fast enough for 1,700 areas x 500 draws.
  num <- rowsum(draws * weights, group = f, reorder = TRUE)
  den <- rowsum(weights, group = f, reorder = TRUE)[, 1]
  post_area <- num / den                      # areas x draws

  data.frame(
    area = rownames(post_area),
    est = apply(post_area, 1, mean),
    median = apply(post_area, 1, stats::median),
    sd = apply(post_area, 1, stats::sd),
    lower = apply(post_area, 1, stats::quantile, 0.025),
    upper = apply(post_area, 1, stats::quantile, 0.975),
    n_cells = as.vector(table(f)),
    row.names = NULL)
}

msg("aggregating to communes")
com <- aggregate_draws(coords$ADM3_PCODE, w, draws)
names(com)[1] <- "ADM3_PCODE"

com <- com |>
  dplyr::left_join(sf::st_drop_geometry(adm3), by = "ADM3_PCODE") |>
  dplyr::mutate(cv = sd / est,
                ci_width = upper - lower,
                n_stunted = est * pop_u5)

msg("communes estimated: ", nrow(com), " of ", nrow(adm3))
msg("prevalence range: ", round(100 * min(com$est), 1), "% to ",
    round(100 * max(com$est), 1), "%")
msg("median CI width: ", round(100 * stats::median(com$ci_width), 1), " points")
msg("communes with CV > 0.30 (conventionally unpublishable): ",
    sum(com$cv > 0.30, na.rm = TRUE))

utils::write.csv(com, file.path(DIR$processed, "commune_stunting.csv"),
                 row.names = FALSE)
sf::st_write(dplyr::left_join(adm3, com[, setdiff(names(com), names(sf::st_drop_geometry(adm3)))],
                              by = "ADM3_PCODE"),
             file.path(DIR$processed, "commune_stunting.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

# ===========================================================================
# 3. Higher levels, for benchmarking
# ===========================================================================
lookup <- sf::st_drop_geometry(adm3)[, c("ADM3_PCODE", "ADM2_PCODE", "ADM1_PCODE")]
map2 <- stats::setNames(lookup$ADM2_PCODE, lookup$ADM3_PCODE)
map1 <- stats::setNames(lookup$ADM1_PCODE, lookup$ADM3_PCODE)

msg("aggregating to districts and regions")
d2 <- aggregate_draws(unname(map2[coords$ADM3_PCODE]), w, draws)
names(d2)[1] <- "ADM2_PCODE"
d1 <- aggregate_draws(unname(map1[coords$ADM3_PCODE]), w, draws)
names(d1)[1] <- "ADM1_PCODE"

utils::write.csv(d2, file.path(DIR$processed, "aggregated_adm2.csv"), row.names = FALSE)
utils::write.csv(d1, file.path(DIR$processed, "aggregated_adm1.csv"), row.names = FALSE)

# National figure, weighting communes by their under-5 population. Compared in
# 09 against the design-based national direct estimate.
nat <- stats::weighted.mean(com$est, com$pop_u5, na.rm = TRUE)
msg("model-implied national stunting: ", round(100 * nat, 1), "%")

# ===========================================================================
# 4. Highland vs lowland, the headline contrast
# ===========================================================================
if (file.exists(OUT$cov_commune)) {
  cc <- utils::read.csv(OUT$cov_commune)
  if ("elevation" %in% names(cc)) {
    com2 <- dplyr::left_join(com, cc[, c("ADM3_PCODE", "elevation")], by = "ADM3_PCODE")
    com2$zone <- ifelse(com2$elevation > 800, "Highland (>800 m)", "Lowland")
    summ <- com2 |>
      dplyr::group_by(zone) |>
      dplyr::summarise(communes = dplyr::n(),
                       pop_u5 = sum(pop_u5, na.rm = TRUE),
                       stunting = stats::weighted.mean(est, pop_u5, na.rm = TRUE),
                       children_stunted = sum(n_stunted, na.rm = TRUE),
                       .groups = "drop")
    utils::write.csv(summ, file.path(DIR$tables, "08_highland_lowland.csv"),
                     row.names = FALSE)
    print(as.data.frame(summ), row.names = FALSE)
  }
}

# ===========================================================================
# 5. Targeting table: where the stunted children actually are
# ===========================================================================
# Prevalence alone is a poor targeting rule: a very high rate in an empty
# commune affects few children. Ranking by absolute burden is what a health
# ministry allocating basic health centres needs.
top <- com |>
  dplyr::arrange(dplyr::desc(n_stunted)) |>
  dplyr::mutate(cum_share = cumsum(n_stunted) / sum(n_stunted, na.rm = TRUE)) |>
  dplyr::select(ADM3_PCODE, ADM3_EN, ADM2_EN, ADM1_EN, est, lower, upper,
                pop_u5, n_stunted, cum_share)
utils::write.csv(top, file.path(DIR$tables, "08_commune_burden_ranking.csv"),
                 row.names = FALSE)
n50 <- which(top$cum_share >= 0.5)[1]
msg(n50, " communes (", round(100 * n50 / nrow(top), 1),
    "% of all communes) contain half of Madagascar's stunted children")

msg("08_aggregate.R complete")
