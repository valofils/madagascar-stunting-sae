# ---------------------------------------------------------------------------
# 02b_covariate_diagnostics.R
#
# Cheap, standalone checks on the commune covariate stack. Reads only the CSV
# written by 02, so it runs in seconds and needs no DHS data and no rasters.
#
# Purpose: confirm the covariates behave before they are trusted in a model,
# and characterise the highland zone the paradox is about. Two things worth
# knowing in advance of the DHS:
#   - whether the hypothesis blocks are collinear enough to be inseparable
#     (if elevation and cold and cropland move together perfectly, the
#     decomposition in 07 cannot attribute the paradox between them), and
#   - what the highland zone actually looks like relative to the lowlands.
#
# Outputs: outputs/tables/02b_*.csv, outputs/figures/02b_*.png
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("dplyr", "tidyr", "ggplot2", "sf", "scales", "patchwork")

stopifnot(file.exists(OUT$cov_commune))
cc <- utils::read.csv(OUT$cov_commune)
msg("commune covariate stack: ", nrow(cc), " x ", ncol(cc))

HIGHLAND_M <- 800

# ===========================================================================
# 1. Integrity checks
# ===========================================================================
fr <- grep("^frac_", names(cc), value = TRUE)
tot <- rowSums(cc[fr], na.rm = TRUE)
chk <- function(cond, what) {
  if (!isTRUE(cond)) warning("CHECK FAILED: ", what, call. = FALSE) else msg("ok: ", what)
}
chk(all(abs(tot - 1) < 0.01), "land-cover fractions sum to 1 in every commune")
chk(all(cc$elevation >= 0 & cc$elevation < 3000), "elevation is in a plausible range")
chk(all(cc$pop_u5 > 0), "every commune has a positive under-5 population")
chk(!anyNA(cc$elevation), "no missing elevation")

# Coverage: which covariates actually made it into the stack.
cover <- data.frame(variable = names(cc),
                    pct_missing = round(100 * colMeans(is.na(cc)), 2))
missing_planned <- setdiff(
  c("cattle_density", "travel_time", "nightlights"), names(cc))
if (length(missing_planned) > 0)
  msg("NOT YET AVAILABLE (manual download, see 02_covariates_helpers.R): ",
      paste(missing_planned, collapse = ", "))
utils::write.csv(cover, file.path(DIR$tables, "02b_covariate_coverage.csv"),
                 row.names = FALSE)

# ===========================================================================
# 2. Highland vs lowland profile
# ===========================================================================
cc$zone <- ifelse(cc$elevation > HIGHLAND_M, "Highland (>800 m)", "Lowland")

VARS <- intersect(c("elevation", "ruggedness", "temp_min_cold", "temp_mean",
                    "temp_seasonality", "precip_annual", "precip_seasonality",
                    "frac_crop", "frac_tree", "frac_grass", "frac_built",
                    "frac_water_perm", "pop_dens", "cattle_density",
                    "travel_time", "travel_time_healthcare", "nightlights"),
                  names(cc))

# Population-weighted, because the question is what the average CHILD is
# exposed to, not what the average commune polygon looks like.
# The weight is carried in its own column (wt_u5): summarise() evaluates its
# arguments in order, so reusing the name pop_u5 for the group total would
# replace the weight vector with a scalar before across() ever sees it.
profile <- cc |>
  dplyr::mutate(wt_u5 = pop_u5) |>
  dplyr::group_by(zone) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(VARS),
                                 ~ stats::weighted.mean(.x, wt_u5, na.rm = TRUE)),
                   communes = dplyr::n(),
                   pop_u5 = sum(wt_u5),
                   .groups = "drop") |>
  dplyr::relocate(zone, communes, pop_u5)
utils::write.csv(profile, file.path(DIR$tables, "02b_highland_profile.csv"),
                 row.names = FALSE)
msg("population-weighted highland / lowland profile:")
print(as.data.frame(profile), row.names = FALSE, digits = 4)

msg("share of Madagascar's under-5 population living above ", HIGHLAND_M, " m: ",
    round(100 * sum(cc$pop_u5[cc$zone != "Lowland"]) / sum(cc$pop_u5), 1), "%")

# ===========================================================================
# 3. Collinearity between hypothesis blocks
# ===========================================================================
# If the H1 and H2 covariates are near-perfectly correlated, no model can
# apportion the paradox between cold stress and diet - the decomposition in 07
# would be reporting the order the blocks were entered, not evidence. This
# table is the honest precondition for that analysis.
cm <- stats::cor(cc[VARS], use = "pairwise.complete.obs")
utils::write.csv(round(cm, 3), file.path(DIR$tables, "02b_covariate_correlations.csv"))

pairs_hi <- which(abs(cm) > 0.8 & upper.tri(cm), arr.ind = TRUE)
if (nrow(pairs_hi) > 0) {
  hi <- data.frame(a = rownames(cm)[pairs_hi[, 1]],
                   b = colnames(cm)[pairs_hi[, 2]],
                   r = round(cm[pairs_hi], 3))
  msg("covariate pairs with |r| > 0.8 (watch these in the decomposition):")
  print(hi, row.names = FALSE)
  utils::write.csv(hi, file.path(DIR$tables, "02b_high_correlations.csv"),
                   row.names = FALSE)
} else {
  msg("no covariate pair exceeds |r| = 0.8 - blocks are separable")
}

# Correlation of each covariate with elevation: how much of the covariate set
# is really just altitude in disguise.
with_elev <- data.frame(variable = VARS,
                        r_with_elevation = round(cm["elevation", VARS], 3))
with_elev <- with_elev[order(-abs(with_elev$r_with_elevation)), ]
utils::write.csv(with_elev, file.path(DIR$tables, "02b_correlation_with_elevation.csv"),
                 row.names = FALSE)
print(with_elev, row.names = FALSE)

# ===========================================================================
# 4. Figures
# ===========================================================================
adm3 <- sf::st_read(OUT$adm3, quiet = TRUE)
adm1 <- sf::st_read(OUT$adm1, quiet = TRUE)
# adm3 already carries pop_dens, pop_u5 and friends from 01, so join onto the
# geometry and the key alone - otherwise dplyr silently suffixes the duplicated
# columns to .x/.y and the plotting aesthetics stop resolving.
g <- dplyr::left_join(adm3[, "ADM3_PCODE"], cc[, c("ADM3_PCODE", VARS, "zone")],
                      by = "ADM3_PCODE")

theme_map <- ggplot2::theme_void(base_size = 8) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                 legend.key.width = ggplot2::unit(0.3, "cm"))

one <- function(var, ttl, opt = "viridis", trans = "identity") {
  ggplot2::ggplot(g) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[var]]), colour = NA) +
    ggplot2::geom_sf(data = adm1, fill = NA, colour = "white", linewidth = 0.15) +
    ggplot2::scale_fill_viridis_c(option = opt, name = NULL, trans = trans,
                                  labels = scales::label_number()) +
    ggplot2::labs(title = ttl) + theme_map
}

# Population density spans four orders of magnitude (Antananarivo at ~20,000
# per km2 against a national median near 54), so on a linear scale every
# commune but the capital renders as the same dark colour. Log10 is the only
# way this panel carries information.
panel <- one("elevation", "Elevation (m)", "cividis") +
  one("temp_min_cold", "Min temp, coldest month (C)", "mako") +
  one("frac_crop", "Cropland fraction", "viridis") +
  one("pop_dens", "Population density (per km2, log scale)", "magma",
      trans = "log10") +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "The highland signature: cold, cropped and crowded",
    subtitle = paste("The central plateau is at once the highest, the coldest",
                     "and by far the most cultivated part of the country. It is",
                     "also densely settled, though the eastern littoral is too."),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13)))
save_fig(panel, "02b_highland_signature.png", width = 16, height = 7)

# Population density is extremely skewed (Antananarivo), so plot it on a log
# scale or the whole distribution collapses into one bar.
long <- cc |>
  dplyr::select(zone, pop_u5, dplyr::all_of(
    intersect(c("temp_min_cold", "frac_crop", "precip_seasonality", "ruggedness"),
              VARS))) |>
  tidyr::pivot_longer(-c(zone, pop_u5))

p_dist <- ggplot2::ggplot(long, ggplot2::aes(value, fill = zone, weight = pop_u5)) +
  ggplot2::geom_density(alpha = 0.5, colour = NA) +
  ggplot2::facet_wrap(~name, scales = "free", ncol = 2) +
  ggplot2::scale_fill_manual(values = c("Highland (>800 m)" = "#B2182B",
                                        "Lowland" = "#2166AC"), name = NULL) +
  ggplot2::labs(x = NULL, y = "density (under-5 population weighted)",
                title = "How highland and lowland exposures differ",
                subtitle = "Weighted by under-5 population, not by area") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(legend.position = "bottom")
save_fig(p_dist, "02b_highland_lowland_distributions.png", width = 8, height = 6)

msg("02b_covariate_diagnostics.R complete")
