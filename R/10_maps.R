# ---------------------------------------------------------------------------
# 10_maps.R
#
# Final figures: the commune stunting map with its uncertainty, the comparison
# against the two precedents, and the paradox panel.
#
# Design choices worth stating, because they change what a reader concludes:
#   - Prevalence uses a sequential palette on a fixed scale across every map,
#     so panels are visually comparable.
#   - Uncertainty is shown as its own map, not hidden. A commune map with 1,700
#     units and no uncertainty layer invites over-reading of small differences.
#   - The burden map (count of stunted children) is shown next to prevalence,
#     because they imply different targeting and are routinely confused.
#
# Inputs : data/processed/commune_stunting.gpkg and friends
# Outputs: outputs/figures/10_*.png
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("dplyr", "sf", "ggplot2", "scales", "patchwork")

f_com <- file.path(DIR$processed, "commune_stunting.gpkg")
stopifnot(file.exists(f_com))
com <- sf::st_read(f_com, quiet = TRUE)
adm1 <- sf::st_read(OUT$adm1, quiet = TRUE)

theme_map <- ggplot2::theme_void(base_size = 9) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey30"),
    legend.key.width = ggplot2::unit(0.35, "cm"),
    legend.key.height = ggplot2::unit(1.1, "cm"))

region_lines <- ggplot2::geom_sf(data = adm1, fill = NA, colour = "white",
                                 linewidth = 0.25)

# A single fixed prevalence scale for every prevalence map in the output.
PREV_LIMITS <- range(com$est, na.rm = TRUE)
prev_scale <- ggplot2::scale_fill_viridis_c(
  option = "rocket", direction = -1, limits = PREV_LIMITS,
  labels = scales::percent, name = "stunting\nprevalence")

# ===========================================================================
# 1. Headline map
# ===========================================================================
p_main <- ggplot2::ggplot(com) +
  ggplot2::geom_sf(ggplot2::aes(fill = est), colour = NA) +
  region_lines + prev_scale +
  ggplot2::labs(
    title = "Child stunting prevalence by commune, Madagascar 2021",
    subtitle = paste0("Bayesian geostatistical small-area estimates, ",
                      nrow(com), " communes | DHS 2021 (EDSMD-V) with ",
                      "earth-observation covariates")) +
  theme_map
save_fig(p_main, "10_commune_stunting.png", width = 7, height = 9)

# ===========================================================================
# 2. Prevalence, uncertainty and burden side by side
# ===========================================================================
p_unc <- ggplot2::ggplot(com) +
  ggplot2::geom_sf(ggplot2::aes(fill = upper - lower), colour = NA) +
  region_lines +
  ggplot2::scale_fill_viridis_c(option = "mako", direction = -1,
                                labels = scales::percent,
                                name = "95% CI\nwidth") +
  ggplot2::labs(title = "Uncertainty",
                subtitle = "Width of the 95% credible interval") + theme_map

p_burden <- ggplot2::ggplot(com) +
  ggplot2::geom_sf(ggplot2::aes(fill = log10(pmax(n_stunted, 1))), colour = NA) +
  region_lines +
  ggplot2::scale_fill_viridis_c(option = "inferno", direction = -1,
                                name = "log10\nchildren") +
  ggplot2::labs(title = "Burden",
                subtitle = "Number of stunted children under 5") + theme_map

p_prev <- ggplot2::ggplot(com) +
  ggplot2::geom_sf(ggplot2::aes(fill = est), colour = NA) +
  region_lines + prev_scale +
  ggplot2::labs(title = "Prevalence",
                subtitle = "Share of children under 5 with HAZ < -2") + theme_map

panel <- p_prev + p_unc + p_burden + patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Prevalence, uncertainty and burden are three different maps",
    subtitle = paste("High prevalence and high burden do not coincide:",
                     "targeting on one is not targeting on the other."),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13)))
save_fig(panel, "10_prevalence_uncertainty_burden.png", width = 15, height = 8)

# ===========================================================================
# 3. The paradox panel
# ===========================================================================
f_cov <- OUT$cov_commune
if (file.exists(f_cov)) {
  cc <- utils::read.csv(f_cov)
  cm <- dplyr::left_join(com, cc[, c("ADM3_PCODE", "elevation")], by = "ADM3_PCODE")

  p_elev <- ggplot2::ggplot(cm) +
    ggplot2::geom_sf(ggplot2::aes(fill = elevation), colour = NA) +
    region_lines +
    ggplot2::scale_fill_viridis_c(option = "cividis", name = "metres") +
    ggplot2::labs(title = "Elevation",
                  subtitle = "The central highland plateau") + theme_map

  p_par <- ggplot2::ggplot(cm) +
    ggplot2::geom_sf(ggplot2::aes(fill = est), colour = NA) +
    region_lines + prev_scale +
    ggplot2::labs(title = "Stunting",
                  subtitle = "Highest where the land is most productive") + theme_map

  sc <- ggplot2::ggplot(sf::st_drop_geometry(cm),
                        ggplot2::aes(elevation, est)) +
    ggplot2::geom_point(ggplot2::aes(size = pop_u5), alpha = 0.25,
                        colour = "grey20") +
    ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                         colour = "firebrick") +
    ggplot2::geom_vline(xintercept = 800, linetype = "dashed") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_size_continuous(guide = "none") +
    ggplot2::labs(x = "Commune mean elevation (m)", y = "Stunting prevalence",
                  title = "The gradient",
                  subtitle = "Dashed line: 800 m, edge of the Hauts Plateaux") +
    ggplot2::theme_minimal(base_size = 9)

  par_panel <- p_elev + p_par + sc + patchwork::plot_layout(nrow = 1) +
    patchwork::plot_annotation(
      title = "The fertile-highland paradox",
      subtitle = paste("Madagascar's most agriculturally productive region",
                       "carries its highest burden of child stunting."),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13)))
  save_fig(par_panel, "10_paradox_panel.png", width = 15, height = 8)
}

# ===========================================================================
# 4. Comparison of estimators at district level
# ===========================================================================
f_fh <- file.path(DIR$processed, "fh_adm2.csv")
f_a2 <- file.path(DIR$processed, "aggregated_adm2.csv")
if (file.exists(f_fh) && file.exists(f_a2)) {
  adm2 <- sf::st_read(OUT$adm2, quiet = TRUE)
  fh <- utils::read.csv(f_fh)
  a2 <- utils::read.csv(f_a2)

  cmp <- adm2 |>
    dplyr::left_join(fh, by = "ADM2_PCODE") |>
    dplyr::left_join(a2[, c("ADM2_PCODE", "est")], by = "ADM2_PCODE")

  lims <- range(c(cmp$direct, cmp$fh_eo_est, cmp$summer_est, cmp$est), na.rm = TRUE)
  one <- function(var, ttl, sub) {
    ggplot2::ggplot(cmp) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data[[var]]), colour = "grey40",
                       linewidth = 0.08) +
      ggplot2::scale_fill_viridis_c(option = "rocket", direction = -1,
                                    limits = lims, labels = scales::percent,
                                    name = "stunting") +
      ggplot2::labs(title = ttl, subtitle = sub) + theme_map
  }
  cmp_panel <- one("direct", "Direct", "Design-based, DHS 2021") +
    one("fh_eo_est", "Fay-Herriot", "Area-level, non-spatial") +
    one("summer_est", "SUMMER BYM2", "Area-level, spatial") +
    one("est", "Geostatistical", "SPDE surface, aggregated") +
    patchwork::plot_layout(nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = "Four estimators of district stunting",
      subtitle = paste("The direct estimator is the noisiest;",
                       "the modelled surfaces agree on the spatial pattern."),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13)))
  save_fig(cmp_panel, "10_estimator_comparison.png", width = 18, height = 7)
}

# ===========================================================================
# 5. Ranked burden chart
# ===========================================================================
f_rank <- file.path(DIR$tables, "08_commune_burden_ranking.csv")
if (file.exists(f_rank)) {
  rk <- utils::read.csv(f_rank)
  top <- utils::head(rk, 30)
  top$label <- paste0(top$ADM3_EN, " (", top$ADM2_EN, ")")
  p_rank <- ggplot2::ggplot(top, ggplot2::aes(stats::reorder(label, n_stunted),
                                              n_stunted)) +
    ggplot2::geom_col(ggplot2::aes(fill = est)) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_c(option = "rocket", direction = -1,
                                  labels = scales::percent, name = "prevalence") +
    ggplot2::labs(x = NULL, y = "Stunted children under 5",
                  title = "The 30 communes carrying the largest burden",
                  subtitle = "Colour shows prevalence: the biggest burdens are not the highest rates") +
    ggplot2::theme_minimal(base_size = 9)
  save_fig(p_rank, "10_top_burden_communes.png", width = 8.5, height = 8)
}

msg("10_maps.R complete")
