# ---------------------------------------------------------------------------
# 03_dhs_direct.R
#
# Read the 2021 Madagascar DHS (EDSMD-V) recodes, build the child-level
# analysis file (HAZ, stunting, covariates), and produce design-based DIRECT
# estimates of stunting prevalence with their sampling variances.
#
# These direct estimates are the input to every model that follows:
#   - 04 SUMMER smoothed-direct (BYM2) consumes the logit-scale direct
#     estimates and their design variances,
#   - 06 Fay-Herriot consumes the same,
#   - 09 benchmarks the modelled results back to them.
#
# Inputs : data/raw/dhs/  MDKR*.DTA (children), MDPR*.DTA (household members),
#                         MDIR*.DTA (women), MDGE*.shp (GPS)
# Outputs: data/interim/dhs_child_haz.rds
#          data/processed/direct_estimates_adm{1,2}.csv
#          outputs/tables/03_*.csv, outputs/figures/03_*.png
#
# DHS microdata are NOT redistributable: data/raw/dhs is gitignored.
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("haven", "dplyr", "tidyr", "rlang", "survey", "srvyr", "sf", "ggplot2")

options(survey.lonely.psu = "adjust")   # single-PSU strata: centre at the mean

# ===========================================================================
# 1. Locate the recode files
# ===========================================================================
find_dta <- function(pattern, what) {
  f <- list.files(DIR$dhs, pattern = pattern, full.names = TRUE,
                  recursive = TRUE, ignore.case = TRUE)
  if (length(f) == 0)
    stop("DHS ", what, " recode not found under ", DIR$dhs, "\n",
         "  expected a file matching: ", pattern, call. = FALSE)
  if (length(f) > 1) msg("note: several ", what, " files, using ", basename(f[1]))
  f[1]
}

kr_file <- find_dta("^MDKR.*[.]DTA$", "children (KR)")
msg("reading ", basename(kr_file))
kr <- haven::read_dta(kr_file)
msg("KR records: ", nrow(kr))

# ===========================================================================
# 2. Build the child-level analysis file
# ===========================================================================
# DHS variable map (standard recode VII):
#   v001 cluster    v002 household   v005 women's weight (x 1e-6)
#   v022 sample stratum   v023 sample domain   v024 region   v025 urban/rural
#   b5   child alive      b8/hw1 age in months  b4 sex
#   hw70 height-for-age Z, WHO 2006 standard, stored x 100
#   hw71 weight-for-age, hw72 weight-for-height
#   v106 mother's education   v190 wealth index quintile
#   v437 mother's weight (x10 kg)  v438 mother's height (x10 cm)
#   v212 age of mother at first birth   bord birth order

pick <- function(df, ...) {
  cand <- c(...)
  hit <- cand[cand %in% names(df)]
  if (length(hit) == 0) return(rep(NA_real_, nrow(df)))
  as.numeric(df[[hit[1]]])
}

child <- dplyr::tibble(
  cluster   = pick(kr, "v001"),
  household = pick(kr, "v002"),
  wt        = pick(kr, "v005") / 1e6,
  strata    = pick(kr, "v022", "v023"),
  region    = pick(kr, "v024"),
  urban     = as.integer(pick(kr, "v025") == 1),
  alive     = pick(kr, "b5"),
  # hw1 is the anthropometry age in months and is the one hw70 was computed
  # from; b19 (or b8*12) can disagree by a month for children measured late.
  age_month = pick(kr, "hw1", "b19"),
  sex       = pick(kr, "b4"),
  haz_raw   = pick(kr, "hw70"),
  waz_raw   = pick(kr, "hw71"),
  whz_raw   = pick(kr, "hw72"),
  mother_edu   = pick(kr, "v106"),
  wealth_q     = pick(kr, "v190"),
  mother_bmi   = pick(kr, "v445") / 100,
  mother_age1b = pick(kr, "v212"),
  birth_order  = pick(kr, "bord"),
  bcg          = pick(kr, "h2"),
  diarrhea_2w  = pick(kr, "h11")
)

# HAZ is stored multiplied by 100. Values 9996-9999 are DHS flags
# (996 height out of range, 997 inconsistent, 998 not measured, 999 missing),
# and WHO flags |HAZ| > 6 as biologically implausible.
child <- child |>
  dplyr::mutate(
    haz = dplyr::if_else(haz_raw > 9000 | is.na(haz_raw), NA_real_, haz_raw / 100),
    waz = dplyr::if_else(waz_raw > 9000 | is.na(waz_raw), NA_real_, waz_raw / 100),
    whz = dplyr::if_else(whz_raw > 9000 | is.na(whz_raw), NA_real_, whz_raw / 100),
    haz = dplyr::if_else(abs(haz) > 6, NA_real_, haz),          # WHO flag
    stunted = as.integer(haz < -2),
    stunted_severe = as.integer(haz < -3)
  )

n_all <- nrow(child)
child <- child |>
  dplyr::filter(alive == 1, !is.na(age_month), age_month < 60, !is.na(haz))

msg("children 0-59 months with a valid HAZ: ", nrow(child), " of ", n_all,
    " KR records (", round(100 * nrow(child) / n_all, 1), "%)")

# ---- Join cluster geography ------------------------------------------------
ge_files <- list.files(DIR$dhs, pattern = "GE.*[.]shp$", full.names = TRUE,
                       recursive = TRUE, ignore.case = TRUE)
if (length(ge_files) == 0) {
  msg("WARNING: no GE shapefile - no commune/district linkage, ",
      "direct estimates limited to DHS region (v024).")
  child$ADM1_PCODE <- NA_character_
  child$ADM2_PCODE <- NA_character_
  child$ADM3_PCODE <- NA_character_
} else {
  # 02_covariates.R already resolved every cluster to its commune (including
  # snapping the ones displaced offshore); reuse that rather than redoing it.
  if (file.exists(OUT$cov_cluster)) {
    clu <- utils::read.csv(OUT$cov_cluster)
    child <- dplyr::left_join(
      child, clu[, c("DHSCLUST", "ADM1_PCODE", "ADM2_PCODE", "ADM3_PCODE")],
      by = c("cluster" = "DHSCLUST"))
  } else {
    stop("Run 02_covariates.R first: it resolves DHS clusters to communes.",
         call. = FALSE)
  }
  msg("children linked to a commune: ", sum(!is.na(child$ADM3_PCODE)),
      " of ", nrow(child))
}

saveRDS(child, OUT$dhs_child)
msg("wrote ", basename(OUT$dhs_child))

# ===========================================================================
# 3. Survey design
# ===========================================================================
# Two-stage stratified design: strata = region x urban/rural (v022), PSU =
# cluster (v001), weight = v005/1e6. Ignoring this understates the variance by
# roughly the design effect, which for DHS anthropometry is typically 1.5-2.5.
des <- srvyr::as_survey_design(child, ids = cluster, strata = strata,
                               weights = wt, nest = TRUE)

deff_national <- survey::svymean(~stunted, des, deff = TRUE)
msg("national stunting: ", round(100 * coef(deff_national), 1), "% (SE ",
    round(100 * survey::SE(deff_national), 2), ", DEFF ",
    round(survey::deff(deff_national), 2), ")")

# ===========================================================================
# 4. Direct estimates by area
# ===========================================================================
# For every area we need, on top of the prevalence:
#   - the design-based variance (Fay-Herriot / SUMMER input),
#   - the same on the logit scale, where the normal approximation behaves far
#     better for proportions near 0 or 1 (Mercer et al. 2015).
direct_by <- function(design, area_var) {
  a <- rlang::sym(area_var)
  est <- design |>
    srvyr::group_by(!!a) |>
    srvyr::summarise(
      n_children = srvyr::unweighted(dplyr::n()),
      n_clusters = srvyr::unweighted(dplyr::n_distinct(cluster)),
      direct     = srvyr::survey_mean(stunted, vartype = "var", na.rm = TRUE),
      mean_haz   = srvyr::survey_mean(haz, vartype = "var", na.rm = TRUE)
    ) |>
    dplyr::rename(haz_var = mean_haz_var)

  # Delta-method transfer of the variance to the logit scale:
  #   Var(logit p) = Var(p) / (p (1 - p))^2
  est |>
    dplyr::mutate(
      direct = pmin(pmax(direct, 1e-6), 1 - 1e-6),
      logit_direct = qlogis(direct),
      logit_var = direct_var / (direct * (1 - direct))^2,
      se = sqrt(direct_var),
      cv = se / direct,
      ci_low = pmax(0, direct - 1.96 * se),
      ci_high = pmin(1, direct + 1.96 * se)
    )
}

## --- Region (DHS domain: the design supports these directly) ---------------
direct_adm1 <- direct_by(des, "region")

# v024 is a labelled Stata integer; recover its labels and match them to the
# official ADM1 pcodes so downstream joins are on pcode, not on region name.
region_lab <- attr(kr$v024, "labels")
if (!is.null(region_lab)) {
  direct_adm1$region_name <- names(region_lab)[match(direct_adm1$region, region_lab)]
}
adm1 <- sf::st_read(OUT$adm1, quiet = TRUE)
norm_nm <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z]", "", x)
  x
}
direct_adm1$ADM1_PCODE <- adm1$ADM1_PCODE[match(norm_nm(direct_adm1$region_name),
                                                norm_nm(adm1$ADM1_EN))]
unmatched <- direct_adm1$region_name[is.na(direct_adm1$ADM1_PCODE)]
if (length(unmatched) > 0)
  msg("WARNING: unmatched DHS region names -> fix by hand: ",
      paste(unmatched, collapse = ", "))

utils::write.csv(direct_adm1, OUT$direct_adm1, row.names = FALSE)

## --- District (adm2) -------------------------------------------------------
# The 2021 DHS was NOT powered at district level, so many districts hold only a
# handful of clusters. These estimates are deliberately noisy: they exist to be
# smoothed by the models in 04/06, and their instability is the whole point of
# doing small-area estimation.
if (any(!is.na(child$ADM2_PCODE))) {
  des2 <- srvyr::as_survey_design(dplyr::filter(child, !is.na(ADM2_PCODE)),
                                  ids = cluster, strata = strata,
                                  weights = wt, nest = TRUE)
  direct_adm2 <- direct_by(des2, "ADM2_PCODE")
  utils::write.csv(direct_adm2, OUT$direct_adm2, row.names = FALSE)

  msg("districts with a direct estimate: ", nrow(direct_adm2),
      " | median clusters/district: ", stats::median(direct_adm2$n_clusters),
      " | districts with <3 clusters: ", sum(direct_adm2$n_clusters < 3))
  msg("median CV of the district direct estimates: ",
      round(stats::median(direct_adm2$cv, na.rm = TRUE), 3),
      "  (a CV above 0.30 is normally considered unpublishable)")
}

## --- Commune (adm3): coverage diagnostic, not an estimate ------------------
# ~1,700 communes vs ~600 DHS clusters: most communes contain no cluster at
# all. This table quantifies exactly that gap, which is the motivation for the
# whole modelling exercise.
if (any(!is.na(child$ADM3_PCODE))) {
  pop <- utils::read.csv(OUT$commune_pop)
  cover <- child |>
    dplyr::filter(!is.na(ADM3_PCODE)) |>
    dplyr::group_by(ADM3_PCODE) |>
    dplyr::summarise(n_children = dplyr::n(),
                     n_clusters = dplyr::n_distinct(cluster), .groups = "drop")
  cover <- dplyr::left_join(pop[, c("ADM3_PCODE", "ADM3_EN", "pop_u5")],
                            cover, by = "ADM3_PCODE") |>
    dplyr::mutate(n_children = tidyr::replace_na(n_children, 0L),
                  n_clusters = tidyr::replace_na(n_clusters, 0L))
  utils::write.csv(cover, file.path(DIR$tables, "03_commune_dhs_coverage.csv"),
                   row.names = FALSE)
  msg("communes containing at least one DHS cluster: ",
      sum(cover$n_clusters > 0), " of ", nrow(cover),
      " (", round(100 * mean(cover$n_clusters > 0), 1), "%)")
}

# ===========================================================================
# 5. First look at the paradox
# ===========================================================================
# Mean HAZ against elevation, at cluster level. If H1 has any purchase, this
# should already slope downward before any modelling.
if (file.exists(OUT$cov_cluster)) {
  clu <- utils::read.csv(OUT$cov_cluster)
  cl_haz <- child |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(mean_haz = mean(haz, na.rm = TRUE),
                     stunting = mean(stunted, na.rm = TRUE),
                     n = dplyr::n(), .groups = "drop") |>
    dplyr::left_join(clu, by = c("cluster" = "DHSCLUST"))

  if ("elevation" %in% names(cl_haz)) {
    p <- ggplot2::ggplot(cl_haz, ggplot2::aes(elevation, mean_haz)) +
      ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.35) +
      ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                           colour = "firebrick") +
      ggplot2::geom_hline(yintercept = -2, linetype = "dashed") +
      ggplot2::labs(
        x = "Cluster mean elevation (m)", y = "Cluster mean HAZ", size = "children",
        title = "The fertile-highland paradox, before any modelling",
        subtitle = "DHS 2021 clusters: mean height-for-age against elevation") +
      ggplot2::theme_minimal(base_size = 10)
    save_fig(p, "03_haz_vs_elevation_clusters.png", width = 7, height = 5)
  }
  utils::write.csv(cl_haz, file.path(DIR$tables, "03_cluster_haz.csv"),
                   row.names = FALSE)
}

msg("03_dhs_direct.R complete")
