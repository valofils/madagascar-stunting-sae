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
  diarrhea_2w  = pick(kr, "h11"),
  fever_2w     = pick(kr, "h22"),
  ari_2w       = pick(kr, "h31"),
  breastfeeding = pick(kr, "m4"),
  water_source = pick(kr, "v113"),
  toilet_type  = pick(kr, "v116"),
  toilet_share = pick(kr, "v160")
)

# ===========================================================================
# 2b. Constructed exposures for H2 (diet quality) and H3 (infection load)
# ===========================================================================
# The covariate stack cannot carry either hypothesis on its own. 02b showed the
# highlands have MORE cattle and BETTER access than the lowlands, so livestock
# density is functioning as a wealth proxy, not a diet measure; H2 has to be
# tested with what children actually ate. Likewise H3 needs a direct infection
# signal, not only settlement density.
#
# The v414 series is the most country-customised part of the DHS recode, so the
# mapping below was VERIFIED against the variable labels in MDKR81FL.DTA on
# 2026-09-08 rather than assumed from the generic recode manual. All eight WHO
# groups resolve, and the country-specific slots (v414a-d, v414t, v414u, v414w)
# are empty in Madagascar and correctly unused. Confirmed codings:
#   m4 == 95 is "still breastfeeding" (the value labels also carry 93/94/96/97/98)
#   v414* are 0/1/8 with 8 = "don't know", which is treated as "not given",
#     the DHS convention.
# One genuine limitation: Madagascar has "other fruits" (v414l) but no separate
# "other vegetables" item, so the eighth WHO group is fruit-only and MDD may be
# very slightly understated relative to surveys that carry both.
#
# The mapping is still written as CANDIDATE lists and the script still prints
# which variables it resolved, so a different survey or a re-release cannot
# silently change the definition underneath the results.

# WHO 2021 minimum dietary diversity, 8 food groups, children 6-23 months.
FOOD_GROUPS <- list(
  breastmilk        = character(0),                 # handled separately via m4
  grains_roots      = c("v414e", "v414f"),
  legumes_nuts      = c("v414o"),
  dairy             = c("v411", "v411a", "v414p", "v414v"),
  flesh_foods       = c("v414h", "v414m", "v414n"),
  eggs              = c("v414g"),
  vitA_fruit_veg    = c("v414i", "v414j", "v414k"),
  other_fruit_veg   = c("v414l")
)

# DHS codes these as 0/1 with 8 = "don't know"; anything not a clean 1 is
# treated as "not given" rather than missing, which is the DHS convention.
ate <- function(df, vars) {
  vars <- vars[vars %in% names(df)]
  if (length(vars) == 0) return(list(value = rep(NA_integer_, nrow(df)), used = character(0)))
  m <- vapply(vars, function(v) as.integer(as.numeric(df[[v]]) == 1), integer(nrow(df)))
  list(value = as.integer(rowSums(m, na.rm = TRUE) > 0), used = vars)
}

fg_found <- list()
fg_mat <- matrix(0L, nrow = nrow(kr), ncol = 0)
for (g in names(FOOD_GROUPS)) {
  if (g == "breastmilk") {
    # m4 == 95 is "still breastfeeding" in the standard recode.
    val <- as.integer(pick(kr, "m4") == 95)
    val[is.na(val)] <- 0L
    fg_found[[g]] <- "m4 (==95)"
  } else {
    r <- ate(kr, FOOD_GROUPS[[g]])
    val <- r$value
    fg_found[[g]] <- if (length(r$used)) paste(r$used, collapse = "+") else "NONE FOUND"
    if (all(is.na(val))) val <- rep(NA_integer_, nrow(kr))
  }
  fg_mat <- cbind(fg_mat, val)
  colnames(fg_mat)[ncol(fg_mat)] <- g
}

n_missing_groups <- sum(vapply(fg_found, function(x) identical(x, "NONE FOUND"), logical(1)))
msg("dietary diversity food groups resolved:")
for (g in names(fg_found)) msg("    ", g, ": ", fg_found[[g]])
if (n_missing_groups > 0)
  msg("WARNING: ", n_missing_groups, " of 8 food groups had no matching variable. ",
      "MDD is NOT comparable to the published indicator - check the EDSMD-V ",
      "country-specific recode documentation before using it.")

child$diet_diversity <- as.integer(rowSums(fg_mat, na.rm = TRUE))
# MDD is defined only for children 6-23 months; outside that window the food
# questions are not asked consistently and the score is meaningless.
in_window <- child$age_month >= 6 & child$age_month <= 23
child$diet_diversity[!in_window | is.na(in_window)] <- NA_integer_
child$mdd <- as.integer(child$diet_diversity >= 5)   # WHO 2021 threshold: 5 of 8

utils::write.csv(
  data.frame(food_group = names(fg_found), variables_used = unlist(fg_found)),
  file.path(DIR$tables, "03_dietary_diversity_mapping.csv"), row.names = FALSE)

# ---- WASH (H3) ------------------------------------------------------------
# JMP (WHO/UNICEF Joint Monitoring Programme) service-ladder definitions, with
# the code lists checked against the v113/v116 value labels in this file on
# 2026-09-08. Codes 96 ("other") and 97 ("not a dejure resident") are NA, not 0.
#
# IMPROVED_WATER: piped (11-14, including piped to a neighbour), tube well or
#   borehole (21), protected dug well (31), protected spring (41), rainwater
#   (51), and delivered water - tanker (61), cart (62), bottled (71) - which the
#   JMP has counted as improved since 2017.
IMPROVED_WATER <- c(11, 12, 13, 14, 21, 31, 41, 51, 61, 62, 71)
# IMPROVED_TOILET: flush to sewer (11), septic tank (12) or pit (13); VIP (21);
#   pit latrine with slab (22); composting (41). Deliberately EXCLUDES 14
#   ("flush to somewhere else"), which the JMP classifies as unimproved, and 15
#   ("flush, don't know where"), which is ambiguous; both are negligible here
#   (9 and 6 records).
IMPROVED_TOILET <- c(11, 12, 13, 21, 22, 41)

child$improved_water <- as.integer(child$water_source %in% IMPROVED_WATER)
child$improved_water[is.na(child$water_source) | child$water_source >= 96] <- NA_integer_

# Two DISTINCT indicators, kept separate because they answer different
# questions and conflating them is the usual way this gets misreported:
#   improved_sanitation - the facility type is improved, sharing ignored
#   basic_sanitation    - improved AND not shared with other households, which
#                         is the JMP "at least basic" definition
# For H3 (faecal-oral exposure) the shared/unshared distinction matters, so
# basic_sanitation is the one the infection block uses.
child$improved_sanitation <- as.integer(child$toilet_type %in% IMPROVED_TOILET)
child$improved_sanitation[is.na(child$toilet_type) | child$toilet_type >= 96] <- NA_integer_
child$basic_sanitation <- as.integer(child$improved_sanitation == 1 &
                                       (is.na(child$toilet_share) | child$toilet_share != 1))
child$basic_sanitation[is.na(child$improved_sanitation)] <- NA_integer_
child$open_defecation <- as.integer(child$toilet_type %in% c(30, 31))
child$open_defecation[is.na(child$toilet_type) | child$toilet_type >= 96] <- NA_integer_

for (v in c("water_source", "toilet_type")) {
  tab <- as.data.frame(table(child[[v]], useNA = "ifany"))
  names(tab) <- c("code", "n")
  utils::write.csv(tab, file.path(DIR$tables, paste0("03_codes_", v, ".csv")),
                   row.names = FALSE)
}
msg("WASH (unweighted, children's households): improved water ",
    round(100 * mean(child$improved_water, na.rm = TRUE), 1),
    "% | improved sanitation ",
    round(100 * mean(child$improved_sanitation, na.rm = TRUE), 1),
    "% | basic sanitation ",
    round(100 * mean(child$basic_sanitation, na.rm = TRUE), 1),
    "% | open defecation ",
    round(100 * mean(child$open_defecation, na.rm = TRUE), 1), "%")

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

# DHS weights (v005/1e6) are normalised to sum to the SAMPLE size, not to the
# population, so deff = TRUE cannot form the simple-random-sample comparison and
# returns NA with a "sample size greater than population size" warning.
# deff = "replace" computes the design effect from the sample itself, which is
# the appropriate choice for self-weighting-within-stratum survey weights.
deff_national <- survey::svymean(~stunted, des, deff = "replace")
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

## --- DHS reporting region (the design domain) ------------------------------
# v024 IS the domain the survey was powered for, so these are the estimates the
# published EDSMD-V report contains and the ones the model is benchmarked to.
#
# It is deliberately NOT forced onto a single ADM1_PCODE. The DHS regions are
# not a relabelling of the 2025 adm1 layer: DHS splits Analamanga into
# "antananarivo" (the six arrondissements) and "analamanga" (the rest), and
# keeps the pre-2021 "vatovavy fitovinany" whole where the adm1 layer splits it
# into Vatovavy and Fitovinany. Roughly a fifth of the children sit in those two
# cases, so name matching would drop them. 03b_dhs_region_crosswalk.R resolves
# the relation at DISTRICT level - which is where the real boundary lies - and
# 08/09 aggregate the model over districts sharing a dhs_region to compare.
direct_adm1 <- direct_by(des, "region")
names(direct_adm1)[names(direct_adm1) == "region"] <- "dhs_region"

region_lab <- attr(kr$v024, "labels")
if (!is.null(region_lab))
  direct_adm1$dhs_region_name <- names(region_lab)[match(direct_adm1$dhs_region,
                                                         region_lab)]

xwalk_f <- file.path(DIR$processed, "dhs_region_crosswalk.csv")
if (file.exists(xwalk_f)) {
  xw <- utils::read.csv(xwalk_f)
  cover <- xw |>
    dplyr::group_by(dhs_region) |>
    dplyr::summarise(n_districts = dplyr::n(),
                     districts = paste(sort(ADM2_PCODE), collapse = ";"),
                     .groups = "drop")
  direct_adm1 <- dplyr::left_join(direct_adm1, cover, by = "dhs_region")
  missing_x <- direct_adm1$dhs_region_name[is.na(direct_adm1$n_districts)]
  if (length(missing_x) > 0)
    msg("WARNING: DHS regions absent from the crosswalk: ",
        paste(missing_x, collapse = ", "))
} else {
  msg("NOTE: run 03b_dhs_region_crosswalk.R to link DHS regions to districts.")
}

utils::write.csv(direct_adm1, OUT$direct_adm1, row.names = FALSE)
msg("direct estimates for ", nrow(direct_adm1), " DHS reporting regions")

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
