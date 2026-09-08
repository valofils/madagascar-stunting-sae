# ---------------------------------------------------------------------------
# 03b_dhs_region_crosswalk.R
#
# The DHS reporting regions are NOT the same set as the 2025 adm1 layer, so
# every comparison between a survey estimate and a modelled aggregate needs an
# explicit crosswalk. Two genuine mismatches, both affecting a large share of
# the sample:
#
#   1. DHS splits Analamanga in two: v024 = 10 "antananarivo" (the capital) and
#      v024 = 11 "analamanga" (the rest of the region). The adm1 layer has only
#      Analamanga (MG11). The split is a DISTRICT boundary - the capital is the
#      six arrondissements of Antananarivo Renivohitra (ADM2_PCODE MG11101A-F) -
#      so the crosswalk has to live at adm2, not adm1.
#   2. DHS keeps the pre-2021 Vatovavy-Fitovinany as one region (v024 = 23),
#      while the adm1 layer applies the split into Vatovavy (MG23) and
#      Fitovinany (MG35).
#
# Together those two cases cover roughly a fifth of the children, so name
# matching alone would silently drop them.
#
# Method. The crosswalk is derived EMPIRICALLY rather than hand-written: the
# GPS displacement for this survey is restricted to admin2 (stated in
# MDGE81FL/GPS_Displacement_README.txt), which means a displaced cluster is
# guaranteed to fall in its true district. Cluster -> district is therefore
# exact, and cross-tabulating district against the DHS region recovers the
# mapping from the data itself. Districts holding no cluster are then filled by
# region-name matching, and the whole thing is checked for being a clean
# partition.
#
# Inputs : data/raw/dhs/MDGE81FL/*.shp, MDKR81DT/*.DTA, data/processed/adm2*.gpkg
# Outputs: data/processed/dhs_region_crosswalk.csv   (ADM2_PCODE -> DHS region)
#          outputs/tables/03b_*.csv
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("sf", "dplyr", "haven")

# ===========================================================================
# 1. Inputs
# ===========================================================================
ge_file <- list.files(DIR$dhs, pattern = "GE.*[.]shp$", full.names = TRUE,
                      recursive = TRUE, ignore.case = TRUE)[1]
kr_file <- list.files(DIR$dhs, pattern = "^MDKR.*[.]DTA$", full.names = TRUE,
                      recursive = TRUE, ignore.case = TRUE)[1]
stopifnot(!is.na(ge_file), !is.na(kr_file))

ge <- sf::st_read(ge_file, quiet = TRUE)
ge <- ge[!(ge$LATNUM == 0 & ge$LONGNUM == 0), ]
adm2 <- sf::st_read(OUT$adm2, quiet = TRUE)
adm3 <- sf::st_read(OUT$adm3, quiet = TRUE)

kr <- haven::read_dta(kr_file, col_select = c("v001", "v024"))
reg_lab <- attr(kr$v024, "labels")
cl_region <- kr |>
  dplyr::transmute(DHSCLUST = as.numeric(v001), dhs_region = as.numeric(v024)) |>
  dplyr::distinct()

# A cluster must belong to exactly one DHS region.
dupes <- cl_region$DHSCLUST[duplicated(cl_region$DHSCLUST)]
if (length(dupes) > 0)
  warning("clusters with more than one v024 value: ", paste(dupes, collapse = ", "))

msg("clusters with GPS: ", nrow(ge), " | clusters in KR: ", nrow(cl_region))

# ===========================================================================
# 2. Cluster -> district (exact, because displacement is admin2-restricted)
# ===========================================================================
j <- sf::st_join(ge["DHSCLUST"], adm2[, c("ADM1_PCODE", "ADM2_PCODE")],
                 join = sf::st_intersects, left = TRUE)
j <- sf::st_drop_geometry(j)

# Coastal displacement can still push a point just offshore of every polygon;
# snap those to the nearest district rather than losing the cluster.
miss <- which(is.na(j$ADM2_PCODE))
if (length(miss) > 0) {
  nearest <- sf::st_nearest_feature(ge[miss, ], adm2)
  j[miss, c("ADM1_PCODE", "ADM2_PCODE")] <-
    sf::st_drop_geometry(adm2)[nearest, c("ADM1_PCODE", "ADM2_PCODE")]
  msg("snapped ", length(miss), " clusters that fell outside every district")
}

clu <- dplyr::inner_join(j, cl_region, by = "DHSCLUST")
clu$dhs_region_name <- names(reg_lab)[match(clu$dhs_region, reg_lab)]
msg("clusters matched to a district and a DHS region: ", nrow(clu))

# ===========================================================================
# 3. Derive the district -> DHS region map
# ===========================================================================
tab <- clu |>
  dplyr::count(ADM2_PCODE, dhs_region, dhs_region_name, name = "n_clusters")

# A district assigned to more than one DHS region means either the adm2 layer
# disagrees with the DHS admin2, or displacement crossed a district boundary
# after all. Either way it must be seen, not averaged away.
conflict <- tab |>
  dplyr::count(ADM2_PCODE, name = "n_regions") |>
  dplyr::filter(n_regions > 1)

if (nrow(conflict) > 0) {
  det <- tab |> dplyr::filter(ADM2_PCODE %in% conflict$ADM2_PCODE) |>
    dplyr::left_join(sf::st_drop_geometry(adm2)[, c("ADM2_PCODE", "ADM2_EN")],
                     by = "ADM2_PCODE") |>
    dplyr::arrange(ADM2_PCODE, dplyr::desc(n_clusters))
  msg("WARNING: ", nrow(conflict), " district(s) map to more than one DHS region; ",
      "taking the majority and logging the detail")
  print(as.data.frame(det), row.names = FALSE)
  utils::write.csv(det, file.path(DIR$tables, "03b_district_region_conflicts.csv"),
                   row.names = FALSE)
}

observed <- tab |>
  dplyr::group_by(ADM2_PCODE) |>
  dplyr::slice_max(n_clusters, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(ADM2_PCODE, dhs_region, dhs_region_name)

msg("districts resolved from data: ", nrow(observed), " of ", nrow(adm2))

# ===========================================================================
# 4. Fill the districts that hold no cluster
# ===========================================================================
# These get their DHS region from the region name of their adm1 parent. The two
# known irregular cases are handled explicitly rather than left to the matcher.
norm <- function(x) gsub("[^A-Z]", "", toupper(trimws(as.character(x))))

adm1_of <- stats::setNames(adm2$ADM1_PCODE, adm2$ADM2_PCODE)
adm1_nm <- sf::st_drop_geometry(sf::st_read(OUT$adm1, quiet = TRUE))
adm1_name_of <- stats::setNames(adm1_nm$ADM1_EN, adm1_nm$ADM1_PCODE)

by_name <- stats::setNames(as.numeric(reg_lab), norm(names(reg_lab)))
# Vatovavy and Fitovinany are one DHS region (23); the adm1 layer splits them.
by_name[["VATOVAVY"]] <- 23
by_name[["FITOVINANY"]] <- 23
# Ambatosoa was split out of Analanjirofo after the survey frame was drawn.
by_name[["AMBATOSOA"]] <- unname(by_name[["ANALANJIROFO"]])

full <- data.frame(ADM2_PCODE = adm2$ADM2_PCODE,
                   ADM2_EN = adm2$ADM2_EN,
                   ADM1_PCODE = unname(adm1_of[adm2$ADM2_PCODE]),
                   stringsAsFactors = FALSE)
full$ADM1_EN <- unname(adm1_name_of[full$ADM1_PCODE])
full <- dplyr::left_join(full, observed, by = "ADM2_PCODE")

fill <- is.na(full$dhs_region)
full$dhs_region[fill] <- unname(by_name[norm(full$ADM1_EN[fill])])

# The six Antananarivo arrondissements are the DHS "antananarivo" region (10);
# everything else in Analamanga is DHS "analamanga" (11). Applied only where the
# data did not already settle it.
ANTANANARIVO_CITY <- grepl("^MG11101", full$ADM2_PCODE)
still_na <- is.na(full$dhs_region)
full$dhs_region[ANTANANARIVO_CITY & still_na] <- 10
full$dhs_region[full$ADM1_PCODE == "MG11" & !ANTANANARIVO_CITY & still_na] <- 11

full$dhs_region_name <- names(reg_lab)[match(full$dhs_region, reg_lab)]
full$source <- ifelse(!is.na(full$dhs_region) & full$ADM2_PCODE %in% observed$ADM2_PCODE,
                      "observed (has clusters)", "filled by region name")

# ===========================================================================
# 5. Checks
# ===========================================================================
chk <- function(cond, what) {
  if (!isTRUE(cond)) warning("CHECK FAILED: ", what, call. = FALSE) else msg("ok: ", what)
}
chk(!anyNA(full$dhs_region), "every district has a DHS region")
chk(nrow(full) == nrow(adm2), "one row per district")
chk(all(sort(unique(clu$dhs_region)) %in% sort(unique(full$dhs_region))),
    "every sampled DHS region appears in the crosswalk")

unmapped <- full[is.na(full$dhs_region), ]
if (nrow(unmapped) > 0) {
  msg("districts with no DHS region:")
  print(unmapped[, c("ADM2_PCODE", "ADM2_EN", "ADM1_EN")], row.names = FALSE)
}

utils::write.csv(full, file.path(DIR$processed, "dhs_region_crosswalk.csv"),
                 row.names = FALSE)

# Commune-level version, for population-weighting model output to DHS regions.
com_x <- sf::st_drop_geometry(adm3)[, c("ADM3_PCODE", "ADM2_PCODE")] |>
  dplyr::left_join(full[, c("ADM2_PCODE", "dhs_region", "dhs_region_name")],
                   by = "ADM2_PCODE")
utils::write.csv(com_x, file.path(DIR$processed, "dhs_region_crosswalk_commune.csv"),
                 row.names = FALSE)

summ <- full |>
  dplyr::count(dhs_region, dhs_region_name, name = "n_districts") |>
  dplyr::arrange(dhs_region)
utils::write.csv(summ, file.path(DIR$tables, "03b_dhs_region_summary.csv"),
                 row.names = FALSE)
msg("crosswalk written: ", nrow(full), " districts -> ", nrow(summ), " DHS regions")
print(as.data.frame(summ), row.names = FALSE)

msg("03b_dhs_region_crosswalk.R complete")
