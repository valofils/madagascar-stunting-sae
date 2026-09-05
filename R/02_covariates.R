# ---------------------------------------------------------------------------
# 02_covariates.R
#
# Assemble the earth-observation covariate stack that encodes the four paradox
# hypotheses, then summarise it twice:
#   (a) per commune  -> prediction frame for the area-level models (runs NOW)
#   (b) per DHS cluster -> model frame for the geostatistical model
#       (needs data/raw/dhs/*GE*.shp; skipped cleanly until that lands)
#
# Hypothesis -> covariate mapping
#   H1 altitude / cold stress .. elevation, ruggedness, min temperature of the
#                                coldest month, temperature seasonality
#   H2 rice monoculture ........ cropland fraction, precipitation, permanent
#                                water fraction, cattle density (protein access)
#   H3 infection load .......... population density, built-up fraction,
#                                permanent-water fraction (paddy proxy)
#   H4 care / access ........... travel time to cities, night-time lights
#
# Raster acquisition lives in 02_covariates_helpers.R so that this script and
# 05_model_spde.R build their covariates from exactly the same sources.
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("sf", "terra", "exactextractr", "dplyr")
source("R/02_covariates_helpers.R")

terra::terraOptions(progress = 0)

stopifnot(file.exists(OUT$adm3))
adm3 <- sf::st_read(OUT$adm3, quiet = TRUE)

# ===========================================================================
# 1. Generic extraction
# ===========================================================================
# `zones` is any sf polygon layer: communes, or buffered DHS cluster points.
# Both go through this same function so the model frame and the prediction
# frame are constructed identically.
extract_continuous <- function(zones, id_col) {
  out <- sf::st_drop_geometry(zones)[, id_col, drop = FALSE]
  for (nm in names(RASTER_GETTERS)) {
    r <- tryCatch(RASTER_GETTERS[[nm]](), error = function(e) {
      msg("SKIPPED ", nm, ": ", conditionMessage(e)); NULL
    })
    if (is.null(r)) next
    msg("extracting ", nm, " (", terra::nlyr(r), " layer(s))")
    z <- sf::st_transform(zones, terra::crs(r))
    vals <- as.data.frame(exactextractr::exact_extract(r, z, "mean",
                                                       progress = FALSE))
    names(vals) <- names(r)
    out <- cbind(out, vals)
  }
  out
}

# Coarse rasters have no value at all for zones smaller than one cell that sit
# in water-masked areas - GLW4 cattle is 5 arc-minutes (~10 km), so the small
# offshore communes (Nosy Komba, Nosy Be) come back NA. Left alone, a single
# NA commune propagates: 04 and 06 drop any covariate that is missing anywhere,
# so one 26 km2 island would remove cattle density from the whole analysis.
# Fill from the nearest zone that does have a value, and say how many.
fill_from_nearest <- function(df, zones, vars) {
  stopifnot(nrow(df) == nrow(zones))
  for (v in intersect(vars, names(df))) {
    miss <- which(is.na(df[[v]]))
    if (length(miss) == 0) next
    have <- which(!is.na(df[[v]]))
    if (length(have) == 0) { msg("all values missing for ", v, " - left as NA"); next }
    nearest <- sf::st_nearest_feature(zones[miss, ], zones[have, ])
    df[[v]][miss] <- df[[v]][have[nearest]]
    msg("filled ", length(miss), " missing ", v, " from the nearest zone")
  }
  df
}

extract_landcover <- function(zones, id_col) {
  vrt <- tryCatch(get_worldcover_vrt(), error = function(e) NULL)
  if (is.null(vrt)) {
    msg("SKIPPED land cover: no WorldCover tiles available")
    return(NULL)
  }
  msg("extracting WorldCover class fractions (slow: 10 m grid)")
  z <- sf::st_transform(zones, terra::crs(vrt))
  frac <- exactextractr::exact_extract(
    vrt, z,
    function(value, coverage_fraction) {
      tot <- sum(coverage_fraction, na.rm = TRUE)
      vapply(WORLDCOVER_CLASSES,
             function(k) sum(coverage_fraction[value == k], na.rm = TRUE) / tot,
             numeric(1))
    },
    progress = FALSE, summarize_df = FALSE)
  frac <- as.data.frame(t(simplify2array(frac)))
  names(frac) <- paste0("frac_", names(WORLDCOVER_CLASSES))
  cbind(sf::st_drop_geometry(zones)[, id_col, drop = FALSE], frac)
}

# ===========================================================================
# 2. Commune-level covariates  (runs now, no DHS needed)
# ===========================================================================
msg("=== commune-level extraction (", nrow(adm3), " communes) ===")
cov_com <- extract_continuous(adm3, "ADM3_PCODE")
cov_com <- fill_from_nearest(cov_com, adm3, setdiff(names(cov_com), "ADM3_PCODE"))
lc_com <- extract_landcover(adm3, "ADM3_PCODE")
if (!is.null(lc_com)) cov_com <- dplyr::left_join(cov_com, lc_com, by = "ADM3_PCODE")

# Population variables from 01 (H3 infection load).
pop <- utils::read.csv(OUT$commune_pop)
cov_com <- dplyr::left_join(
  cov_com,
  pop[, c("ADM3_PCODE", "ADM1_PCODE", "ADM2_PCODE", "ADM3_EN",
          "area_km2", "lon", "lat", "pop_total", "pop_u5", "pop_dens")],
  by = "ADM3_PCODE")

utils::write.csv(cov_com, OUT$cov_commune, row.names = FALSE)
msg("wrote ", basename(OUT$cov_commune), ": ", nrow(cov_com), " communes x ",
    ncol(cov_com), " columns")

# Report coverage so a silently-missing raster cannot slip into the models.
cover <- data.frame(
  variable = names(cov_com),
  pct_missing = round(100 * colMeans(is.na(cov_com)), 2),
  row.names = NULL)
utils::write.csv(cover, file.path(DIR$tables, "02_covariate_coverage.csv"),
                 row.names = FALSE)
gaps <- cover[cover$pct_missing > 0, ]
if (nrow(gaps) > 0) {
  msg("covariates with missing values:")
  print(gaps, row.names = FALSE)
}

# ===========================================================================
# 3. DHS cluster covariates  (waits for the GE shapefile)
# ===========================================================================
# DHS GPS coordinates are displaced: up to 2 km urban, 5 km rural, with 1% of
# rural clusters displaced up to 10 km. Reading a raster at the exact point
# would therefore sample the wrong place. Standard practice (Perez-Heydrich et
# al. 2013) is to average over a buffer matching the displacement radius.
CLUSTER_BUFFER_M <- c(urban = 2000, rural = 5000)

ge_files <- list.files(DIR$dhs, pattern = "GE.*[.]shp$", full.names = TRUE,
                       recursive = TRUE, ignore.case = TRUE)

if (length(ge_files) == 0) {
  msg("No DHS GE shapefile under ", DIR$dhs, " - cluster extraction skipped.")
  msg("Drop the MDGE*FL folder there and re-run this script.")
} else {
  ge <- sf::st_read(ge_files[1], quiet = TRUE)
  msg("=== DHS cluster extraction (", nrow(ge), " clusters) ===")

  # DHS records an unknown position as exactly (0, 0); those cannot be used.
  bad <- ge$LATNUM == 0 & ge$LONGNUM == 0
  if (any(bad)) msg("dropping ", sum(bad), " clusters with missing coordinates")
  ge <- ge[!bad, ]

  ge_eq <- sf::st_transform(ge, CRS_EQA)
  ge_eq$buffer_m <- ifelse(toupper(ge_eq$URBAN_RURA) == "U",
                           CLUSTER_BUFFER_M[["urban"]], CLUSTER_BUFFER_M[["rural"]])
  buf <- sf::st_buffer(ge_eq, dist = ge_eq$buffer_m)

  cov_clu <- extract_continuous(buf, "DHSCLUST")
  cov_clu <- fill_from_nearest(cov_clu, buf, setdiff(names(cov_clu), "DHSCLUST"))
  lc_clu <- extract_landcover(buf, "DHSCLUST")
  if (!is.null(lc_clu)) cov_clu <- dplyr::left_join(cov_clu, lc_clu, by = "DHSCLUST")

  pts <- sf::st_transform(ge, CRS_GEO)
  xy <- sf::st_coordinates(pts)
  cov_clu$lon <- xy[, "X"]
  cov_clu$lat <- xy[, "Y"]
  cov_clu$urban <- as.integer(toupper(ge$URBAN_RURA) == "U")

  # Resolve each cluster to its commune once, here, so 03 and 05 can just join.
  j <- sf::st_join(pts["DHSCLUST"],
                   adm3[, c("ADM1_PCODE", "ADM2_PCODE", "ADM3_PCODE")],
                   join = sf::st_intersects, left = TRUE)
  cov_clu <- dplyr::left_join(cov_clu, sf::st_drop_geometry(j), by = "DHSCLUST")

  # Displacement can push a coastal cluster offshore; snap those to the nearest
  # commune rather than losing them.
  miss <- which(is.na(cov_clu$ADM3_PCODE))
  if (length(miss) > 0) {
    msg("snapping ", length(miss), " clusters displaced outside every commune")
    nearest <- sf::st_nearest_feature(pts[miss, ], adm3)
    cov_clu[miss, c("ADM1_PCODE", "ADM2_PCODE", "ADM3_PCODE")] <-
      sf::st_drop_geometry(adm3)[nearest, c("ADM1_PCODE", "ADM2_PCODE", "ADM3_PCODE")]
  }

  utils::write.csv(cov_clu, OUT$cov_cluster, row.names = FALSE)
  msg("wrote ", basename(OUT$cov_cluster), ": ", nrow(cov_clu), " clusters")
}

msg("02_covariates.R complete")
