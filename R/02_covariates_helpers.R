# ---------------------------------------------------------------------------
# 02_covariates_helpers.R
#
# Raster acquisition, one function per source. Sourcing this file has no side
# effects beyond creating cache directories: nothing is downloaded until a
# getter is actually called. Both 02_covariates.R (commune and cluster
# extraction) and 05_model_spde.R (prediction grid) source this file, so the
# model frame and the prediction frame are guaranteed to come from identical
# rasters - the single most common source of silent error in this kind of
# pipeline.
#
# Every getter returns a terra SpatRaster cropped to Madagascar, with layer
# names already set to the covariate names used in the models.
# ---------------------------------------------------------------------------

if (!exists("DIR")) source("R/00_setup.R")
need("terra")

MDG_EXT <- terra::ext(42.5, 51.5, -26.5, -11.5)

rdir <- function(sub) {
  d <- file.path(DIR$rasters, sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
crop_mdg <- function(r) terra::crop(r, MDG_EXT, snap = "out")

# --- H1: elevation and terrain ---------------------------------------------
get_elevation <- function() {
  f <- file.path(rdir("elevation"), "mdg_elevation_30s.tif")
  if (!file.exists(f)) {
    need("geodata")
    msg("downloading SRTM elevation")
    e <- geodata::elevation_30s(country = "MDG", path = rdir("elevation"))
    terra::writeRaster(crop_mdg(e), f, overwrite = TRUE)
  }
  r <- terra::rast(f)
  names(r) <- "elevation"
  # Terrain Ruggedness Index: mean absolute elevation difference to the eight
  # neighbours. Separates "high and flat" (the highland plateau, where the
  # paradox lives) from "high and broken" (the eastern escarpment).
  tri <- terra::terrain(r, v = "TRI", neighbors = 8)
  names(tri) <- "ruggedness"
  c(r, tri)
}

# --- H1/H2: WorldClim bioclimatic variables --------------------------------
get_worldclim <- function() {
  f <- file.path(rdir("worldclim"), "mdg_bioclim.tif")
  if (!file.exists(f)) {
    need("geodata")
    msg("downloading WorldClim bioclim")
    # worldclim_country is always 30 arc-seconds; there is no res argument.
    b <- crop_mdg(geodata::worldclim_country("MDG", var = "bio",
                                             path = rdir("worldclim")))
    # bio1  annual mean temperature
    # bio4  temperature seasonality        (highland annual swing)
    # bio6  min temperature, coldest month (H1 cold stress: the key one)
    # bio12 annual precipitation
    # bio15 precipitation seasonality      (hunger-season exposure)
    b <- b[[c(1, 4, 6, 12, 15)]]
    names(b) <- c("temp_mean", "temp_seasonality", "temp_min_cold",
                  "precip_annual", "precip_seasonality")
    terra::writeRaster(b, f, overwrite = TRUE)
  }
  terra::rast(f)
}

# --- H2: livestock density (FAO GLW4 cattle) -------------------------------
# Proxy for access to animal-source foods, the nutrient missing from a
# rice-dominated highland diet.
get_livestock <- function() {
  f <- file.path(rdir("livestock"), "mdg_cattle.tif")
  if (!file.exists(f)) {
    src <- file.path(rdir("livestock"), "GLW4_cattle_global.tif")
    if (!file.exists(src))
      stop("Cattle density raster missing.\n",
           "  Download GLW4 cattle (2020, dasymetric) from\n",
           "  https://data.apps.fao.org/catalog/  and save it as:\n  ", src,
           call. = FALSE)
    terra::writeRaster(crop_mdg(terra::rast(src)), f, overwrite = TRUE)
  }
  r <- terra::rast(f); names(r) <- "cattle_density"; r
}

# --- H4: travel time to cities (Weiss et al. 2018) -------------------------
get_accessibility <- function() {
  f <- file.path(rdir("access"), "mdg_travel_time_cities.tif")
  if (!file.exists(f))
    stop("Travel-time raster missing.\n",
         "  Download 'Accessibility to Cities 2015' from https://data.malariaatlas.org/,\n",
         "  crop to Madagascar and save as:\n  ", f, call. = FALSE)
  r <- terra::rast(f); names(r) <- "travel_time"; r
}

# --- H4: night-time lights (VIIRS annual composite 2021) -------------------
get_nightlights <- function() {
  f <- file.path(rdir("viirs"), "mdg_viirs_2021.tif")
  if (!file.exists(f))
    stop("VIIRS night-lights raster missing.\n",
         "  Export NOAA/VIIRS/DNB/ANNUAL_V22 (band 'average', 2021) for\n",
         "  Madagascar from Earth Engine and save as:\n  ", f, call. = FALSE)
  r <- terra::rast(f); names(r) <- "nightlights"; r
}

# --- H3: WorldPop population density ---------------------------------------
get_popdens <- function() {
  f <- file.path(rdir("worldpop"), "mdg_ppp_2020_UNadj_constrained.tif")
  if (!file.exists(f)) stop("Run 01_boundaries_pop.R first (WorldPop cache).",
                            call. = FALSE)
  r <- terra::rast(f); names(r) <- "pop_count"; r
}

# Continuous rasters, extracted identically for communes, clusters and grid.
RASTER_GETTERS <- list(
  elevation = get_elevation,
  worldclim = get_worldclim,
  livestock = get_livestock,
  access    = get_accessibility,
  viirs     = get_nightlights,
  popdens   = get_popdens
)

# --- H2/H3: ESA WorldCover 10 m land-cover class fractions -----------------
WORLDCOVER_CLASSES <- c(tree = 10, shrub = 20, grass = 30, crop = 40,
                        built = 50, bare = 60, water_perm = 80, wetland = 90,
                        mangrove = 95)

get_worldcover_tiles <- function() {
  base <- paste0("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/",
                 "map/ESA_WorldCover_10m_2021_v200_")
  tiles <- expand.grid(lat = seq(-27, -12, by = 3), lon = seq(42, 51, by = 3))
  ids <- unique(sprintf("S%02dE%03d", abs(tiles$lat), tiles$lon))
  files <- character(0)
  for (id in ids) {
    f <- file.path(rdir("worldcover"), paste0(id, "_Map.tif"))
    # Ocean-only tiles simply do not exist in the archive; a 404 is expected.
    ok <- tryCatch({ download_if_missing(paste0(base, id, "_Map.tif"), f,
                                         quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (ok) files <- c(files, f)
  }
  msg("WorldCover tiles available: ", length(files), " of ", length(ids))
  files
}

get_worldcover_vrt <- function() {
  tiles <- get_worldcover_tiles()
  if (length(tiles) == 0) return(NULL)
  terra::vrt(tiles, file.path(rdir("worldcover"), "mdg_worldcover.vrt"),
             overwrite = TRUE)
}
