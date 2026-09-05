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

# --- H2: livestock density (Gridded Livestock of the World 4, cattle) ------
# Proxy for access to animal-source foods, the nutrient missing from a
# rice-dominated highland diet. GLW4 is published on Harvard Dataverse; we take
# the DASYMETRIC product (Da), which redistributes census counts using
# suitability covariates, rather than the areal-weighted (Aw) one, because Aw
# simply spreads district totals uniformly and would wash out exactly the
# within-district contrast this project is about.
# Resolution is 5 arc-minutes (~10 km), coarser than the other covariates: it
# resolves district-scale variation but not commune-scale variation, and the
# effect table should not over-read a fine-grained cattle signal.
GLW4_CATTLE_URL <- "https://dataverse.harvard.edu/api/access/datafile/6769711"

get_livestock <- function() {
  f <- file.path(rdir("livestock"), "mdg_cattle.tif")
  if (!file.exists(f)) {
    src <- download_if_missing(GLW4_CATTLE_URL,
                               file.path(rdir("livestock"), "GLW4_cattle_global.tif"))
    terra::writeRaster(crop_mdg(terra::rast(src)), f, overwrite = TRUE)
  }
  r <- terra::rast(f); names(r) <- "cattle_density"; r
}

# --- H4: travel time (Malaria Atlas Project WCS) ---------------------------
# Two complementary accessibility surfaces, both at 30 arc-seconds (~1 km):
#   travel_time            Weiss et al. (2018), minutes to the nearest city
#                          (>=50,000 inhabitants), nominal year 2015.
#   travel_time_healthcare Weiss et al. (2020), motorized minutes to the
#                          nearest health facility, nominal year 2020.
# The second is the more direct H4 measure for this project, since the commune
# is the catchment unit for Madagascar's basic health centres, and it is closer
# in time to the 2021 DHS. Both are kept: distance to markets and distance to
# care are different exposures and need not move together.
#
# The MAP GeoServer serves these over WCS 2.0.1 and honours a bounding-box
# subset, so we fetch only the Madagascar window rather than the global grid.
MAP_WCS <- "https://data.malariaatlas.org/geoserver/ows"

map_wcs_url <- function(coverage_id) {
  paste0(MAP_WCS, "?service=WCS&version=2.0.1&request=GetCoverage",
         "&coverageId=", coverage_id,
         "&format=image%2Fgeotiff",
         "&subset=Lat(-26.5,-11.5)&subset=Long(42.5,51.5)")
}

get_accessibility <- function() {
  f <- file.path(rdir("access"), "mdg_travel_time_cities.tif")
  download_if_missing(map_wcs_url("Accessibility__201501_Global_Travel_Time_to_Cities"), f)
  r <- terra::rast(f); names(r) <- "travel_time"; r
}

get_accessibility_health <- function() {
  f <- file.path(rdir("access"), "mdg_travel_time_healthcare.tif")
  download_if_missing(
    map_wcs_url("Accessibility__202001_Global_Motorized_Travel_Time_to_Healthcare"), f)
  r <- terra::rast(f); names(r) <- "travel_time_healthcare"; r
}

# --- H4: night-time lights, 2021 -------------------------------------------
# The obvious source (NOAA/EOG annual VNL V2) moved behind an OAuth account, so
# this uses the harmonized DMSP-VIIRS series of Li et al. (Sci Data 2020,
# extended to 2024), which is openly hosted on figshare and is VIIRS-derived
# from 2014 onward.
#
# Consequence worth knowing: values are harmonized DMSP-like digital numbers on
# a 0-63 scale, NOT VIIRS radiances in nW/cm2/sr. That scale saturates over
# bright urban cores, so nightlights here is a usable rural/peri-urban economic
# activity gradient but must not be read as a linear intensity measure in
# Antananarivo.
NTL_2021_URL <- "https://ndownloader.figshare.com/files/57065294"

get_nightlights <- function() {
  f <- file.path(rdir("viirs"), "mdg_nightlights_2021.tif")
  if (!file.exists(f)) {
    src <- download_if_missing(
      NTL_2021_URL,
      file.path(rdir("viirs"), "Harmonized_DN_NTL_2021_simVIIRS_global.tif"))
    terra::writeRaster(crop_mdg(terra::rast(src)), f, overwrite = TRUE)
  }
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
  elevation     = get_elevation,
  worldclim     = get_worldclim,
  livestock     = get_livestock,
  access        = get_accessibility,
  access_health = get_accessibility_health,
  viirs         = get_nightlights,
  popdens       = get_popdens
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
