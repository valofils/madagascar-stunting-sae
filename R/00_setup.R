# ---------------------------------------------------------------------------
# 00_setup.R  --  Paths, packages, constants, small helpers.
# Sourced at the top of every other script. Safe to source repeatedly.
# ---------------------------------------------------------------------------

## --- Project root ----------------------------------------------------------
# Scripts are expected to run with the working directory at the project root
# (or anywhere below it). We walk up until we see a marker that identifies it.
find_root <- function(start = getwd()) {
  markers <- c("R/00_setup.R", "CLAUDE.md")
  p <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (any(file.exists(file.path(p, markers)))) return(p)
    parent <- dirname(p)
    if (identical(parent, p))
      stop("Project root not found above ", start,
           " (looked for: ", paste(markers, collapse = ", "), ")")
    p <- parent
  }
}
PROJ_ROOT <- find_root()

path_ <- function(...) file.path(PROJ_ROOT, ...)

DIR <- list(
  raw        = path_("data", "raw"),
  boundaries = path_("data", "raw", "boundaries"),
  dhs        = path_("data", "raw", "dhs"),
  rasters    = path_("data", "raw", "rasters"),
  interim    = path_("data", "interim"),
  processed  = path_("data", "processed"),
  figures    = path_("outputs", "figures"),
  tables     = path_("outputs", "tables")
)
invisible(lapply(DIR, dir.create, recursive = TRUE, showWarnings = FALSE))

## --- Packages --------------------------------------------------------------
# need(): attach a package, failing with an actionable message rather than a
# bare "there is no package called ...".
need <- function(...) {
  for (p in c(...)) {
    ok <- suppressPackageStartupMessages(
      requireNamespace(p, quietly = TRUE) && require(p, character.only = TRUE, quietly = TRUE)
    )
    if (!ok) {
      stop("Package '", p, "' is required but not installed.\n",
           if (p == "INLA")
             "  install.packages('INLA', repos = c(INLA = 'https://inla.r-inla-download.org/R/stable'))"
           else paste0("  install.packages('", p, "')"),
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

## --- Coordinate reference systems -----------------------------------------
CRS_GEO <- "EPSG:4326"   # storage / DHS GPS / raster extraction
# Lambert azimuthal equal-area centred on Madagascar: use for areas, distances,
# adjacency and any metric buffer. Madagascar spans ~12S-26S, 43E-51E.
CRS_EQA <- "+proj=laea +lat_0=-19 +lon_0=47 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

## --- Admin layer file names (source: BNGRC/OCHA via PAM, 2025 edition) ------
BND <- list(
  adm1 = file.path(DIR$boundaries, "03_Boundaries_2025", "mdg_bnd_adm1_region_pam_2025.shp"),
  adm2 = file.path(DIR$boundaries, "03_Boundaries_2025", "mdg_bnd_adm2_district_pam_2025.shp"),
  adm3 = file.path(DIR$boundaries, "03_Boundaries_2025", "mdg_bnd_adm3_com_pam_2025.shp"),
  adm4 = file.path(DIR$boundaries, "03_Boundaries_2025", "mdg_bnd_adm4_bngrc_ocha_20181031.shp")
)

## --- Processed-output file names (single source of truth for all scripts) ---
OUT <- list(
  adm1        = file.path(DIR$processed, "adm1_regions.gpkg"),
  adm2        = file.path(DIR$processed, "adm2_districts.gpkg"),
  adm3        = file.path(DIR$processed, "adm3_communes.gpkg"),
  commune_pop = file.path(DIR$processed, "commune_population.csv"),
  adj_graph   = file.path(DIR$processed, "adm3_adjacency.graph"),  # INLA/SUMMER format
  adj_nb      = file.path(DIR$processed, "adm3_nb.rds"),
  cov_commune = file.path(DIR$processed, "commune_covariates.csv"),
  cov_cluster = file.path(DIR$processed, "cluster_covariates.csv"),
  dhs_child   = file.path(DIR$interim,   "dhs_child_haz.rds"),
  direct_adm2 = file.path(DIR$processed, "direct_estimates_adm2.csv"),
  direct_adm1 = file.path(DIR$processed, "direct_estimates_adm1.csv")
)

## --- Helpers ---------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# Download only if absent (or zero-length). Returns the destination path.
download_if_missing <- function(url, dest, quiet = FALSE) {
  if (file.exists(dest) && file.size(dest) > 0) {
    if (!quiet) msg("cached: ", basename(dest))
    return(invisible(dest))
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!quiet) msg("downloading: ", basename(dest))
  tmp <- paste0(dest, ".part")
  old <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  utils::download.file(url, tmp, mode = "wb", quiet = quiet, method = "libcurl")
  if (!file.exists(tmp) || file.size(tmp) == 0) stop("Download produced an empty file: ", url)
  file.rename(tmp, dest)
  invisible(dest)
}

# Save a ggplot with consistent defaults.
save_fig <- function(plot, name, width = 7, height = 8, dpi = 300) {
  f <- file.path(DIR$figures, name)
  ggplot2::ggsave(f, plot, width = width, height = height, dpi = dpi, bg = "white")
  msg("figure written: ", name)
  invisible(f)
}

options(stringsAsFactors = FALSE)
