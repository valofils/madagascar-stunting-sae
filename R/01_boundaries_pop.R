# ---------------------------------------------------------------------------
# 01_boundaries_pop.R
#
# Harmonise the official BNGRC/OCHA administrative boundaries (adm1 region,
# adm2 district, adm3 commune), attach WorldPop population (total and under-5)
# to every commune, and build the commune adjacency graph used later by the
# BYM2 / ICAR spatial models.
#
# Inputs : data/raw/boundaries/03_Boundaries_2025/*.shp
#          WorldPop 100 m constrained rasters (downloaded on first run)
# Outputs: data/processed/adm{1,2,3}_*.gpkg
#          data/processed/commune_population.csv
#          data/processed/adm3_nb.rds, adm3_adjacency.graph
#          outputs/figures/01_*.png, outputs/tables/01_*.csv
#
# Requires no DHS data - safe to execute now.
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("sf", "terra", "exactextractr", "dplyr", "ggplot2", "spdep")

sf::sf_use_s2(TRUE)

# ===========================================================================
# 1. Read and harmonise the administrative layers
# ===========================================================================
read_adm <- function(path, label) {
  stopifnot(file.exists(path))
  x <- sf::st_read(path, quiet = TRUE)
  # adm2 ships in World Mercator, the others in geographic coordinates.
  x <- sf::st_transform(x, CRS_GEO)
  # Self-intersections are common in these merged fokontany polygons and they
  # make s2 predicates (adjacency, point-in-polygon) fail, so repair up front.
  bad <- !sf::st_is_valid(x)
  if (any(bad, na.rm = TRUE)) {
    msg(label, ": repairing ", sum(bad, na.rm = TRUE), " invalid geometries")
    x <- sf::st_make_valid(x)
  }
  msg(label, ": ", nrow(x), " features")
  x
}

adm1 <- read_adm(BND$adm1, "adm1 region")
adm2 <- read_adm(BND$adm2, "adm2 district")
adm3 <- read_adm(BND$adm3, "adm3 commune")

# Keep only the identifying hierarchy. The PAM programme columns that ride
# along in these files (CARI_ECMEN, PECMAM, ...) are not part of this analysis.
adm1 <- adm1[, c("ADM1_PCODE", "ADM1_EN", "PROV_CODE", "OLD_PROVIN")]
adm2 <- adm2[, c("ADM1_PCODE", "ADM1_EN", "ADM2_PCODE", "ADM2_EN")]
adm3 <- adm3[, c("ADM1_PCODE", "ADM1_EN", "ADM2_PCODE", "ADM2_EN",
                 "ADM3_PCODE", "ADM3_EN")]

# ===========================================================================
# 1b. Repair the pcode hierarchy
# ===========================================================================
# The three layers were edited at different dates and disagree in three
# specific, reproducible ways. Each is fixed explicitly below rather than
# silently dropped, because all of them touch real inhabited areas.
#
#   (i)   Antananarivo's six arrondissements carry an 11-character ADM2_PCODE
#         in adm3 (MG11101001A) but a 8-character one in adm2 (MG11101A..F).
#   (ii)  The 2021 split of Vatovavy-Fitovinany is half-applied: adm2/adm3 code
#         Fitovinany as MG23A, while the adm1 layer codes it MG35. Two communes
#         of Manakara Atsimo also still carry the pre-split region code MG23
#         although their region name already reads Fitovinany.
#   (iii) Region Ambatosoa (MG34, split from Analanjirofo) exists in adm1 but
#         has no districts or communes in adm2/adm3.
#
# Rule adopted: adm2 is authoritative for the commune -> district -> region
# chain (it is the layer the communes were actually dissolved from), and
# region NAMES are used to reconcile adm1's codes with it.

# Digits MUST be kept. Antananarivo's arrondissements are named "1er
# Arrondissement", "2e Arrondissement", ... and stripping digits collapses 2e
# through 6e to the same string, so match() silently maps five distinct
# districts onto whichever one comes first. Keep alphanumerics only.
norm_nm <- function(x) gsub("[^A-Z0-9]", "", toupper(trimws(as.character(x))))

## (i) arrondissement crosswalk, matched on name inside Antananarivo city ----
arr2 <- sf::st_drop_geometry(adm2)[grepl("^MG11101", adm2$ADM2_PCODE), ]
is_arr <- grepl("^MG11101", adm3$ADM2_PCODE) & !(adm3$ADM2_PCODE %in% adm2$ADM2_PCODE)
if (any(is_arr)) {
  new_code <- arr2$ADM2_PCODE[match(norm_nm(adm3$ADM2_EN[is_arr]), norm_nm(arr2$ADM2_EN))]
  if (anyNA(new_code))
    stop("Arrondissement crosswalk failed: unmatched names.", call. = FALSE)
  # Each arrondissement must receive its OWN code. A collapse here would merge
  # districts of the capital without any error being raised, so assert it.
  if (anyDuplicated(new_code))
    stop("Arrondissement crosswalk collapsed ", sum(duplicated(new_code)),
         " districts onto a shared code - check norm_nm() is not stripping digits.",
         call. = FALSE)
  msg("repaired ", sum(is_arr), " communes: arrondissement ADM2_PCODE re-coded ",
      "to the adm2 form (e.g. MG11101001A -> MG11101A)")
  adm3$ADM2_PCODE[is_arr] <- new_code
}

## (ii) inherit the region code from the district, not from adm3's own column
adm2_lookup <- stats::setNames(adm2$ADM1_PCODE, adm2$ADM2_PCODE)
adm1_from_adm2 <- unname(adm2_lookup[adm3$ADM2_PCODE])
n_wrong <- sum(adm3$ADM1_PCODE != adm1_from_adm2, na.rm = TRUE)
if (n_wrong > 0) {
  bad_rows <- which(adm3$ADM1_PCODE != adm1_from_adm2)
  utils::write.csv(sf::st_drop_geometry(adm3)[bad_rows,
                     c("ADM3_PCODE", "ADM3_EN", "ADM2_EN", "ADM1_PCODE")],
                   file.path(DIR$tables, "01_repaired_region_codes.csv"),
                   row.names = FALSE)
  msg("repaired ", n_wrong, " communes whose ADM1_PCODE contradicted their district")
}
adm3$ADM1_PCODE <- adm1_from_adm2

## (ii cont.) reconcile adm2/adm3 region codes with the adm1 layer, by name --
# Fitovinany is MG23A in adm2/adm3 and MG35 in adm1; we adopt the adm1 code so
# that every layer in data/processed shares one region vocabulary.
adm1_key <- stats::setNames(adm1$ADM1_PCODE, norm_nm(adm1$ADM1_EN))
recode_region <- function(x, nm) {
  target <- unname(adm1_key[norm_nm(nm)])
  changed <- !is.na(target) & target != x
  if (any(changed)) {
    tab <- unique(data.frame(from = x[changed], to = target[changed],
                             name = nm[changed]))
    for (i in seq_len(nrow(tab)))
      msg("harmonised region code ", tab$from[i], " -> ", tab$to[i],
          " (", tab$name[i], ") to match the adm1 layer")
    x[changed] <- target[changed]
  }
  x
}
adm2$ADM1_PCODE <- recode_region(adm2$ADM1_PCODE, adm2$ADM1_EN)
adm3$ADM1_PCODE <- recode_region(adm3$ADM1_PCODE, adm3$ADM1_EN)

# Keep ADM1_EN/ADM2_EN consistent with the repaired codes.
adm3$ADM1_EN <- unname(stats::setNames(adm1$ADM1_EN, adm1$ADM1_PCODE)[adm3$ADM1_PCODE])
adm3$ADM2_EN <- unname(stats::setNames(adm2$ADM2_EN, adm2$ADM2_PCODE)[adm3$ADM2_PCODE])

## (iii) report regions that no commune belongs to ---------------------------
empty_regions <- adm1[!adm1$ADM1_PCODE %in% adm3$ADM1_PCODE, ]
if (nrow(empty_regions) > 0) {
  msg("NOTE: ", nrow(empty_regions), " adm1 region(s) contain no commune in the ",
      "2025 adm3 layer: ", paste(empty_regions$ADM1_EN, collapse = ", "))
  msg("      (their territory is still covered, under the parent region they ",
      "were split from - no population is lost, but results cannot be ",
      "reported for them separately)")
  utils::write.csv(sf::st_drop_geometry(empty_regions)[, c("ADM1_PCODE", "ADM1_EN")],
                   file.path(DIR$tables, "01_regions_without_communes.csv"),
                   row.names = FALSE)
}

# ---- Integrity checks (must all pass after the repairs above) --------------
chk <- function(cond, what) {
  if (!isTRUE(cond)) warning("CHECK FAILED: ", what, call. = FALSE) else msg("ok: ", what)
  invisible(cond)
}
chk(!anyDuplicated(adm3$ADM3_PCODE), "commune pcodes are unique")
chk(!anyDuplicated(adm2$ADM2_PCODE), "district pcodes are unique")
chk(all(adm3$ADM2_PCODE %in% adm2$ADM2_PCODE), "every commune's district exists in adm2")
chk(all(adm2$ADM1_PCODE %in% adm1$ADM1_PCODE), "every district's region exists in adm1")
chk(all(adm3$ADM1_PCODE %in% adm1$ADM1_PCODE), "every commune's region exists in adm1")
chk(!anyNA(adm3$ADM1_PCODE) && !anyNA(adm3$ADM2_PCODE), "no missing parent codes")

# Commune names are NOT unique nationally (many repeated Ambohimanga,
# Antanimena, ...), so ADM3_PCODE is the only safe join key. Record how bad it is.
msg("note: ", sum(duplicated(adm3$ADM3_EN)),
    " commune names are non-unique -> always join on ADM3_PCODE")

# ===========================================================================
# 2. Geometry-derived attributes
# ===========================================================================
adm3_eq <- sf::st_transform(adm3, CRS_EQA)
adm3$area_km2 <- as.numeric(sf::st_area(adm3_eq)) / 1e6

# Point on surface rather than centroid: guarantees the point falls inside the
# polygon even for crescent-shaped coastal communes.
pos <- sf::st_point_on_surface(sf::st_geometry(adm3_eq))
pos <- sf::st_coordinates(sf::st_transform(pos, CRS_GEO))
adm3$lon <- pos[, "X"]
adm3$lat <- pos[, "Y"]

# ===========================================================================
# 3. WorldPop population, aggregated to commune
# ===========================================================================
# Constrained 2020 release: built-settlement-constrained, so population is not
# smeared over uninhabited land - important for the sparse arid south-west.
WP_POP <- "https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/MDG"
WP_AGE <- "https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_Constrained/2020/MDG"
wp_dir <- file.path(DIR$rasters, "worldpop")

get_wp <- function(base, file) download_if_missing(paste0(base, "/", file),
                                                   file.path(wp_dir, file))

f_pop <- get_wp(WP_POP, "mdg_ppp_2020_UNadj_constrained.tif")
f_age <- vapply(c("mdg_f_0_2020_constrained.tif", "mdg_f_1_2020_constrained.tif",
                  "mdg_m_0_2020_constrained.tif", "mdg_m_1_2020_constrained.tif"),
                function(f) get_wp(WP_AGE, f), character(1))

# WorldPop age bands: "_0_" = 0 to <1 year, "_1_" = 1 to <5 years.
# Under-5 = f_0 + f_1 + m_0 + m_1, the DHS anthropometry denominator.
msg("summing WorldPop under-5 bands")
u5 <- sum(terra::rast(unname(f_age)), na.rm = TRUE)
names(u5) <- "pop_u5"

msg("zonal statistics over ", nrow(adm3), " communes (a few minutes)")
adm3$pop_total <- exactextractr::exact_extract(terra::rast(f_pop), adm3, "sum",
                                               progress = FALSE)
adm3$pop_u5 <- exactextractr::exact_extract(u5, adm3, "sum", progress = FALSE)

adm3$pop_dens <- adm3$pop_total / adm3$area_km2
adm3$u5_share <- adm3$pop_u5 / adm3$pop_total

msg("national population  : ", format(round(sum(adm3$pop_total, na.rm = TRUE)), big.mark = " "))
msg("national under-5     : ", format(round(sum(adm3$pop_u5, na.rm = TRUE)), big.mark = " "))
msg("communes with zero population: ",
    sum(adm3$pop_total <= 0 | is.na(adm3$pop_total)))

# ===========================================================================
# 4. Commune adjacency graph (for BYM2 / ICAR)
# ===========================================================================
# queen = TRUE: communes touching at a single point count as neighbours, which
# keeps the graph better connected along Madagascar's ragged coastline.
msg("building adjacency graph")
nb <- spdep::poly2nb(adm3, row.names = adm3$ADM3_PCODE, queen = TRUE)

n_iso <- sum(spdep::card(nb) == 0)
nc_before <- spdep::n.comp.nb(nb)$nc
msg("neighbours: mean ", round(mean(spdep::card(nb)), 2),
    " | isolated: ", n_iso, " | components: ", nc_before)

# Offshore islands (Nosy Be, Sainte-Marie, the Barren isles) form their own
# components. An ICAR/BYM2 prior on a disconnected graph is improper: each
# component gets its own free level, so the model cannot borrow strength into
# the islands and the intercept is only weakly identified. Rather than relying
# on INLA's constraint bookkeeping, we make the graph connected by joining each
# component to the nearest polygon of another component, and record every edge
# we invented so the choice is auditable.
cent <- sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(adm3_eq)))

add_edge <- function(nb, i, j) {
  nb[[i]] <- sort(unique(c(as.integer(nb[[i]][nb[[i]] > 0L]), j)))
  nb[[j]] <- sort(unique(c(as.integer(nb[[j]][nb[[j]] > 0L]), i)))
  nb
}

added <- list()
repeat {
  comp <- spdep::n.comp.nb(nb)
  if (comp$nc <= 1) break
  # Grow the largest component by absorbing whichever component is closest.
  sizes <- table(comp$comp.id)
  main <- as.integer(names(sizes)[which.max(sizes)])
  in_main <- which(comp$comp.id == main)
  outside <- which(comp$comp.id != main)

  d <- outer(seq_along(outside), seq_along(in_main), function(a, b) {
    ia <- outside[a]; ib <- in_main[b]
    sqrt((cent[ia, 1] - cent[ib, 1])^2 + (cent[ia, 2] - cent[ib, 2])^2)
  })
  k <- which(d == min(d), arr.ind = TRUE)[1, ]
  i <- outside[k[1]]; j <- in_main[k[2]]
  nb <- add_edge(nb, i, j)
  added[[length(added) + 1]] <- data.frame(
    from_pcode = adm3$ADM3_PCODE[i], from_name = adm3$ADM3_EN[i],
    to_pcode = adm3$ADM3_PCODE[j], to_name = adm3$ADM3_EN[j],
    distance_km = round(min(d) / 1000, 1))
}

if (length(added) > 0) {
  added <- do.call(rbind, added)
  utils::write.csv(added, file.path(DIR$tables, "01_adjacency_edges_added.csv"),
                   row.names = FALSE)
  msg("added ", nrow(added), " artificial edge(s) to connect island communes ",
      "(see outputs/tables/01_adjacency_edges_added.csv)")
}
attr(nb, "region.id") <- adm3$ADM3_PCODE

nc_after <- spdep::n.comp.nb(nb)$nc
chk(nc_after == 1, "commune adjacency graph is connected (required by ICAR)")

saveRDS(nb, OUT$adj_nb)
spdep::nb2INLA(OUT$adj_graph, nb)
msg("adjacency written: ", nc_before, " components before patching, ",
    nc_after, " after")

# ===========================================================================
# 5. Write outputs
# ===========================================================================
sf::st_write(adm1, OUT$adm1, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(adm2, OUT$adm2, delete_dsn = TRUE, quiet = TRUE)
sf::st_write(adm3, OUT$adm3, delete_dsn = TRUE, quiet = TRUE)

utils::write.csv(sf::st_drop_geometry(adm3), OUT$commune_pop, row.names = FALSE)
msg("wrote ", basename(OUT$commune_pop))

# District-level population, used to population-weight commune predictions up
# to the DHS sampling domains in 08_aggregate.R.
adm2_pop <- sf::st_drop_geometry(adm3) |>
  dplyr::group_by(ADM1_PCODE, ADM1_EN, ADM2_PCODE, ADM2_EN) |>
  dplyr::summarise(pop_total = sum(pop_total, na.rm = TRUE),
                   pop_u5 = sum(pop_u5, na.rm = TRUE),
                   n_communes = dplyr::n(), .groups = "drop")
utils::write.csv(adm2_pop, file.path(DIR$tables, "01_district_population.csv"),
                 row.names = FALSE)

# ===========================================================================
# 6. Diagnostic maps
# ===========================================================================
theme_map <- ggplot2::theme_void(base_size = 9) +
  ggplot2::theme(legend.position = "right",
                 plot.title = ggplot2::element_text(face = "bold", size = 11))

p1 <- ggplot2::ggplot(adm3) +
  ggplot2::geom_sf(ggplot2::aes(fill = log10(pmax(pop_u5, 1))), colour = NA) +
  ggplot2::geom_sf(data = adm1, fill = NA, colour = "white", linewidth = 0.2) +
  ggplot2::scale_fill_viridis_c(option = "magma", name = "log10\nunder-5") +
  ggplot2::labs(title = "Under-5 population by commune",
                subtitle = "WorldPop 2020, 100 m constrained, UN-adjusted") +
  theme_map
save_fig(p1, "01_pop_u5_commune.png")

p2 <- ggplot2::ggplot(adm3) +
  ggplot2::geom_sf(ggplot2::aes(fill = log10(pmax(pop_dens, 0.1))), colour = NA) +
  ggplot2::scale_fill_viridis_c(option = "viridis", name = "log10\nper km2") +
  ggplot2::labs(title = "Population density by commune",
                subtitle = "The dense central highlands stand out as the bright core") +
  theme_map
save_fig(p2, "01_pop_density_commune.png")

msg("01_boundaries_pop.R complete")
