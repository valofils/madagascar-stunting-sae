# CLAUDE.md — Geospatial Small-Area Estimation of Child Stunting in Madagascar

## Objective

Produce commune-level estimates of child stunting (HAZ < −2 SD) in Madagascar from the 2021 DHS (EDSMD-V), using a spatially explicit, earth-observation–enriched small-area estimation (SAE) pipeline — and use the model to explain the **fertile-highland paradox** (stunting is highest in the agro-ecologically richest central highlands).

Working language of this repo: English (code, comments, docs). Analyst works FR/EN. Stack: R only for modelling (Stata cannot run SUMMER/INLA). Windows + VS Code.

## Research question & novelty

Why is stunting maximal in Madagascar's fertile central highlands, and can a spatial, EO-enriched model on the 2021 DHS both improve precision and explain this poverty–nutrition decoupling, at commune level?

Positioning against the two direct precedents (do not just reproduce them):

- **World Bank WP 10627** (Matekenya, Mulangu, Newhouse, 2023): commune-level SAE of poverty (EPM 2021-22) + stunting (DHS 2021) via unit-level EBP (`povmap`/`emdi`), borrowing from the 2018 census. Two documented weaknesses we exploit: (1) the stunting model is very weak (marginal R² ≈ 0.075, conditional R² ≈ 0.113) because it used only census-shared household variables and no environmental/geospatial covariates; (2) the authors explicitly could not explain the highland paradox.
- **IHME Local Burden of Disease** (Nature 2020): Bayesian model-based geostatistics, 5×5 km stunting surfaces for all of Africa/LMICs incl. Madagascar. But: data only up to ~2017-2019 (not the post-COVID 2021 DHS), continental covariate calibration (not Malagasy agro-ecology), aggregation to admin-2 (district, not commune), and purely predictive (no explanation of the paradox).

**Our contribution:** newest data (2021 DHS), Madagascar-specific agro-ecological covariates, commune granularity (the actual targeting unit for basic health centres), a spatial model, and — the core intellectual payoff — an interpretable decomposition of the highland paradox.

## Paradox hypotheses (drive the covariate set)

1. **Altitude / cold stress** — highlands at ~1,400–1,500 m: hypoxia, respiratory infection, higher energy expenditure (lower HAZ, cf. Andean/Ethiopian literature).
2. **Rice monoculture** — caloric sufficiency but low dietary diversity, few animal-source foods.
3. **Infection load** — settlement density + paddy water → environmental enteropathy.
4. **Care practices / early pregnancy** (contextual DHS variables).

## Data

### DHS 2021 Madagascar (EDSMD-V) — status: APPROVED 2026-09-08 (Survey + GPS)

Both approvals granted to valofils@gmail.com on 2026-09-08: Survey data (archive@dhsprogram.com) and GPS data (gpsrequests@dhsprogram.com). Project registered as *"Geospatial small-area estimation of child stunting in Madagascar: explaining the fertile-highland paradox using the 2021 DHS"*.

Download from <https://dhsprogram.com/data/dataset_admin/login_main.cfm> (Stata format), unzip, and place the extracted folders under `data/raw/dhs/`. Madagascar is country code `MD`, 2021 is DHS phase 8, so the files are:

| Recode | Zip | Extracted file the pipeline looks for | Used for |
|---|---|---|---|
| **KR** Children's | `MDKR81DT.ZIP` | `MDKR81FL.DTA` | HAZ, stunting, dietary diversity, WASH, maternal variables |
| **GE** GPS | `MDGE81FL.ZIP` | `MDGE81FL.shp` | cluster coordinates for covariate linkage and the SPDE model |
| **PR** Household member | `MDPR81DT.ZIP` | `MDPR81FL.DTA` | optional cross-check on household WASH |
| **IR** Women's | `MDIR81DT.ZIP` | `MDIR81FL.DTA` | optional; ANC and maternal detail not carried into KR |

KR and GE are the two that actually gate the pipeline. The scripts glob recursively and case-insensitively (`^MDKR.*[.]DTA$`, `GE.*[.]shp$`), so the exact folder layout under `data/raw/dhs/` does not matter.

**Terms of use.** Microdata must not be redistributed, shared, or embedded in any tool or dashboard. `data/raw/dhs/` is gitignored, as are `*.DTA`. Any resulting publication PDF goes to references@dhsprogram.com.

**GPS displacement**, restated from the approval letter because it is a modelling constraint, not a footnote: urban clusters are displaced up to 2 km, rural up to 5 km, and 1% of rural clusters up to 10 km, restricted to stay within the second administrative level where possible. DHS state explicitly that *measuring direct distance from a GPS location to another site is not appropriate*. This is why `02_covariates.R` reads every covariate over a displacement-matched buffer rather than at the point, and why no commune-level result should be read as if cluster positions were exact.

### Boundaries & population — DONE

- **Admin boundaries** — `03_Boundaries_2025.rar`, extracted to `data/raw/boundaries/`. adm1 region (24), adm2 district (120), adm3 commune (1,701), adm4 fokontany (17,470). Source: BNGRC/OCHA via PAM, 2025 edition.
- **Population raster** — WorldPop 2020 100 m constrained, total (UN-adjusted) and under-5 (f_0 + f_1 + m_0 + m_1). Downloaded automatically by `01_boundaries_pop.R`.

### EO / gridded covariates

Acquisition lives in `R/02_covariates_helpers.R`, one getter per source, all cached under `data/raw/rasters/`.

| Source | Covariates | Status |
|---|---|---|
| SRTM (via `geodata`) | elevation, ruggedness (TRI) | automatic |
| WorldClim 2.1 | temp mean/seasonality/min-cold, precip annual/seasonality | automatic |
| ESA WorldCover 10 m | class fractions: crop, built, water, tree, grass … | automatic |
| WorldPop | population count/density | automatic |
| GLW4 (Harvard Dataverse) | cattle density, dasymetric, 5 arc-min | automatic |
| Weiss et al. 2018 (MAP WCS) | travel time to cities, 2015 | automatic |
| Weiss et al. 2020 (MAP WCS) | motorized travel time to healthcare, 2020 | automatic |
| Li et al. harmonized DMSP–VIIRS | night-time lights, 2021 | automatic |

All covariate rasters now download automatically; there is no manual step. Each getter caches under `data/raw/rasters/` and is shared between `02` (commune/cluster) and `05` (prediction grid). If a source ever goes away, the getter fails loudly and the covariate is dropped and reported in `outputs/tables/02b_covariate_coverage.csv` rather than silently becoming NA.

Three source substitutions worth recording, because they change how results must be read:

- **Night lights.** NOAA/EOG's annual VNL V2 moved behind an OAuth account, so the pipeline uses the harmonized DMSP–VIIRS series of Li et al. (*Sci Data* 2020, extended through 2024), openly hosted on figshare and VIIRS-derived from 2014 on. Values are harmonized **DMSP-like digital numbers on a 0–63 scale, not VIIRS radiances**. That scale saturates over bright cores (Antananarivo reaches DN 56), so this is a usable rural/peri-urban activity gradient but not a linear intensity measure in the capital.
- **Cattle.** GLW4 dasymetric (`Da`), not areal-weighted (`Aw`): `Aw` spreads district totals uniformly and would erase exactly the within-district contrast of interest. Resolution is 5 arc-minutes (~10 km), coarser than every other covariate — it resolves district-scale but not commune-scale variation, and the effect table must not over-read a fine-grained cattle signal.
- **Travel time.** Both the 2015 travel-time-to-cities surface (as specified) and the 2020 motorized travel-time-to-healthcare surface are retained. The second is the more direct H4 measure here, since the commune is the catchment unit for basic health centres, and it is a year closer to the DHS. Distance to markets and distance to care are different exposures and need not move together.

## Methodology

1. **Direct estimates** — `03`: stunting prevalence per area from KR, fully respecting the DHS design (strata `v022`, PSU `v001`, weights `v005/1e6`) via `survey`/`srvyr`; design variances transferred to the logit scale by the delta method.
2. **Primary model** — `04` SUMMER BYM2 area-level smoothed-direct; `05` SPDE/INLA geostatistical model on cluster GPS + covariates for the continuous surface.
3. **Comparison model** — `06` area-level Fay-Herriot in `emdi`, fitted twice (household-type covariates only vs + EO) to quantify what the EO covariates buy against the WB specification.
4. **Explanatory model (paradox core)** — `07`: nested sequence of spatial HAZ models; the shrinkage of the highland coefficient as each hypothesis block enters is the decomposition.
5. **Aggregation** — `08`: population-weighted aggregation of the surface to commune level, **over posterior draws** so credible intervals respect the spatial correlation between grid cells.
6. **Validation** — `09`: spatial *block* cross-validation (not random k-fold), WAIC/DIC/CPO/PIT, benchmarking to DHS regional estimates, coherence check against IHME.

**Deliverables:** commune stunting surface with credible intervals; covariate effect table (paradox decomposition); comparison maps vs WB & IHME.

## Repository structure

```
.
├── CLAUDE.md
├── README.md
├── run_all.R                # pipeline runner
├── renv.lock                # written by R/00b_renv_init.R
├── data/
│   ├── raw/                 # DHS (gitignored), boundaries, rasters
│   ├── interim/
│   └── processed/
├── R/
│   ├── 00_setup.R           # paths, packages, CRS, helpers
│   ├── 00b_renv_init.R      # opt-in renv setup
│   ├── 01_boundaries_pop.R  # boundaries + WorldPop + adjacency   [DONE]
│   ├── 02_covariates_helpers.R  # raster getters, shared by 02 and 05
│   ├── 02_covariates.R      # EO stack + commune/cluster extraction   [DONE]
│   ├── 02b_covariate_diagnostics.R  # checks, highland profile, collinearity
│   ├── 03_dhs_direct.R      # read KR/PR/IR, design, direct estimates
│   ├── 04_model_summer.R    # BYM2 area-level
│   ├── 05_model_spde.R      # SPDE/INLA continuous surface
│   ├── 06_model_fh.R        # Fay-Herriot comparison (emdi)
│   ├── 07_explain_paradox.R # HAZ covariate model + decomposition
│   ├── 08_aggregate.R       # population-weighted commune aggregation
│   ├── 09_validate.R        # spatial block CV, WAIC, benchmarking
│   └── 10_maps.R            # output maps & figures
├── outputs/
│   ├── figures/
│   └── tables/
└── .gitignore               # excludes data/raw, DHS files, large rasters
```

## Environment & reproducibility

- R 4.5.1 with: `SUMMER`, `INLA` (own repo), `survey`, `srvyr`, `emdi`, `sae`, `sf`, `terra`, `exactextractr`, `spdep`, `geodata`, `dplyr`, `ggplot2`, `haven`.
- `renv` is **opt-in** via `R/00b_renv_init.R`, not wired into the analysis scripts — see the rationale in that file. Run `renv::snapshot()` after any package install.
- If any Python is used for GEE extraction (`rgee` or the GEE Python API), keep a frozen `requirements.txt`; use `python -m pip`; venv on Windows.
- Git: commit incrementally; push source files only (`.R`, `.md`, `.py`) — no notebooks, no zip archives. Portfolio: github.com/valofils.
- `.gitignore` excludes all DHS microdata and large rasters (DHS terms forbid redistribution).

## Known data issues (handled in `01`, do not "re-fix")

The three boundary layers were edited at different dates and disagree. `01_boundaries_pop.R` repairs them explicitly and writes an audit trail to `outputs/tables/`:

1. Antananarivo's six arrondissements carry an 11-character `ADM2_PCODE` in adm3 (`MG11101001A`) but 8 characters in adm2 (`MG11101A`–`F`). Crosswalked by name.
2. The 2021 split of Vatovavy-Fitovinany is half-applied: adm2/adm3 code Fitovinany as `MG23A`, adm1 codes it `MG35`. Two communes of Manakara Atsimo still carry the pre-split `MG23`. Region codes are harmonised to the adm1 layer by name; adm2 is treated as authoritative for the commune → district → region chain.
3. Region **Ambatosoa** (`MG34`) exists in adm1 but has no districts or communes in adm2/adm3, so results cannot be reported for it separately. No population is lost — its territory sits under Analanjirofo.
4. Commune **names are not unique** (166 duplicates). Always join on `ADM3_PCODE`.
5. The commune adjacency graph has 4 components (offshore islands). `01` adds the minimum number of artificial edges to connect it, because an ICAR/BYM2 prior on a disconnected graph is improper. Edges added are logged to `outputs/tables/01_adjacency_edges_added.csv`.

## Analytical constraints established by `02b` (read before interpreting `07`)

Measured on the 1,701-commune covariate stack (population-weighted where relevant).

### The highlands are not the poor, remote periphery — the paradox is sharper than it looks

Highland (>800 m) vs lowland, weighted by under-5 population:

| | Highland | Lowland | direction |
|---|---|---|---|
| Elevation (m) | 1,254 | 232 | |
| Min temp, coldest month (°C) | 9.0 | 14.4 | **colder** |
| Cropland fraction | 0.292 | 0.075 | **3.9× more cultivated** |
| Tree cover fraction | 0.098 | 0.305 | less forest |
| Cattle per 10 km cell | 2,184 | 1,531 | **more cattle** |
| Travel time to cities (min) | 200 | 390 | **half the distance** |
| Travel time to healthcare (min) | 55 | 81 | **better served** |
| Night lights (DN) | 6.42 | 0.51 | **12× brighter** |
| Population density (per km²) | 2,368 | 157 | denser |

This is the single most important thing the covariates say, and it constrains
the whole explanation. On every conventional deprivation axis — market access,
health-facility access, economic activity, livestock wealth, agricultural land —
**the highlands are better off than the lowlands**, yet they carry the higher
stunting burden. So:

- **H4 (access) cannot be the mediator in the expected direction.** Adjusting
  for travel time should make the highland penalty *larger*, not smaller. If the
  decomposition shows H4 "explaining" part of the gap, check the sign before
  reporting it.
- **H2 needs restating.** Cattle are *more* abundant in the highlands, so a
  simple "no animals, no animal-source food" story fails. The plausible
  mechanism is that Malagasy zebu function largely as stored wealth and ritual
  capital rather than as a dietary source for young children — presence is not
  consumption. Cattle density measures the wrong construct for H2 and should be
  interpreted as a wealth proxy unless a consumption variable (DHS dietary
  diversity, IR recode) is brought in to carry H2 properly.
- What survives as candidate explanations is therefore **H1 (altitude/cold)** and
  **H3 (infection load in a cold, dense, paddy-irrigated, heavily cultivated
  landscape)** — plus diet *quality* rather than diet *quantity*.

### Two collinearity limits

- Elevation correlates **−0.93** with mean temperature and **−0.90** with the
  coldest-month minimum. Altitude and cold are effectively the same variable in
  Madagascar. `07` **cannot** separate "cold stress" from "altitude per se";
  block A must be reported as a single altitude–temperature construct, not as
  evidence for a thermal mechanism specifically.
- `frac_built` and `nightlights` correlate **0.95**. They are one variable, not
  two, and they sit in different hypothesis blocks (H3 and H4 respectively).
  Entering both inflates the apparent contribution of whichever block is entered
  first. Keep one per model, or report the pair jointly.

Cropland (0.52), cattle (0.30), population density (0.14) and travel time to
healthcare (−0.04) are all weakly enough related to elevation that H2, H3 and H4
remain separable from H1. 40% of Madagascar's under-5 population lives above
800 m, so the highland contrast is between two large groups.

## Working conventions for Claude Code

- Proceed autonomously — do not ask for per-task "allow" confirmation; batch the work.
- **Never edit an `.R` file while `Rscript` is executing it.** R streams the file from disk, so an edit mid-run corrupts the parse and produces a baffling syntax error at a line you did not touch.
- Write `03_*`–`10_*` as runnable skeletons with clear file-path placeholders so they execute the moment the DHS files land.

## Key references (NLM/Vancouver)

Matekenya D, Mulangu FM, Newhouse D. Malnourished but not destitute: the spatial interplay between nutrition and poverty in Madagascar. Washington (DC): World Bank; 2023. Policy Research Working Paper No. 10627.

Local Burden of Disease Child Growth Failure Collaborators. Mapping child growth failure across low- and middle-income countries. Nature. 2020;577(7789):231-4.

Molina I, Rao JNK. Small area estimation of poverty indicators. Can J Stat. 2010;38(3):369-85.

Kreutzmann AK, Pannier S, Rojas-Perilla N, Schmid T, Templ M, Tzavidis N. The R package emdi for estimating and mapping regionally disaggregated indicators. J Stat Softw. 2019;91(7):1-33.

Wakefield J, Okonek T, Pedersen J. Small area estimation for disease prevalence mapping. Int Stat Rev. 2020;88(2):398-418.

Lindgren F, Rue H, Lindström J. An explicit link between Gaussian fields and Gaussian Markov random fields: the SPDE approach. J R Stat Soc Series B Stat Methodol. 2011;73(4):423-98.

Riebler A, Sørbye SH, Simpson D, Rue H. An intuitive Bayesian spatial model for disease mapping that accounts for scaling. Stat Methods Med Res. 2016;25(4):1145-65.

Mercer LD, Wakefield J, Pantazis A, Lutambi AM, Masanja H, Clark S. Space-time smoothing of complex survey data: small area estimation for child mortality. Ann Appl Stat. 2015;9(4):1889-905.

Pérez-Heydrich C, Warren JL, Burgert CR, Emch ME. Guidelines on the use of DHS GPS data. Calverton (MD): ICF International; 2013. Spatial Analysis Reports No. 8.

Weiss DJ, Nelson A, Gibson HS, et al. A global map of travel time to cities to assess inequalities in accessibility in 2015. Nature. 2018;553(7688):333-6.
