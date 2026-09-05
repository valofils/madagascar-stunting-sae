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

### DHS 2021 Madagascar (EDSMD-V) — status: request PENDING approval

Registered at dhsprogram.com; Survey + GPS both requested for Madagascar. Files to download once approved, into `data/raw/dhs/`:

- **KR** (Children's Recode) — anthropometry / HAZ. Primary outcome source.
- **PR** (Household Member Recode) — WASH, household context.
- **IR** (Individual/Women's Recode) — maternal education, dietary diversity, ANC.
- **GE** (Geographic / GPS cluster coordinates) — required for covariate linkage & spatial models. Coordinates are DHS-displaced; never used to identify households.

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
| FAO GLW4 | cattle density | **manual** — see below |
| Weiss et al. 2018 | travel time to cities | **manual** |
| VIIRS annual composite | night-time lights | **manual** |

Manual rasters stop the getter with an explicit message naming the URL and the exact target path. The pipeline runs without them; the affected covariates are simply dropped and reported in `outputs/tables/02_covariate_coverage.csv`.

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

## Analytical constraint established by `02b` (read before interpreting `07`)

Measured on the 1,701-commune covariate stack:

- Elevation correlates **−0.93** with mean temperature and **−0.90** with the
  coldest-month minimum. Altitude and cold are effectively the same variable in
  Madagascar. The decomposition in `07` therefore **cannot** separate "H1 cold
  stress" from "altitude per se" — block A should be reported as a single
  altitude–temperature construct, not as evidence for a thermal mechanism
  specifically.
- Elevation correlates only **0.52** with cropland fraction, **0.14** with
  population density and **0.13** with built-up fraction. H2 and H3 *are*
  separable from H1, so the interesting attribution — how much of the highland
  penalty is diet, how much is infection load — is identifiable.
- 40% of Madagascar's under-5 population lives above 800 m, so the highland
  contrast is between two large groups, not a small subgroup against the rest.

Population-weighted highland (>800 m) vs lowland: coldest-month minimum 9.0 vs
14.4 °C, cropland fraction 0.29 vs 0.075, tree cover 0.098 vs 0.305.

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
