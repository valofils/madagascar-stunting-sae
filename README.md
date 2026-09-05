# Geospatial Small-Area Estimation of Child Stunting in Madagascar

Commune-level estimates of child stunting (height-for-age Z < −2) for Madagascar,
from the 2021 DHS (EDSMD-V), using a spatially explicit small-area estimation
pipeline enriched with earth-observation covariates — and an attempt to explain
the **fertile-highland paradox**: stunting is highest in the country's most
agriculturally productive region.

---

## The question

Madagascar's central highlands are its agricultural heart, yet they carry its
heaviest burden of child stunting. Poverty and undernutrition, which normally
track each other closely, come apart here. The World Bank's 2023 commune-level
SAE (WP 10627) documented the pattern but reported a stunting model with a
marginal R² of about 0.075, and stated plainly that it could not explain the
paradox. IHME's Local Burden of Disease surfaces map child growth failure at 5 km
but stop before the 2021 DHS, calibrate covariates continentally rather than to
Malagasy agro-ecology, and aggregate to district rather than commune.

This project takes a different route: model the paradox rather than only map it.

## Approach

Four hypotheses drive the covariate set, and the analysis is built so that each
can be tested against the others:

| | Hypothesis | Operationalised as |
|---|---|---|
| H1 | Altitude / cold stress | elevation, ruggedness, min temperature of coldest month, temperature seasonality |
| H2 | Rice monoculture, low dietary diversity | cropland fraction, precipitation and its seasonality, cattle density |
| H3 | Infection load / environmental enteropathy | population density, built-up fraction, permanent-water fraction, recent diarrhoea |
| H4 | Care practices, maternal condition, access | maternal education, wealth, BMI, age at first birth, birth order, travel time to cities |

The decomposition in `R/07_explain_paradox.R` fits a nested sequence of spatial
models for HAZ and tracks a single number — the coefficient on *highland*
(elevation > 800 m) — as each hypothesis block enters. The share of the raw
highland penalty that a block absorbs is that block's explanatory contribution.
Every model carries an SPDE spatial field, so a highland effect that survives is
an altitude signal rather than a relabelled "central Madagascar" indicator.

This is a mediation-style decomposition, not a causal identification. A covariate
that absorbs the penalty is a candidate mechanism, not a proven cause.

## Pipeline

```bash
Rscript run_all.R          # everything the available data supports
Rscript run_all.R 1 2      # selected stages only
```

| Stage | Script | What it does | Needs DHS |
|---|---|---|---|
| 1 | `01_boundaries_pop.R` | harmonise adm1–adm3, WorldPop population, commune adjacency graph | no |
| 2 | `02_covariates.R` | EO raster stack → commune and DHS-cluster covariates | partly |
| 2b | `02b_covariate_diagnostics.R` | covariate integrity checks, highland profile, collinearity | no |
| 3 | `03_dhs_direct.R` | design-based direct estimates (`survey`/`srvyr`) | yes |
| 4 | `04_model_summer.R` | SUMMER BYM2 smoothed-direct, district level | yes |
| 5 | `05_model_spde.R` | binomial SPDE/INLA geostatistical surface, 1 km | yes |
| 6 | `06_model_fh.R` | Fay–Herriot comparison in `emdi` | yes |
| 7 | `07_explain_paradox.R` | the paradox decomposition | yes |
| 8 | `08_aggregate.R` | population-weighted aggregation to commune, over posterior draws | yes |
| 9 | `09_validate.R` | spatial block CV, WAIC/CPO/PIT, benchmarking | yes |
| 10 | `10_maps.R` | final figures | yes |

Stage 1 is complete. Stage 2 runs now for communes and completes for clusters
once the DHS GPS file lands. Stages 3–10 are written and ready; they fail with an
explicit message naming the missing file until the microdata arrive.

## Three decisions worth knowing about

**Aggregation is over posterior draws, not the posterior mean.** Averaging a mean
surface into a commune gives the right point estimate and the wrong interval:
neighbouring grid cells share a spatial field, so their errors are correlated and
do not average away. `08_aggregate.R` therefore aggregates 500 posterior draws.

**Cross-validation is spatially blocked, not random.** Random k-fold is optimistic
for a spatially correlated model, because a held-out cluster usually has a
neighbour still in training. `09_validate.R` holds out contiguous 60 km tiles, so
the model has to extrapolate the distance it will really face when predicting an
unsampled commune.

**Covariates are read over a displacement-matched buffer, not at the point.** DHS
GPS coordinates are displaced up to 2 km (urban) or 5 km (rural). Sampling a
raster at the published point reads the wrong place; `02_covariates.R` averages
over a buffer matching the displacement radius.

## Data

- **DHS 2021 Madagascar (EDSMD-V)** — request pending. Place the KR/PR/IR recodes
  and the GE shapefile under `data/raw/dhs/`. DHS terms forbid redistribution;
  the directory is gitignored.
- **Boundaries** — BNGRC/OCHA via PAM, 2025 edition: 24 regions, 120 districts,
  1,701 communes, 17,470 fokontany. Supplied as `03_Boundaries_2025.rar`.
- **WorldPop 2020** 100 m constrained, total and under-5 — downloaded automatically.
- **EO covariates** — SRTM, WorldClim, ESA WorldCover downloaded automatically;
  FAO GLW4 cattle, Weiss travel time and VIIRS night lights need a manual
  download (the scripts print the URL and target path, and run without them).

Three of the boundary layers disagree with each other in specific, reproducible
ways — an incomplete 2021 region split, inconsistent Antananarivo arrondissement
codes, and a region with no communes. `01_boundaries_pop.R` repairs all of them
explicitly and writes an audit trail to `outputs/tables/`. See the *Known data
issues* section of [CLAUDE.md](CLAUDE.md) before "fixing" anything there.

## Environment

R 4.5.1 on Windows. Key packages: `INLA`, `SUMMER`, `emdi`, `survey`, `srvyr`,
`sf`, `terra`, `exactextractr`, `spdep`, `geodata`.

`renv` setup is opt-in and deliberately not wired into the analysis scripts:

```bash
Rscript R/00b_renv_init.R
```

The rationale is documented at the top of that file — `renv::init()` switches R to
a project-private library, which is a decision rather than something an analysis
script should do as a side effect.

## Outputs

- `data/processed/commune_stunting.gpkg` — commune estimates with credible
  intervals, benchmarked and unbenchmarked
- `outputs/tables/07_paradox_decomposition.csv` — the decomposition
- `outputs/tables/08_commune_burden_ranking.csv` — communes ranked by absolute
  number of stunted children, the form a targeting decision actually needs
- `outputs/figures/` — commune map, uncertainty and burden panel, paradox panel,
  estimator comparison

## References

Matekenya D, Mulangu FM, Newhouse D. *Malnourished but not destitute: the spatial
interplay between nutrition and poverty in Madagascar.* World Bank Policy Research
Working Paper 10627; 2023.

Local Burden of Disease Child Growth Failure Collaborators. Mapping child growth
failure across low- and middle-income countries. *Nature.* 2020;577(7789):231-4.

Lindgren F, Rue H, Lindström J. An explicit link between Gaussian fields and
Gaussian Markov random fields: the SPDE approach. *J R Stat Soc B.* 2011;73(4):423-98.

Riebler A, Sørbye SH, Simpson D, Rue H. An intuitive Bayesian spatial model for
disease mapping that accounts for scaling. *Stat Methods Med Res.* 2016;25(4):1145-65.

Wakefield J, Okonek T, Pedersen J. Small area estimation for disease prevalence
mapping. *Int Stat Rev.* 2020;88(2):398-418.

Full reference list in [CLAUDE.md](CLAUDE.md).
