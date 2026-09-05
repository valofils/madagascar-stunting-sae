# ---------------------------------------------------------------------------
# 06_model_fh.R
#
# Area-level Fay-Herriot model in emdi, fitted at district level. This is the
# comparison arm: it is the classical, non-spatial SAE estimator, and it is
# the family the World Bank working paper (Matekenya et al. 2023) used.
#
# Its role here is evidential, not predictive. Three contrasts are produced:
#   (1) FH with household-type covariates only, echoing the WB specification
#       - expected to be weak, as they reported (marginal R2 ~ 0.075);
#   (2) FH with the earth-observation covariates added;
#   (3) FH vs the spatial BYM2 fit from 04.
# The gap between (1) and (2) is the measurable value of the EO covariates;
# the gap between (2) and (3) is the value of modelling space explicitly.
#
# Inputs : data/processed/direct_estimates_adm2.csv, commune_covariates.csv,
#          summer_adm2.csv
# Outputs: data/processed/fh_adm2.csv, outputs/tables/06_*.csv
# ---------------------------------------------------------------------------

source("R/00_setup.R")
need("dplyr", "sf", "emdi", "ggplot2")

set.seed(20210)

# ===========================================================================
# 1. Model frame
# ===========================================================================
direct <- utils::read.csv(OUT$direct_adm2)
adm2 <- sf::st_read(OUT$adm2, quiet = TRUE)
cov_com <- utils::read.csv(OUT$cov_commune)

EO_COVARS <- c("elevation", "ruggedness", "temp_min_cold", "temp_seasonality",
               "precip_annual", "precip_seasonality", "frac_crop", "frac_built",
               "frac_water_perm", "cattle_density", "travel_time", "nightlights")
# Stand-ins for the census/household variables the WB model was limited to.
# They are all derivable without any EO data, which is the point of the contrast.
BASE_COVARS <- c("pop_dens", "urban_share")

have_eo <- intersect(EO_COVARS, names(cov_com))

X2 <- cov_com |>
  dplyr::mutate(w = pmax(pop_u5, 1e-6),
                urban_share = if ("frac_built" %in% names(cov_com)) frac_built else NA_real_) |>
  dplyr::group_by(ADM2_PCODE) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(c(have_eo, "pop_dens", "urban_share")),
                                 ~ stats::weighted.mean(.x, w, na.rm = TRUE)),
                   pop_u5 = sum(pop_u5, na.rm = TRUE), .groups = "drop")

frame <- adm2 |>
  sf::st_drop_geometry() |>
  dplyr::select(ADM2_PCODE, ADM2_EN, ADM1_PCODE, ADM1_EN) |>
  dplyr::left_join(direct, by = "ADM2_PCODE") |>
  dplyr::left_join(X2, by = "ADM2_PCODE")

# emdi's fh() needs one row per area, the direct estimate, and its sampling
# variance. Areas without a direct estimate are carried as out-of-sample and
# predicted from the synthetic part of the model.
frame <- frame |>
  dplyr::mutate(dplyr::across(dplyr::all_of(c(have_eo, BASE_COVARS)),
                              ~ as.numeric(scale(.x))))

# The Fay-Herriot variance argument must be strictly positive and finite.
# Districts with a single cluster have a degenerate (zero or NA) design
# variance and cannot enter the fit; they are predicted instead.
frame$usable <- is.finite(frame$logit_direct) & is.finite(frame$logit_var) &
  frame$logit_var > 0 & !is.na(frame$n_clusters) & frame$n_clusters >= 2

msg("districts total ", nrow(frame), " | usable in the FH fit ", sum(frame$usable),
    " | predicted out of sample ", sum(!frame$usable))

# ===========================================================================
# 2. Smoothed sampling variances
# ===========================================================================
# The FH model assumes the sampling variances are KNOWN. With 3-8 clusters per
# district they are themselves very noisy, and plugging in the raw values
# transfers that noise into the shrinkage weights. The standard remedy is a
# generalised variance function: regress log(design variance) on log(sample
# size) and use the fitted value.
vf_dat <- frame[frame$usable, ]
vf <- stats::lm(log(logit_var) ~ log(n_children), data = vf_dat)
frame$logit_var_smooth <- exp(stats::predict(vf, newdata = frame))
msg("variance smoothing: slope on log(n) = ", round(stats::coef(vf)[2], 3),
    " (theory says about -1), R2 = ", round(summary(vf)$r.squared, 3))
utils::capture.output(summary(vf),
                      file = file.path(DIR$tables, "06_variance_function.txt"))

# ===========================================================================
# 3. Fit the two specifications
# ===========================================================================
fit_fh <- function(covars, label) {
  covars <- covars[covars %in% names(frame)]
  covars <- covars[vapply(covars, function(v) !anyNA(frame[[v]]), logical(1))]
  if (length(covars) == 0) { msg("no usable covariates for ", label); return(NULL) }

  fml <- stats::as.formula(paste("logit_direct ~", paste(covars, collapse = " + ")))
  d <- frame[, c("ADM2_PCODE", "logit_direct", "logit_var_smooth", covars)]
  d <- d[stats::complete.cases(d[, covars]), ]

  msg("fitting FH [", label, "] with ", length(covars), " covariates on ",
      sum(!is.na(d$logit_direct)), " sampled districts")

  f <- tryCatch(
    emdi::fh(fixed = fml, vardir = "logit_var_smooth", combined_data = d,
             domains = "ADM2_PCODE", method = "ml", MSE = TRUE, B = c(0, 50),
             transformation = "no", eff_smpsize = NULL),
    error = function(e) { msg("FH [", label, "] failed: ", conditionMessage(e)); NULL })
  if (is.null(f)) return(NULL)

  r2 <- tryCatch(summary(f)$model$model_select, error = function(e) NULL)
  if (!is.null(r2)) {
    utils::write.csv(as.data.frame(r2),
                     file.path(DIR$tables, paste0("06_fh_", label, "_fit.csv")))
    msg("  [", label, "] marginal R2 ", round(r2$AdjR2 %||% NA_real_, 3))
  }
  utils::capture.output(summary(f),
    file = file.path(DIR$tables, paste0("06_fh_", label, "_summary.txt")))
  f
}

fh_base <- fit_fh(BASE_COVARS, "base")
fh_eo <- fit_fh(c(BASE_COVARS, have_eo), "eo")

# ===========================================================================
# 4. Collect predictions, back-transform to the prevalence scale
# ===========================================================================
collect <- function(f, label) {
  if (is.null(f)) return(NULL)
  e <- f$ind
  m <- f$MSE
  out <- data.frame(ADM2_PCODE = as.character(e$Domain),
                    logit_est = e$FH,
                    logit_mse = if (!is.null(m)) m$FH else NA_real_)
  # Delta method back to the probability scale:
  #   p = plogis(theta);  Var(p) = (dp/dtheta)^2 Var(theta) = (p(1-p))^2 Var(theta)
  out$est <- stats::plogis(out$logit_est)
  out$var <- (out$est * (1 - out$est))^2 * out$logit_mse
  out$se <- sqrt(out$var)
  names(out)[-1] <- paste0(label, "_", names(out)[-1])
  out
}

res <- frame[, c("ADM2_PCODE", "ADM2_EN", "ADM1_EN", "direct", "se",
                 "n_children", "n_clusters")]
for (p in list(list(fh_base, "fh_base"), list(fh_eo, "fh_eo"))) {
  cc <- collect(p[[1]], p[[2]])
  if (!is.null(cc)) res <- dplyr::left_join(res, cc, by = "ADM2_PCODE")
}

# Bring in the spatial fit from 04 for the three-way comparison.
f_summer <- file.path(DIR$processed, "summer_adm2.csv")
if (file.exists(f_summer)) {
  s <- utils::read.csv(f_summer)[, c("ADM2_PCODE", "summer_est", "summer_var")]
  res <- dplyr::left_join(res, s, by = "ADM2_PCODE")
}

utils::write.csv(res, file.path(DIR$processed, "fh_adm2.csv"), row.names = FALSE)

# ===========================================================================
# 5. Efficiency comparison
# ===========================================================================
# The headline number: how much narrower are the intervals than the direct
# estimator, model by model. Reported as a median across districts.
eff <- data.frame(
  estimator = c("direct", "FH (household-type covariates)",
                "FH (+ earth observation)", "SUMMER BYM2 (spatial)"),
  median_se = c(
    stats::median(res$se, na.rm = TRUE),
    if ("fh_base_se" %in% names(res)) stats::median(res$fh_base_se, na.rm = TRUE) else NA,
    if ("fh_eo_se" %in% names(res)) stats::median(res$fh_eo_se, na.rm = TRUE) else NA,
    if ("summer_var" %in% names(res)) stats::median(sqrt(res$summer_var), na.rm = TRUE) else NA))
eff$se_ratio_vs_direct <- round(eff$median_se / eff$median_se[1], 3)
eff$effective_sample_gain <- round(1 / eff$se_ratio_vs_direct^2, 2)

utils::write.csv(eff, file.path(DIR$tables, "06_efficiency_comparison.csv"),
                 row.names = FALSE)
msg("efficiency comparison:")
print(eff, row.names = FALSE)

msg("06_model_fh.R complete")
