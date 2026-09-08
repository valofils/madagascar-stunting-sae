# ---------------------------------------------------------------------------
# run_all.R  --  Run the pipeline end to end.
#
#   Rscript run_all.R          # everything that can run with the data present
#   Rscript run_all.R 1 2      # only the numbered stages given
#
# Stages that need the DHS microdata fail with a clear message and are
# reported as skipped rather than aborting the run, so the non-DHS half of the
# pipeline can be re-run at any time while the data request is pending.
# ---------------------------------------------------------------------------

source("R/00_setup.R")

STAGES <- c(
  "1"  = "R/01_boundaries_pop.R",
  "2"  = "R/02_covariates.R",
  "2b" = "R/02b_covariate_diagnostics.R",
  "3b" = "R/03b_dhs_region_crosswalk.R",
  "3"  = "R/03_dhs_direct.R",
  "4"  = "R/04_model_summer.R",
  "5"  = "R/05_model_spde.R",
  "6"  = "R/06_model_fh.R",
  "7"  = "R/07_explain_paradox.R",
  "8"  = "R/08_aggregate.R",
  "9"  = "R/09_validate.R",
  "10" = "R/10_maps.R"
)
NEEDS_DHS <- c("3b", "3", "4", "5", "6", "7", "8", "9", "10")   # "1", "2", "2b" run without it

args <- commandArgs(trailingOnly = TRUE)
run <- if (length(args) > 0) args else names(STAGES)

have_dhs <- length(list.files(DIR$dhs, pattern = "[.]DTA$", recursive = TRUE,
                              ignore.case = TRUE)) > 0
if (!have_dhs)
  msg("NOTE: no DHS microdata under ", DIR$dhs,
      " - stages 3-10 will be reported as blocked.")

results <- list()
for (s in run) {
  f <- STAGES[[s]]
  if (is.null(f)) { msg("unknown stage: ", s); next }
  msg("=================== stage ", s, ": ", basename(f), " ===================")
  t0 <- Sys.time()
  ok <- tryCatch({ source(f, echo = FALSE); TRUE },
                 error = function(e) { msg("FAILED: ", conditionMessage(e)); FALSE })
  results[[s]] <- data.frame(
    stage = s, script = basename(f),
    status = if (ok) "ok" else if (!have_dhs && s %in% NEEDS_DHS) "blocked (no DHS)" else "failed",
    minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
}

summary_df <- do.call(rbind, results)
utils::write.csv(summary_df, file.path(DIR$tables, "00_pipeline_status.csv"),
                 row.names = FALSE)
msg("pipeline summary:")
print(summary_df, row.names = FALSE)
