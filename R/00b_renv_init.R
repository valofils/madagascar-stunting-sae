# ---------------------------------------------------------------------------
# 00b_renv_init.R
#
# OPT-IN. Run this once, deliberately, from a fresh R session in the project
# root. It is deliberately NOT sourced by any other script.
#
# Why it is separate: renv::init() switches R to a project-private library.
# Everything this pipeline needs (INLA above all, which is large and installed
# from its own repository) then has to be present in that library or the
# scripts stop working. renv::hydrate() copies what is already installed in
# the user library across, so nothing is re-downloaded, but this is still a
# change to how R resolves packages on this machine and should be a decision
# rather than a side effect of running an analysis script.
#
#   Rscript R/00b_renv_init.R
#
# Afterwards, re-run it (or just renv::snapshot()) after installing anything
# new, so renv.lock stays an accurate record of the environment.
# ---------------------------------------------------------------------------

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

PKGS <- c(
  # spatial
  "sf", "terra", "exactextractr", "spdep", "geodata",
  # survey and small-area estimation
  "survey", "srvyr", "emdi", "sae", "SUMMER", "INLA",
  # data handling
  "dplyr", "tidyr", "stringr", "readr", "haven", "rlang",
  # output
  "ggplot2", "scales", "patchwork", "viridis"
)

message("renv::init(bare = TRUE) - creating the project library")
renv::init(bare = TRUE, restart = FALSE)

message("renv::hydrate() - copying packages from the user library")
# hydrate links or copies from the existing user library rather than
# downloading, which is what keeps this fast and keeps INLA working.
renv::hydrate(packages = PKGS)

missing <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("not yet available, installing: ", paste(missing, collapse = ", "))
  # INLA is not on CRAN and needs its own repository.
  if ("INLA" %in% missing) {
    renv::install("INLA", repos = c(getOption("repos"),
                                    INLA = "https://inla.r-inla-download.org/R/stable"))
    missing <- setdiff(missing, "INLA")
  }
  if (length(missing) > 0) renv::install(missing)
}

renv::snapshot(packages = PKGS, prompt = FALSE)
message("renv.lock written. Re-run renv::snapshot() after any further install.")
