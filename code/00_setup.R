# =============================================================================
# 00_setup.R: Environment Setup and Package Installation
# =============================================================================
#
# Purpose : Install and load all required packages, verify R version and
#           CmdStan installation, define project paths, and set random seed
#           for reproducibility.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
# =============================================================================


###############################################################################
## Section 1 : R Version Check
###############################################################################

cat("=== svyder-replication: Environment Setup ===\n\n")

## Require R >= 4.5
r_version <- getRversion()
if (r_version < "4.5.0") {
  stop(sprintf(
    "R version >= 4.5.0 is required (current: %s). ",
    as.character(r_version)
  ))
}
cat(sprintf("  R version: %s  [OK]\n", as.character(r_version)))


###############################################################################
## Section 2 : Package Installation and Loading
###############################################################################

## --- CRAN packages ---
## These are the core packages required for the replication.

cran_packages <- c(
  "ggplot2",     # plotting (>= 3.4.0 for linewidth support)
  "dplyr",       # data manipulation
  "tidyr",       # data reshaping
  "purrr",       # functional programming
  "patchwork",   # multi-panel figure composition
  "scales",      # axis formatting for ggplot2
  "Matrix",      # nearPD for matrix conditioning
  "mvtnorm",     # multivariate normal (Cholesky correction)
  "xtable",      # LaTeX table generation
  "here",        # project-relative paths
  "remotes"      # GitHub package installation
)

cat("\nChecking CRAN packages ...\n")

installed <- rownames(installed.packages())
to_install <- setdiff(cran_packages, installed)

if (length(to_install) > 0L) {
  cat(sprintf("  Installing %d packages: %s\n",
              length(to_install), paste(to_install, collapse = ", ")))
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  cat("  All CRAN packages already installed  [OK]\n")
}

## Load all CRAN packages
invisible(lapply(cran_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))
cat("  All CRAN packages loaded  [OK]\n")


## --- cmdstanr (Stan interface) ---

cat("\nChecking cmdstanr ...\n")

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  cat("  Installing cmdstanr from mc-stan.org ...\n")
  install.packages("cmdstanr",
                   repos = c("https://mc-stan.org/r-packages/",
                             getOption("repos")))
}

library(cmdstanr)
cat(sprintf("  cmdstanr version: %s  [OK]\n",
            as.character(packageVersion("cmdstanr"))))


## --- CmdStan (Stan compiler) ---

cat("\nChecking CmdStan installation ...\n")

cmdstan_path <- tryCatch(cmdstan_path(), error = function(e) NULL)

if (is.null(cmdstan_path)) {
  cat("  CmdStan not found. Installing (this may take ~10 minutes) ...\n")
  install_cmdstan()
  cmdstan_path <- cmdstan_path()
}

cat(sprintf("  CmdStan path: %s  [OK]\n", cmdstan_path))
cat(sprintf("  CmdStan version: %s\n", cmdstan_version()))


## --- svyder package (from GitHub) ---

cat("\nChecking svyder package ...\n")

if (!requireNamespace("svyder", quietly = TRUE)) {
  cat("  Installing svyder from GitHub ...\n")
  remotes::install_github("joonho112/svyder")
}

library(svyder)
cat(sprintf("  svyder version: %s  [OK]\n",
            as.character(packageVersion("svyder"))))


###############################################################################
## Section 3 : Project Paths
###############################################################################

cat("\nSetting project paths ...\n")

## Use here::here() for project-relative paths. This requires that
## svyder-replication.Rproj exists in the project root.
PROJECT_ROOT <- here::here()

PATHS <- list(
  root       = PROJECT_ROOT,
  code       = file.path(PROJECT_ROOT, "code"),
  helpers    = file.path(PROJECT_ROOT, "code", "helpers"),
  figures    = file.path(PROJECT_ROOT, "code", "03_figures"),
  stan       = file.path(PROJECT_ROOT, "stan"),
  data       = file.path(PROJECT_ROOT, "data"),
  precomputed = file.path(PROJECT_ROOT, "data", "precomputed"),
  output     = file.path(PROJECT_ROOT, "output"),
  fig_out    = file.path(PROJECT_ROOT, "output", "figures"),
  tab_out    = file.path(PROJECT_ROOT, "output", "tables")
)

## Create output directories if they do not exist
for (dir_path in c(PATHS$fig_out, PATHS$tab_out)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

cat(sprintf("  Project root: %s\n", PATHS$root))


###############################################################################
## Section 4 : Random Seed
###############################################################################

## Master random seed for reproducibility.
## All simulation and MCMC scripts use this seed as their base.
MASTER_SEED <- 20260303L
set.seed(MASTER_SEED)

cat(sprintf("  Random seed: %d  [OK]\n", MASTER_SEED))


###############################################################################
## Section 5 : Global Constants
###############################################################################

## DER threshold for three-tier classification (Paper Section 3.2)
DER_THRESHOLD_LO <- 0.80   # Below: Tier I (Survey-Robust)
DER_THRESHOLD_HI <- 1.20   # Above: Tier III (Survey-Inflated)

## Default prior SD for fixed effects (Paper Section 2)
BETA_PRIOR_SD <- 5.0       # beta ~ N(0, 25)


###############################################################################
## Section 6 : Session Info
###############################################################################

cat("\n=== Setup complete ===\n\n")
cat("Session info:\n")
print(sessionInfo())
