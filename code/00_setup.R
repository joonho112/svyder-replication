# =============================================================================
# 00_setup.R: Environment Setup and Package Installation
# =============================================================================
#
# Purpose : Install and load packages for the selected replication track,
#           define project paths, and set a random seed for reproducibility.
#           Track B is the default public path and does not install Stan.
# Paper   : Lee, J., Williams, M. R., & Savitsky, T. D. (2026). Design Effect
#           Ratios for Bayesian Survey Models: A Diagnostic Framework for
#           Identifying Survey-Sensitive Parameters. arXiv preprint.
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
# =============================================================================


###############################################################################
## Section 1 : R Version Check
###############################################################################

cat("=== svyder-replication: Environment Setup ===\n\n")

SETUP_TRACK <- toupper(Sys.getenv("SVYDER_TRACK", "B"))
if (!SETUP_TRACK %in% c("A", "B")) {
  stop("SVYDER_TRACK must be either 'A' or 'B' (current: ",
       SETUP_TRACK, ").", call. = FALSE)
}
IS_TRACK_A <- identical(SETUP_TRACK, "A")
cat(sprintf("  Selected track: Track %s\n", SETUP_TRACK))
if (!IS_TRACK_A) {
  cat("  Track B setup skips cmdstanr, CmdStan, and svyder.\n")
  cat("  For Track A, run: Sys.setenv(SVYDER_TRACK = 'A'); ",
      "source('code/00_setup.R')\n", sep = "")
}

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
## Track B is the public reporting/book path. Track A adds the packages needed
## for from-scratch simulation and restricted-data application refits.

track_b_cran_packages <- c(
  "ggplot2",     # plotting (>= 3.4.0 for linewidth support)
  "patchwork",   # multi-panel figure composition
  "scales",      # axis formatting for ggplot2
  "here",        # project-relative paths
  "knitr",       # Quarto/R Markdown table rendering
  "rmarkdown"    # Quarto/R Markdown execution support
)

track_a_extra_cran_packages <- c(
  "dplyr",       # data manipulation
  "tidyr",       # data reshaping
  "purrr",       # functional programming
  "Matrix",      # nearPD for matrix conditioning
  "mvtnorm",     # multivariate normal (Cholesky correction)
  "xtable",      # LaTeX table generation
  "posterior",   # posterior draw handling in the NSECE refit
  "remotes"      # GitHub package installation
)

cran_packages <- if (IS_TRACK_A) {
  unique(c(track_b_cran_packages, track_a_extra_cran_packages))
} else {
  track_b_cran_packages
}

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
cat("  All selected CRAN packages loaded  [OK]\n")


if (IS_TRACK_A) {
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
} else {
  cat("\nSkipping Track A dependencies (cmdstanr, CmdStan, svyder)  [OK]\n")
}


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
  results    = file.path(PROJECT_ROOT, "code", "03_results"),
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

## A base random seed for interactive work in this session. The simulation and
## application scripts set their own seeds internally (see sim_00_config.R and
## the application scripts), so a full rerun does not depend on this value.
MASTER_SEED <- 20260303L
set.seed(MASTER_SEED)

cat(sprintf("  Random seed: %d  [OK]\n", MASTER_SEED))


###############################################################################
## Section 5 : Global Constants
###############################################################################

## Default flagging threshold c0 for the compute-classify-correct workflow
## (Algorithm 1): a parameter is flagged when its DER exceeds this value.
DER_THRESHOLD <- 1.20

## Prior SD for the fixed effects in stan/hlr_weighted.stan.
BETA_PRIOR_SD <- 5.0


###############################################################################
## Section 6 : Session Info
###############################################################################

cat("\n=== Setup complete ===\n\n")
cat("Session info:\n")
print(sessionInfo())
