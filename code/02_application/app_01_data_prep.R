# =============================================================================
# app_01_data_prep.R: NSECE 2019 Data Preparation
# =============================================================================
#
# Purpose : Load NSECE 2019 restricted-use data, construct the Stan-ready
#           data list for a hierarchical logistic regression of IT enrollment
#           on within-state poverty (CWC) and between-state tiered reimbursement
#           policy. Normalize survey weights, compute Kish DEFF diagnostics,
#           and save analysis-ready objects.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Section 5 (Application: NSECE 2019)
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
# Track   : A/B
# Inputs  :
#   Track A: data/restricted/cb_master_2019.rds (NSECE restricted-use data)
#   Track B: data/precomputed/application/stan_data.rds
#                                         analysis_data.rds
# Outputs :
#   data/precomputed/application/stan_data.rds     (Stan data list)
#   data/precomputed/application/analysis_data.rds (enriched data frame)
#   data/precomputed/application/data_meta.rds     (metadata and diagnostics)
#
# Data access:
#   This script requires NSECE 2019 restricted-use data from ICPSR Study 38445.
#   The data must be obtained through the NSECE data access procedures at:
#     https://www.icpsr.umich.edu/web/DSDR/studies/38445
#   See data/DATA_ACCESS.md for detailed instructions.
#
# Key processing steps:
#   1. Load analysis data (N = 6,785 providers, 51 states)
#   2. Group-mean center poverty (CWC) within each state
#   3. Verify tiered reimbursement is constant within each state
#   4. Build design matrix X = [intercept | poverty_cwc | tiered_reim]
#   5. Normalize survey weights: sum(w) = n_j per state
#   6. Compute Kish DEFF diagnostics (global and per-state)
#   7. Save Stan data list and enriched data frame
# =============================================================================


# =============================================================================
# Track selection
# =============================================================================
# Set USE_PRECOMPUTED = TRUE to skip restricted-data processing and load
# pre-computed results directly (Track B). Set FALSE for full pipeline (Track A).

USE_PRECOMPUTED <- TRUE


# =============================================================================
# Setup
# =============================================================================

cat("==============================================================\n")
cat("  Application Step 1: NSECE 2019 Data Preparation\n")
cat("==============================================================\n\n")

## Detect project root using here::here() or fallback
if (requireNamespace("here", quietly = TRUE)) {
  PROJECT_ROOT <- here::here()
} else {
  PROJECT_ROOT <- getwd()
}
cat(sprintf("[Setup] Project root: %s\n\n", PROJECT_ROOT))

## Define paths
DATA_DIR    <- file.path(PROJECT_ROOT, "data")
PRECOMP_DIR <- file.path(DATA_DIR, "precomputed", "application")

## Ensure output directory exists
if (!dir.exists(PRECOMP_DIR)) {
  dir.create(PRECOMP_DIR, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("[Setup] Created output directory: %s\n", PRECOMP_DIR))
}


# =============================================================================
# Track B: Load pre-computed results and exit early
# =============================================================================

if (USE_PRECOMPUTED) {

  cat("--- Track B: Loading pre-computed results ---\n\n")

  stan_data_path    <- file.path(PRECOMP_DIR, "stan_data.rds")
  analysis_path     <- file.path(PRECOMP_DIR, "analysis_data.rds")
  meta_path         <- file.path(PRECOMP_DIR, "data_meta.rds")

  ## Check that pre-computed files exist
  required_files <- c(stan_data_path, analysis_path, meta_path)
  missing_files  <- required_files[!file.exists(required_files)]

  if (length(missing_files) > 0) {
    stop(
      "Pre-computed files not found:\n",
      paste("  ", missing_files, collapse = "\n"),
      "\n\nTo use Track B, ensure these files are present in:\n  ",
      PRECOMP_DIR,
      "\n\nAlternatively, set USE_PRECOMPUTED <- FALSE and provide ",
      "restricted-use data (Track A)."
    )
  }

  stan_data     <- readRDS(stan_data_path)
  analysis_data <- readRDS(analysis_path)
  data_meta     <- readRDS(meta_path)

  ## Validate loaded data
  stopifnot(
    "N"     %in% names(stan_data), stan_data$N == 6785L,
    "J"     %in% names(stan_data), stan_data$J == 51L,
    "p"     %in% names(stan_data), stan_data$p == 3L,
    length(stan_data$y) == stan_data$N,
    nrow(stan_data$X)   == stan_data$N,
    ncol(stan_data$X)   == stan_data$p,
    length(stan_data$w) == stan_data$N
  )

  cat(sprintf("  stan_data: N = %d, J = %d, p = %d\n",
              stan_data$N, stan_data$J, stan_data$p))
  cat(sprintf("  Outcome prevalence: %.3f\n", mean(stan_data$y)))
  cat(sprintf("  Weight range: [%.4f, %.4f]\n",
              min(stan_data$w), max(stan_data$w)))
  cat(sprintf("  Kish DEFF (global): %.4f\n", data_meta$kish_deff_global))
  cat("\n  [Track B] Pre-computed data loaded successfully.\n")
  cat("  Skipping data processing pipeline.\n\n")
  cat("==============================================================\n")
  cat("  Step 1 Complete (Track B).\n")
  cat("==============================================================\n")

  ## Note: stan_data, analysis_data, data_meta are now in the global

  ## environment for downstream scripts.

} else {

# =============================================================================
# Track A: Full pipeline with restricted data
# =============================================================================

cat("--- Track A: Full data preparation pipeline ---\n\n")

## ---- Required packages ----
suppressPackageStartupMessages({
  library(dplyr)
})


## ---- Path to restricted-use analysis data ----
## This file is derived from the NSECE 2019 CB master data (ICPSR 38445)
## and contains the pre-processed analysis sample.
ANALYSIS_DATA_PATH <- file.path(DATA_DIR, "restricted", "analysis_data.rds")

if (!file.exists(ANALYSIS_DATA_PATH)) {
  stop(
    "Restricted-use analysis data not found at:\n  ", ANALYSIS_DATA_PATH,
    "\n\nThis script requires NSECE 2019 restricted-use data.",
    "\nSee data/DATA_ACCESS.md for data access instructions.",
    "\nFor partial replication, set USE_PRECOMPUTED <- TRUE."
  )
}


# =============================================================================
# Section 1: Load analysis data
# =============================================================================

cat("--- 1. Loading analysis data ---\n")

dat <- readRDS(ANALYSIS_DATA_PATH)

cat(sprintf("  Loaded: %d rows x %d columns\n", nrow(dat), ncol(dat)))

## ---- Define variable names ----
## Mapping from raw variable names to model roles:
##   within_var:  community poverty rate (standardized), varies within states
##   between_var: tiered reimbursement policy, constant within states
##   state_var:   integer state index (1..51)
##   weight_var:  survey weight (will be normalized)
##   outcome_var: binary IT enrollment indicator (0/1)
within_var   <- "comm_pct_poverty_num_std"
between_var  <- "ccdf_TieredReim"
state_var    <- "state_idx"
psu_var      <- "psu_idx"
weight_var   <- "w_tilde"
outcome_var  <- "z"
stratum_var  <- "stratum_idx"


# =============================================================================
# Section 2: Validate required columns
# =============================================================================

cat("\n--- 2. Validating columns ---\n")

required_cols <- c(within_var, between_var, state_var, psu_var,
                   weight_var, outcome_var, stratum_var, "state_name")

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Required columns not found: ", paste(missing_cols, collapse = ", "))
}

## Check for NAs in key columns
na_counts <- sapply(required_cols, function(col) sum(is.na(dat[[col]])))
if (any(na_counts > 0)) {
  na_report <- na_counts[na_counts > 0]
  stop("Missing values in:\n",
       paste(sprintf("  %s: %d NAs", names(na_report), na_report),
             collapse = "\n"))
}
cat("  [PASS] All required columns present, no missing values.\n")

## Extract key vectors
x_raw    <- dat[[within_var]]    # standardized poverty rate
z_btwn   <- dat[[between_var]]   # tiered reimbursement (0/1)
state    <- dat[[state_var]]     # state index 1..51
psu      <- dat[[psu_var]]       # PSU index 1..415
w_tilde  <- dat[[weight_var]]    # survey weights
y_out    <- dat[[outcome_var]]   # binary IT enrollment
stratum  <- dat[[stratum_var]]   # stratum index 1..30

N <- nrow(dat)
J <- length(unique(state))

## Validate indices, outcome, and weights
stopifnot(
  all(state >= 1L & state <= J),
  length(unique(state)) == J,
  all(y_out %in% c(0L, 1L, 0.0, 1.0)),
  all(w_tilde > 0)
)
y_out <- as.integer(y_out)

cat(sprintf("  N = %d observations, J = %d states\n", N, J))
cat(sprintf("  Outcome prevalence: %.3f\n", mean(y_out)))


# =============================================================================
# Section 3: Group-mean center poverty (CWC)
# =============================================================================
# Centering within context (CWC) removes between-state variation from the
# poverty covariate, creating a clean within/between separation. This maps
# directly onto the DER Tier I-a/I-b distinction (Section 3, Definition 2).
# After CWC, the poverty variable captures only within-state variation,
# making it identifiable as a "within-state" covariate whose DER should
# approximate the design effect (Theorem 1).

cat("\n--- 3. Group-mean centering poverty (CWC) ---\n")

state_means_x <- tapply(x_raw, state, mean)
x_cwc <- x_raw - state_means_x[state]

## Verification: state means of centered variable must be ~0
cwc_residuals <- tapply(x_cwc, state, mean)
max_cwc_dev   <- max(abs(cwc_residuals))
stopifnot(max_cwc_dev < 1e-10)

cat(sprintf("  CWC verified: max|state mean of x_cwc| = %.2e [PASS]\n",
            max_cwc_dev))
cat(sprintf("  x_cwc: mean = %.4f, sd = %.4f\n", mean(x_cwc), sd(x_cwc)))


# =============================================================================
# Section 4: Verify between-state covariate
# =============================================================================
# Tiered reimbursement is a state-level policy variable. It must be constant
# within each state. This is a between-state covariate whose information
# is absorbed into the between-group component (high R_k), yielding
# low DER (Tier I-b, protected).

cat("\n--- 4. Verifying tiered reimbursement is constant within states ---\n")

n_unique_between <- tapply(z_btwn, state, function(v) length(unique(v)))
states_with_variation <- which(n_unique_between > 1)

if (length(states_with_variation) > 0) {
  stop(sprintf("'%s' varies within %d state(s)", between_var,
               length(states_with_variation)))
}

z_between_state <- as.numeric(tapply(z_btwn, state, unique))
z_between_i     <- z_between_state[state]

cat(sprintf("  [PASS] Constant within all %d states.\n", J))
cat(sprintf("  Tiered: %d states = 1, %d states = 0\n",
            sum(z_between_state == 1), sum(z_between_state == 0)))


# =============================================================================
# Section 5: Build design matrix
# =============================================================================
# X has 3 columns matching the model in Equation (16) of the paper:
#   [1] intercept       : identifies beta_0
#   [2] poverty_cwc     : within-state covariate (CWC), identifies beta_1
#   [3] tiered_reim     : between-state policy, identifies beta_2

cat("\n--- 5. Building design matrix X (N x 3) ---\n")

X <- cbind(
  intercept   = rep(1.0, N),
  poverty_cwc = x_cwc,
  tiered_reim = z_between_i
)

stopifnot(nrow(X) == N, ncol(X) == 3L, all(X[, 1] == 1.0))

cat(sprintf("  X dimensions: %d x %d\n", nrow(X), ncol(X)))
cat(sprintf("  poverty_cwc: mean = %.4f, sd = %.4f\n",
            mean(X[, 2]), sd(X[, 2])))
cat(sprintf("  tiered_reim: mean = %.4f (proportion 1s)\n", mean(X[, 3])))


# =============================================================================
# Section 6: Normalize survey weights
# =============================================================================
# Survey weights are normalized within each state so that the sum of weights
# equals the state sample size n_j. This "sum-to-N" normalization is the
# standard pseudo-likelihood scaling for Bayesian survey models
# (Savitsky & Toth, 2016; Williams & Savitsky, 2021).
#
# Normalization formula: w_norm_i = w_tilde_i * n_j / sum(w_tilde in state j)
#
# After normalization, the weights capture relative importance within each
# state but do not inflate or deflate the effective sample size.

cat("\n--- 6. Normalizing survey weights (sum-to-N per state) ---\n")

n_j_vec      <- as.integer(table(state))
w_norm       <- numeric(N)
state_wt_sums <- tapply(w_tilde, state, sum)

for (j in seq_len(J)) {
  idx_j <- which(state == j)
  if (state_wt_sums[j] < .Machine$double.eps) {
    stop(sprintf("State %d has zero total weight.", j))
  }
  w_norm[idx_j] <- w_tilde[idx_j] * n_j_vec[j] / state_wt_sums[j]
}

## Verification: within-state sums must equal n_j
w_norm_state_sums <- tapply(w_norm, state, sum)
max_wt_dev <- max(abs(w_norm_state_sums - n_j_vec))
stopifnot(max_wt_dev < 1e-8)

cat(sprintf("  Normalization verified: max|sum(w_norm) - n_j| = %.2e [PASS]\n",
            max_wt_dev))
cat(sprintf("  w_norm: mean = %.4f, range = [%.4f, %.4f]\n",
            mean(w_norm), min(w_norm), max(w_norm)))


# =============================================================================
# Section 7: Kish DEFF diagnostics
# =============================================================================
# The Kish (1965) design effect measures survey weight variability:
#   DEFF_kish = n * sum(w^2) / (sum(w))^2
# This provides a baseline for interpreting DER values. For within-state
# covariates, the DER should approximate DEFF_kish (Theorem 1).

cat("\n--- 7. Computing Kish DEFF diagnostics ---\n")

kish_deff <- function(w) {
  n   <- length(w)
  sw  <- sum(w)
  sw2 <- sum(w^2)
  n * sw2 / sw^2
}

## Global DEFF (using raw weights)
deff_global <- kish_deff(w_tilde)

## Per-state DEFF
deff_state <- vapply(seq_len(J), function(j) {
  w_j <- w_tilde[state == j]
  if (length(w_j) < 2L) return(1.0)
  kish_deff(w_j)
}, numeric(1))

## Weight CV
w_cv_global <- sd(w_tilde) / mean(w_tilde)

cat(sprintf("  Global Kish DEFF: %.4f\n", deff_global))
cat(sprintf("  Global weight CV: %.4f\n", w_cv_global))
cat(sprintf("  Per-state DEFF: mean = %.4f, range = [%.4f, %.4f]\n",
            mean(deff_state), min(deff_state), max(deff_state)))


# =============================================================================
# Section 8: Assemble Stan data list
# =============================================================================

cat("\n--- 8. Assembling Stan data list ---\n")

## State name mapping
state_names_map <- tapply(as.character(dat$state_name), state,
                          function(v) v[1])

## Stan data list: matches the interface of hlr_weighted.stan
stan_data <- list(
  N     = N,                    # total observations (6785)
  J     = J,                    # number of states (51)
  p     = 3L,                   # fixed-effect predictors (intercept + 2)
  y     = as.integer(y_out),    # binary outcome
  X     = X,                    # N x 3 design matrix
  group = as.integer(state),    # state membership (1..51)
  w     = w_norm                # normalized survey weights
)

## Metadata for downstream analysis
data_meta <- list(
  psu_idx          = as.integer(psu),
  stratum_idx      = as.integer(stratum),
  n_psu            = length(unique(psu)),
  n_strata         = length(unique(stratum)),
  state_n          = setNames(n_j_vec, seq_len(J)),
  state_names      = state_names_map,
  kish_deff_global = deff_global,
  kish_deff_state  = setNames(as.numeric(deff_state), seq_len(J)),
  covariate_names  = c("intercept", "poverty_cwc", "tiered_reim"),
  param_types      = c("fe_between", "fe_within", "fe_between"),
  w_cv_global      = w_cv_global
)

## Validation
stopifnot(
  stan_data$N == 6785L,
  stan_data$J == 51L,
  stan_data$p == 3L,
  length(stan_data$y) == stan_data$N,
  nrow(stan_data$X)   == stan_data$N,
  ncol(stan_data$X)   == stan_data$p,
  all(stan_data$y %in% c(0L, 1L)),
  all(stan_data$group >= 1L & stan_data$group <= stan_data$J),
  all(stan_data$w > 0),
  all(stan_data$X[, 1] == 1.0)
)
cat("  [PASS] Stan data structure validated.\n")


# =============================================================================
# Section 9: Build enriched analysis data frame
# =============================================================================

cat("\n--- 9. Building enriched analysis data frame ---\n")

analysis_data <- dat
analysis_data$poverty_cwc <- X[, 2]
analysis_data$w_norm      <- w_norm
analysis_data$tiered_reim <- X[, 3]

cat(sprintf("  Enriched data frame: %d rows x %d columns\n",
            nrow(analysis_data), ncol(analysis_data)))


# =============================================================================
# Section 10: Save outputs
# =============================================================================

cat("\n--- 10. Saving outputs ---\n")

## Stan data
stan_data_path <- file.path(PRECOMP_DIR, "stan_data.rds")
saveRDS(stan_data, stan_data_path)
cat(sprintf("  Saved: %s (%.1f KB)\n", stan_data_path,
            file.info(stan_data_path)$size / 1024))

## Enriched analysis data
analysis_path <- file.path(PRECOMP_DIR, "analysis_data.rds")
saveRDS(analysis_data, analysis_path)
cat(sprintf("  Saved: %s (%.1f KB)\n", analysis_path,
            file.info(analysis_path)$size / 1024))

## Metadata
meta_path <- file.path(PRECOMP_DIR, "data_meta.rds")
saveRDS(data_meta, meta_path)
cat(sprintf("  Saved: %s (%.1f KB)\n", meta_path,
            file.info(meta_path)$size / 1024))


# =============================================================================
# Summary
# =============================================================================

cat("\n==============================================================\n")
cat("  DATA PREPARATION SUMMARY\n")
cat("==============================================================\n")
cat(sprintf("  N = %d providers, J = %d states, p = 3 fixed effects\n", N, J))
cat(sprintf("  Outcome: IT enrollment (z), prevalence = %.3f\n", mean(y_out)))
cat(sprintf("  Covariates:\n"))
cat(sprintf("    [1] intercept\n"))
cat(sprintf("    [2] poverty_cwc  (within-state, CWC)\n"))
cat(sprintf("    [3] tiered_reim  (between-state policy)\n"))
cat(sprintf("  Global Kish DEFF: %.4f\n", deff_global))
cat(sprintf("  State sizes: min = %d, median = %d, max = %d\n",
            min(n_j_vec), as.integer(median(n_j_vec)), max(n_j_vec)))
cat("==============================================================\n")
cat("  Step 1 Complete (Track A).\n")
cat("==============================================================\n")

}  # end Track A
