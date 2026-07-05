# =============================================================================
# app_01_prepare_data.R -- Build the model-ready NSECE arrays (Track A)
# Author : JoonHo Lee (jlee296@ua.edu)
# =============================================================================
# Turn a restricted-use NSECE 2019 analysis extract into the two objects the
# fitting and diagnostic scripts consume:
#
#   phase3_data.rds     a list with the model arrays (N, J, p, y, X, group, w);
#                       here w is the within-state sensitivity convention
#   phase3_analysis.rds an enriched data frame carrying the raw design weight and
#                       the PSU / stratum indices for the design-based target
#
# The derivations are: group-mean-centering the poverty covariate within each
# state (so it is a purely within-state predictor), verifying that the tiered
# reimbursement indicator is constant within states (a between-state predictor),
# assembling the design matrix X = [intercept, poverty_cwc, tiered_reim], and
# forming the within-state pseudo-likelihood weights used by the supplement.
# The headline v2 fit rebuilds global unit-mean weights from
# phase3_analysis$weight in app_02_fit_refit.R.
#
# This is a Track A step and needs the restricted data (see data/DATA_ACCESS.md).
# It writes only into the restricted directory NSECE_RESTRICTED_DIR (default the
# git-ignored data/nsece/); nothing here is shipped with the package. Readers
# reproducing the figures and tables (Track B) do not run this script.
#
# The input extract is expected as an RDS data frame (one row per provider) with
# the columns mapped below. Adapt the column names to your own NSECE extract.
#
# Usage : Rscript code/02_application/app_01_prepare_data.R
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

here_root <- tryCatch(here::here(), error = function(e) getwd())
restricted_dir <- Sys.getenv("NSECE_RESTRICTED_DIR",
                             file.path(here_root, "data", "nsece"))

# The provider-level analysis extract derived from the NSECE restricted files.
extract_path <- file.path(restricted_dir, "nsece_analysis_extract.rds")
if (!file.exists(extract_path)) {
  stop("NSECE analysis extract not found at '", extract_path, "'. This is a ",
       "Track A step; see data/DATA_ACCESS.md. Track B users work from the ",
       "shipped data/precomputed/application/ results and do not run this.",
       call. = FALSE)
}

# ---- Column roles (adapt to your extract) -----------------------------------
within_var  <- "comm_pct_poverty_num_std"  # community poverty rate (standardized)
between_var <- "ccdf_TieredReim"            # tiered reimbursement policy (0/1)
state_var   <- "state_idx"                  # state index 1..J
psu_var     <- "psu_idx"                    # design PSU index
weight_var  <- "weight"                     # raw design weight
outcome_var <- "z"                          # binary Infant-Toddler participation
stratum_var <- "stratum_idx"               # design stratum index

dat <- readRDS(extract_path)

required_cols <- c(within_var, between_var, state_var, psu_var,
                   weight_var, outcome_var, stratum_var, "state_name")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Required columns not found in the extract: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}
if (any(vapply(dat[required_cols], anyNA, logical(1)))) {
  stop("Missing values in one or more required columns.", call. = FALSE)
}

x_raw   <- dat[[within_var]]
z_btwn  <- dat[[between_var]]
state   <- as.integer(dat[[state_var]])
psu     <- as.integer(dat[[psu_var]])
w_raw   <- dat[[weight_var]]
y_out   <- as.integer(dat[[outcome_var]])
stratum <- as.integer(dat[[stratum_var]])

N <- nrow(dat)
J <- length(unique(state))
stopifnot(all(state >= 1L & state <= J), length(unique(state)) == J,
          all(y_out %in% c(0L, 1L)), all(w_raw > 0))

# ---- Group-mean-center the poverty covariate within each state --------------
# Centering within context removes the between-state variation, so the covariate
# carries only within-state contrasts. This is the within-state predictor whose
# DER should approximate the design effect.
state_means_x <- tapply(x_raw, state, mean)
poverty_cwc   <- x_raw - state_means_x[state]
stopifnot(max(abs(tapply(poverty_cwc, state, mean))) < 1e-10)

# ---- Verify the tiered reimbursement indicator is a state-level covariate ---
n_unique_between <- tapply(z_btwn, state, function(v) length(unique(v)))
if (any(n_unique_between > 1)) {
  stop("'", between_var, "' varies within some states; it must be state-level.",
       call. = FALSE)
}
tiered_state <- as.numeric(tapply(z_btwn, state, unique))
tiered_reim  <- tiered_state[state]

# ---- Design matrix X = [intercept, poverty_cwc, tiered_reim] ----------------
X <- cbind(intercept   = rep(1.0, N),
           poverty_cwc = poverty_cwc,
           tiered_reim = tiered_reim)
stopifnot(nrow(X) == N, ncol(X) == 3L, all(X[, 1] == 1.0))

# ---- Within-state pseudo-likelihood weights (sum to n_j per state) ----------
# The raw design weight is preserved in phase3_analysis for the headline global
# unit-mean convention. phase3_data$w stores the within-state alternative used
# by app_03_der_supplement.R; app_02_fit_refit.R does not use it for the v2 fit.
n_j_vec <- as.integer(table(state))
state_wt_sums <- tapply(w_raw, state, sum)
w_within <- numeric(N)
for (j in seq_len(J)) {
  idx_j <- which(state == j)
  w_within[idx_j] <- w_raw[idx_j] * n_j_vec[j] / state_wt_sums[j]
}
stopifnot(max(abs(tapply(w_within, state, sum) - n_j_vec)) < 1e-8)

# ---- Kish design-effect diagnostics -----------------------------------------
kish_deff <- function(w) length(w) * sum(w^2) / sum(w)^2
deff_global <- kish_deff(w_raw)

# ---- Assemble and save the two objects --------------------------------------
phase3_data <- list(
  N = N, J = J, p = 3L,
  y = y_out, X = X,
  group = state,
  w = w_within                 # within-state pseudo-likelihood weights
)

phase3_analysis <- dat
phase3_analysis$state_idx   <- state
phase3_analysis$psu_idx     <- psu
phase3_analysis$stratum_idx <- stratum
phase3_analysis$weight      <- w_raw           # raw design weight
phase3_analysis$poverty_cwc <- poverty_cwc
phase3_analysis$tiered_reim <- tiered_reim
phase3_analysis$y           <- y_out

stopifnot(phase3_data$N == 6785L, phase3_data$J == 51L, phase3_data$p == 3L,
          all(phase3_data$group == phase3_analysis$state_idx))

dir.create(restricted_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(phase3_data,     file.path(restricted_dir, "phase3_data.rds"))
saveRDS(phase3_analysis, file.path(restricted_dir, "phase3_analysis.rds"))

cat(sprintf("N = %d providers, J = %d states, %d PSUs in %d strata\n",
            N, J, length(unique(psu)), length(unique(stratum))))
cat(sprintf("Outcome prevalence: %.3f | global Kish DEFF: %.3f\n",
            mean(y_out), deff_global))
cat(sprintf("Wrote phase3_data.rds and phase3_analysis.rds to %s\n",
            restricted_dir))
