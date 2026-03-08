# =============================================================================
# app_03_der_analysis.R: DER Computation, Classification, and Correction
# =============================================================================
#
# Purpose : Compute Design Effect Ratios (DER) for all 54 parameters of the
#           NSECE hierarchical logistic regression. Implements the full
#           compute-classify-correct (CCC) workflow from Algorithm 1.
#           Performs three-tier classification and selective Cholesky
#           correction for flagged parameters.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Section 5 (Application: NSECE 2019)
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
# Track   : A/B
# Inputs  :
#   Track A: data/precomputed/application/stan_data.rds      (from app_01)
#            data/precomputed/application/analysis_data.rds   (from app_01)
#            data/precomputed/application/data_meta.rds       (from app_01)
#            data/precomputed/application/model_summary.rds   (from app_02)
#   Track B: data/precomputed/application/der_results.rds
# Outputs :
#   data/precomputed/application/der_results.rds
#     Contains: H_obs, H_obs_inv, J_cluster, V_sand, DER values,
#               classification, selective/blanket correction results,
#               Sigma_MCMC, per-state DEFF, shrinkage factors, etc.
#
# DER methodology (Definition 1):
#   DER_k = diag(V_sand)_k / diag(Sigma_MCMC)_k
#
# where:
#   V_sand    = H_obs^{-1} J_cluster H_obs^{-1}  (sandwich variance)
#   H_obs     = observed information matrix (negative Hessian + prior)
#   J_cluster = sum_j s_j s_j^T  (clustered score outer product)
#   Sigma_MCMC = empirical covariance of posterior draws
#
# The sandwich variance captures the impact of complex survey design
# (stratification, clustering, unequal weights) on parameter uncertainty.
# When DER_k > 1, the model-based posterior underestimates uncertainty;
# when DER_k < 1, hierarchical structure provides variance reduction
# beyond what the survey design would suggest.
#
# Classification (tau = 1.2):
#   DER > tau   -> flag for correction
#   DER <= tau  -> retain original posterior
#
# Selective correction (Algorithm 1, Step 3):
#   For flagged parameter k:
#     phi_corrected^(m) = phi_hat + sqrt(DER_k) * (phi^(m) - phi_hat)
# =============================================================================


# =============================================================================
# Track selection
# =============================================================================

USE_PRECOMPUTED <- TRUE


# =============================================================================
# Setup
# =============================================================================

cat("==============================================================\n")
cat("  Application Step 3: DER Analysis\n")
cat("==============================================================\n\n")

if (requireNamespace("here", quietly = TRUE)) {
  PROJECT_ROOT <- here::here()
} else {
  PROJECT_ROOT <- getwd()
}
cat(sprintf("[Setup] Project root: %s\n\n", PROJECT_ROOT))

PRECOMP_DIR <- file.path(PROJECT_ROOT, "data", "precomputed", "application")
HELPERS_DIR <- file.path(PROJECT_ROOT, "code", "helpers")

## Fixed constants (must match app_01 and app_02)
P_FIXED        <- 3L          # number of fixed effects
J_STATES       <- 51L         # number of states
N_OBS          <- 6785L       # total observations
D_TOTAL        <- P_FIXED + J_STATES  # = 54 total parameters
BETA_PRIOR_SD  <- 5.0         # prior SD for beta: N(0, 25)
TAU_THRESHOLD  <- 1.2         # classification threshold


# =============================================================================
# Track B: Load pre-computed results
# =============================================================================

if (USE_PRECOMPUTED) {

  cat("--- Track B: Loading pre-computed DER results ---\n\n")

  der_results_path <- file.path(PRECOMP_DIR, "der_results.rds")

  if (!file.exists(der_results_path)) {
    stop(
      "Pre-computed DER results not found:\n  ", der_results_path,
      "\n\nTo use Track B, ensure this file is present.",
      "\nSet USE_PRECOMPUTED <- FALSE for Track A."
    )
  }

  der_results <- readRDS(der_results_path)

  ## Validate essential fields
  expected_fields <- c("der", "der_beta", "der_theta",
                       "classification", "comparison",
                       "V_sand", "sigma_mcmc",
                       "beta_hat", "theta_hat")
  missing_fields <- setdiff(expected_fields, names(der_results))
  if (length(missing_fields) > 0) {
    stop("der_results.rds missing fields: ",
         paste(missing_fields, collapse = ", "))
  }

  ## Print summary
  cat("  DER Summary (54 parameters):\n")
  cat(sprintf("    Fixed effects:\n"))
  cov_names <- c("intercept", "poverty_cwc", "tiered_reim")
  for (k in seq_len(P_FIXED)) {
    cat(sprintf("      %-14s  DER = %.4f\n",
                cov_names[k], der_results$der_beta[k]))
  }
  cat(sprintf("    Random effects:\n"))
  cat(sprintf("      Mean DER:  %.4f\n", mean(der_results$der_theta)))
  cat(sprintf("      Range:     [%.4f, %.4f]\n",
              min(der_results$der_theta), max(der_results$der_theta)))

  n_flagged <- sum(der_results$classification$flagged)
  cat(sprintf("\n    Flagged (DER > %.1f): %d / %d parameters\n",
              TAU_THRESHOLD, n_flagged, D_TOTAL))

  cat("\n  [Track B] Pre-computed DER results loaded successfully.\n")
  cat("==============================================================\n")
  cat("  Step 3 Complete (Track B).\n")
  cat("==============================================================\n")

} else {

# =============================================================================
# Track A: Full DER analysis pipeline
# =============================================================================

cat("--- Track A: Full DER analysis pipeline ---\n\n")

suppressPackageStartupMessages({
  library(Matrix)
})


# =============================================================================
# Section 1: Load input data
# =============================================================================

cat("--- 1. Loading input data ---\n")

## 1a. Stan data list
stan_data_path <- file.path(PRECOMP_DIR, "stan_data.rds")
stopifnot(file.exists(stan_data_path))
stan_data <- readRDS(stan_data_path)
stopifnot(stan_data$N == N_OBS, stan_data$J == J_STATES, stan_data$p == P_FIXED)

y     <- as.integer(stan_data$y)
X     <- as.matrix(stan_data$X)
group <- as.integer(stan_data$group)
w     <- as.numeric(stan_data$w)
N     <- stan_data$N
J     <- stan_data$J
p     <- stan_data$p

cat(sprintf("  Stan data: N = %d, J = %d, p = %d\n", N, J, p))

## 1b. Analysis data (for PSU indices)
analysis_path <- file.path(PRECOMP_DIR, "analysis_data.rds")
stopifnot(file.exists(analysis_path))
analysis_df <- readRDS(analysis_path)
psu_idx <- as.integer(analysis_df$psu_idx)
stopifnot(length(psu_idx) == N)
G_psu <- max(psu_idx)
cat(sprintf("  PSU indices: G = %d PSUs\n", G_psu))

## 1c. Data metadata
meta_path <- file.path(PRECOMP_DIR, "data_meta.rds")
stopifnot(file.exists(meta_path))
data_meta <- readRDS(meta_path)

## 1d. Model summary (posterior draws)
summary_path <- file.path(PRECOMP_DIR, "model_summary.rds")
stopifnot(file.exists(summary_path))
model_summary <- readRDS(summary_path)

summ_df     <- model_summary$summary
draws_beta  <- model_summary$draws_beta    # [M x 3]
draws_sigma <- model_summary$draws_sigma   # [M]
draws_theta <- model_summary$draws_theta   # [M x 51]

stopifnot(
  ncol(draws_beta) == P_FIXED,
  ncol(draws_theta) == J_STATES,
  nrow(draws_beta) == nrow(draws_theta)
)

M <- nrow(draws_beta)
cat(sprintf("  Posterior draws: M = %d\n\n", M))


# =============================================================================
# Section 2: Extract posterior means
# =============================================================================

cat("--- 2. Extracting posterior means ---\n")

## Helper: extract posterior mean from summary
get_mean <- function(par_name) {
  idx <- which(summ_df$variable == par_name)
  if (length(idx) == 0L) stop("Parameter '", par_name, "' not found.")
  summ_df$mean[idx[1L]]
}

## Fixed effects
beta_hat <- vapply(seq_len(P_FIXED),
                   function(k) get_mean(paste0("beta[", k, "]")),
                   numeric(1L))
names(beta_hat) <- c("intercept", "poverty_cwc", "tiered_reim")

## Random effects
theta_hat <- vapply(seq_len(J_STATES),
                    function(j) get_mean(paste0("theta[", j, "]")),
                    numeric(1L))

## sigma_theta
sigma_theta_hat <- get_mean("sigma_theta")
stopifnot(sigma_theta_hat > 0)

cat(sprintf("  beta_hat:        [%.4f, %.4f, %.4f]\n",
            beta_hat[1], beta_hat[2], beta_hat[3]))
cat(sprintf("  sigma_theta_hat: %.4f\n", sigma_theta_hat))
cat(sprintf("  theta_hat:       mean = %.4f, sd = %.4f\n\n",
            mean(theta_hat), sd(theta_hat)))


# =============================================================================
# Section 3: Predicted probabilities and working weights
# =============================================================================
# The observed information matrix uses working weights from the logistic
# model: wt_i = q_i * (1 - q_i), where q_i is the predicted probability.
# These are the diagonal entries of the weight matrix W in the Fisher
# information for logistic regression.

cat("--- 3. Computing predicted probabilities ---\n")

## Linear predictor: eta_i = X_i beta + theta_{group[i]}
eta <- as.numeric(X %*% beta_hat) + theta_hat[group]

## Predicted probabilities (inverse logit)
q_hat <- plogis(eta)

## Working weights (Bernoulli variance function)
wt <- q_hat * (1.0 - q_hat)

cat(sprintf("  q_hat: mean = %.4f, range = [%.4f, %.4f]\n",
            mean(q_hat), min(q_hat), max(q_hat)))
cat(sprintf("  wt:    mean = %.4f, range = [%.6f, %.4f]\n\n",
            mean(wt), min(wt), max(wt)))


# =============================================================================
# Section 4: Observed information matrix H_obs
# =============================================================================
# H_obs is the (p + J) x (p + J) observed information matrix, which is the
# negative Hessian of the log pseudo-posterior evaluated at the posterior mean.
# It combines data likelihood information with prior information:
#
#   H_obs = H_data + H_prior
#
# For the hierarchical logistic regression:
#   H_data has a block structure:
#     H_bb = X^T diag(v) X          (p x p, beta-beta block)
#     H_bt = X^T diag(v) Z          (p x J, beta-theta block)
#     H_tt = Z^T diag(v) Z + D_RE   (J x J, theta-theta block)
#   where v_i = w_i * wt_i (survey weight * working weight)
#   and D_RE = diag(1/sigma_theta^2) is the random-effect precision.
#
# The prior contribution adds diag(1/beta_prior_sd^2) to H_bb.

cat("--- 4. Building observed information matrix H_obs ---\n")

d <- P_FIXED + J_STATES  # = 54

## Combined weight: survey weight * working weight
v <- w * wt

## Initialize H_obs
H_obs <- matrix(0.0, d, d)

## Beta-beta block (top-left p x p): X^T diag(v) X + prior precision
H_bb <- crossprod(X, X * v)
## Add prior contribution: beta ~ N(0, beta_prior_sd^2)
H_bb <- H_bb + diag(1.0 / BETA_PRIOR_SD^2, P_FIXED)
H_obs[1:p, 1:p] <- H_bb

## Beta-theta and theta-theta blocks
## For each state j, accumulate the v_i-weighted contributions
for (j in seq_len(J)) {
  idx_j <- which(group == j)
  v_j   <- v[idx_j]
  X_j   <- X[idx_j, , drop = FALSE]

  ## sum_j = sum of v_i for observations in state j
  sum_v_j <- sum(v_j)

  ## Beta-theta block (row: beta, col: theta_j)
  ## H_bt[, j] = sum_{i in j} v_i * X_i
  bt_j <- colSums(X_j * v_j)
  H_obs[1:p, p + j] <- bt_j
  H_obs[p + j, 1:p] <- bt_j  # symmetric

  ## Theta-theta diagonal: sum of v_i in state j + 1/sigma_theta^2
  H_obs[p + j, p + j] <- sum_v_j + 1.0 / sigma_theta_hat^2
}

## Enforce exact symmetry
H_obs <- (H_obs + t(H_obs)) / 2.0

## Verify properties
cat(sprintf("  H_obs: %d x %d\n", nrow(H_obs), ncol(H_obs)))
sym_err <- max(abs(H_obs - t(H_obs)))
cat(sprintf("  Symmetry: max|H - H^T| = %.2e [%s]\n",
            sym_err, ifelse(sym_err < 1e-10, "PASS", "WARNING")))

## Invert H_obs (with nearPD fallback for numerical stability)
H_obs_inv <- tryCatch(
  solve(H_obs),
  error = function(e) {
    warning("H_obs is near-singular. Applying nearPD fallback.")
    M_pd <- as.matrix(nearPD(H_obs, keepDiag = TRUE)$mat)
    solve(M_pd)
  }
)

## Verify inversion quality
resid_inv <- max(abs(H_obs %*% H_obs_inv - diag(d)))
cat(sprintf("  Inversion: max|H * H_inv - I| = %.2e [%s]\n",
            resid_inv, ifelse(resid_inv < 1e-8, "PASS", "WARNING")))
cat(sprintf("  H_obs_inv diag range: [%.4e, %.4e]\n\n",
            min(diag(H_obs_inv)), max(diag(H_obs_inv))))


# =============================================================================
# Section 5: Weighted score residuals
# =============================================================================
# Score residuals are the derivatives of the weighted pseudo-log-likelihood
# with respect to the linear predictor. For Bernoulli with logit link:
#   r_i = w_i * (y_i - q_i)
# These are aggregated by cluster to form the clustered score outer product.

cat("--- 5. Computing weighted score residuals ---\n")

r <- w * (y - q_hat)

cat(sprintf("  r: mean = %.6f (should be ~0), sd = %.4f\n\n",
            mean(r), sd(r)))


# =============================================================================
# Section 6: Clustered score outer product J_cluster
# =============================================================================
# J_cluster = sum_{j=1}^{J} s_j s_j^T is the "meat" of the sandwich
# estimator, where s_j is the (p + J)-dimensional score vector for cluster j.
# The score vector has two parts:
#   s_j[1:p]   = sum_{i in j} r_i * X_i     (beta component)
#   s_j[p + j'] = sum_{i in j, group=j'} r_i (theta component)
#
# When clustering matches the model grouping (state-level), the theta
# component is non-zero only for j' = j.

cat("--- 6. Building clustered score outer product J_cluster ---\n")

J_cluster <- matrix(0.0, d, d)

for (j in seq_len(J)) {
  idx_j <- which(group == j)
  if (length(idx_j) == 0L) next

  s_j <- numeric(d)

  ## Beta components
  if (length(idx_j) == 1L) {
    s_j[1:p] <- r[idx_j] * X[idx_j, ]
  } else {
    s_j[1:p] <- colSums(X[idx_j, , drop = FALSE] * r[idx_j])
  }

  ## Theta component (only entry j, since clustering = grouping)
  s_j[p + j] <- sum(r[idx_j])

  ## Outer product contribution
  J_cluster <- J_cluster + tcrossprod(s_j)
}

## Enforce symmetry
J_cluster <- (J_cluster + t(J_cluster)) / 2.0

cat(sprintf("  J_cluster: %d x %d\n", nrow(J_cluster), ncol(J_cluster)))
cat(sprintf("  Symmetry: max|J - J^T| = %.2e\n\n",
            max(abs(J_cluster - t(J_cluster)))))


# =============================================================================
# Section 7: Sandwich variance
# =============================================================================
# V_sand = H_obs^{-1} J_cluster H_obs^{-1}
# This is the survey-design-adjusted variance for the full parameter vector.
# It accounts for clustering, stratification (through score residuals),
# and unequal weighting.

cat("--- 7. Computing sandwich variance V_sand ---\n")

V_sand <- H_obs_inv %*% J_cluster %*% H_obs_inv
V_sand <- (V_sand + t(V_sand)) / 2.0  # enforce symmetry

diag_Vs <- diag(V_sand)

if (any(diag_Vs <= 0)) {
  warning(sprintf("V_sand has %d non-positive diagonal entries.", sum(diag_Vs <= 0)))
} else {
  cat("  V_sand: all diagonal entries positive [PASS]\n")
}
cat(sprintf("  V_sand diag range: [%.4e, %.4e]\n\n",
            min(diag_Vs), max(diag_Vs)))


# =============================================================================
# Section 8: MCMC posterior covariance
# =============================================================================
# Sigma_MCMC is the empirical covariance of the MCMC draws, which serves as
# the denominator of the DER (Definition 1). This is the model-based posterior
# variance that does NOT account for survey design effects.

cat("--- 8. Computing MCMC posterior covariance ---\n")

## Combine beta and theta draws: [M x d]
draws_all <- cbind(draws_beta, draws_theta)
stopifnot(ncol(draws_all) == d)

sigma_mcmc <- cov(draws_all)
sigma_mcmc <- (sigma_mcmc + t(sigma_mcmc)) / 2.0

diag_mcmc <- diag(sigma_mcmc)
stopifnot(all(diag_mcmc > 0))

cat(sprintf("  Sigma_MCMC: %d x %d, all diag > 0 [PASS]\n", d, d))
cat(sprintf("  Sigma_MCMC diag range: [%.4e, %.4e]\n\n",
            min(diag_mcmc), max(diag_mcmc)))


# =============================================================================
# Section 9: Design Effect Ratios
# =============================================================================
# DER_k = diag(V_sand)_k / diag(Sigma_MCMC)_k
#
# Interpretation (Definition 1):
#   DER > 1 : model-based posterior underestimates uncertainty (design-sensitive)
#   DER = 1 : no design effect on this parameter
#   DER < 1 : hierarchical structure absorbs design effect (protected)

cat("--- 9. Computing Design Effect Ratios ---\n")

## Parameter names
param_names <- c(
  paste0("beta[", seq_len(P_FIXED), "]"),
  paste0("theta[", seq_len(J_STATES), "]")
)

## DER values
der_all        <- diag_Vs / diag_mcmc
names(der_all) <- param_names

der_beta  <- der_all[1:P_FIXED]
der_theta <- der_all[(P_FIXED + 1):D_TOTAL]

## Laplace-based DER (cross-check: sandwich vs model information)
diag_H_inv     <- diag(H_obs_inv)
der_laplace    <- diag_Vs / diag_H_inv
names(der_laplace) <- param_names

covariate_names <- c("intercept", "poverty_cwc", "tiered_reim")

cat("  Fixed-effect DER:\n")
for (k in seq_len(P_FIXED)) {
  cat(sprintf("    %-14s DER = %.4f\n", covariate_names[k], der_beta[k]))
}
cat(sprintf("\n  Random-effect DER:\n"))
cat(sprintf("    Mean: %.4f, Range: [%.4f, %.4f]\n",
            mean(der_theta), min(der_theta), max(der_theta)))
cat(sprintf("    RE with DER < 1: %d / %d (%.1f%%)\n\n",
            sum(der_theta < 1), J_STATES, 100 * sum(der_theta < 1) / J_STATES))


# =============================================================================
# Section 10: Per-state design quantities
# =============================================================================
# Per-state Kish DEFF and shrinkage factors B_j are computed for the
# decomposition theorem (Theorem 2):
#   DER_j ~ B_j * DEFF_j * kappa_j
#
# B_j = n_j^eff * sigma_theta^2 / (1 + n_j^eff * sigma_theta^2)
# is the shrinkage factor measuring how much state j's random effect
# is determined by its own data versus the prior.

cat("--- 10. Computing per-state DEFF and shrinkage factors ---\n")

## Per-state Kish DEFF (computed on normalized weights within states)
deff_j <- vapply(seq_len(J), function(j) {
  w_j <- w[group == j]
  n_j <- length(w_j)
  if (n_j < 2L) return(1.0)
  n_j * sum(w_j^2) / sum(w_j)^2
}, numeric(1))

## Shrinkage factors B_j
## B_j measures the posterior precision attributable to data vs prior
## For each state j:
##   n_j_eff = sum(v_j) (effective sample size in terms of information)
##   B_j = n_j_eff * sigma_theta^2 / (1 + n_j_eff * sigma_theta^2)
B_j <- vapply(seq_len(J), function(j) {
  idx_j   <- which(group == j)
  n_j_eff <- sum(v[idx_j])  # effective information from data
  n_j_eff * sigma_theta_hat^2 / (1 + n_j_eff * sigma_theta_hat^2)
}, numeric(1))

cat(sprintf("  deff_j: mean = %.4f, range = [%.4f, %.4f]\n",
            mean(deff_j), min(deff_j), max(deff_j)))
cat(sprintf("  B_j:    mean = %.4f, range = [%.4f, %.4f]\n\n",
            mean(B_j), min(B_j), max(B_j)))


# =============================================================================
# Section 11: Three-tier classification
# =============================================================================
# Classification scheme (Definition 2):
#   Tier I-a  (Survey-dominated):  within-state FE parameters
#   Tier I-b  (Protected, between): between-state FE parameters
#   Tier II   (Protected, random):  random effects
#
# Correction decision:
#   DER > tau (= 1.2) -> CORRECT (flag for selective adjustment)
#   DER <= tau         -> RETAIN  (keep original posterior)

cat("--- 11. Three-tier classification (tau = %.1f) ---\n", TAU_THRESHOLD)

param_types <- c("fe_between", "fe_within", "fe_between",
                 rep("re", J_STATES))

## Build classification data frame
classification <- data.frame(
  param_name = param_names,
  param_type = param_types,
  der        = as.numeric(der_all),
  stringsAsFactors = FALSE
)

classification$tier <- ifelse(
  param_types == "fe_within",  "I-a",
  ifelse(param_types == "fe_between", "I-b",
         ifelse(param_types == "re", "II", "unknown")))

classification$tier_label <- ifelse(
  classification$tier == "I-a",  "Survey-dominated",
  ifelse(classification$tier == "I-b", "Protected (between)",
         ifelse(classification$tier == "II", "Protected (random effects)",
                "Unknown")))

classification$flagged <- classification$der > TAU_THRESHOLD
classification$action  <- ifelse(classification$flagged, "CORRECT", "retain")

## Report
n_flagged <- sum(classification$flagged)
cat(sprintf("  Flagged (DER > %.1f): %d / %d parameters (%.1f%%)\n",
            TAU_THRESHOLD, n_flagged, D_TOTAL, 100 * n_flagged / D_TOTAL))

if (n_flagged > 0) {
  flagged_rows <- classification[classification$flagged, ]
  cat("  Flagged parameters:\n")
  for (i in seq_len(nrow(flagged_rows))) {
    r <- flagged_rows[i, ]
    cat(sprintf("    %-12s  %s  DER = %.4f\n",
                r$param_name, r$tier, r$der))
  }
}
cat("\n")


# =============================================================================
# Section 12: Selective Cholesky correction
# =============================================================================
# Algorithm 1, Step 3: For each flagged parameter k, the posterior draws
# are re-centered and re-scaled:
#   phi_corrected^(m) = phi_hat + sqrt(DER_k) * (phi^(m) - phi_hat)
#
# This preserves the posterior mean exactly while widening the marginal
# posterior variance to match the sandwich variance. Unflagged parameters
# pass through unchanged.

cat("--- 12. Selective correction ---\n")

point_est <- c(beta_hat, theta_hat)
stopifnot(length(point_est) == D_TOTAL)

## Strategy A: Naive (no correction)
draws_naive <- draws_all

## Strategy B: Selective (DER > tau only)
draws_selective <- draws_all
scale_factors_sel <- rep(1.0, D_TOTAL)

for (k in seq_len(D_TOTAL)) {
  if (classification$flagged[k]) {
    sf <- sqrt(diag_Vs[k] / diag_mcmc[k])  # = sqrt(DER_k)
    scale_factors_sel[k] <- sf
    draws_selective[, k] <- point_est[k] + sf * (draws_all[, k] - point_est[k])
  }
}

## Strategy C: Blanket (correct all)
draws_blanket <- draws_all
scale_factors_blk <- sqrt(diag_Vs / diag_mcmc)

for (k in seq_len(D_TOTAL)) {
  draws_blanket[, k] <- point_est[k] +
    scale_factors_blk[k] * (draws_all[, k] - point_est[k])
}

cat(sprintf("  Selective: %d parameter(s) corrected\n", n_flagged))
if (n_flagged > 0) {
  for (k in which(classification$flagged)) {
    cat(sprintf("    %s: sf = %.4f (sqrt(DER) = %.4f)\n",
                param_names[k], scale_factors_sel[k], sqrt(der_all[k])))
  }
}


# =============================================================================
# Section 13: Credible interval comparison
# =============================================================================

cat("\n--- 13. 90%% credible interval comparison ---\n")

## Compute 90% CIs for each strategy
compute_ci <- function(draws_mat, prob = 0.90) {
  alpha <- (1 - prob) / 2
  probs <- c(alpha, 1 - alpha)
  P <- ncol(draws_mat)
  data.frame(
    param_idx = seq_len(P),
    mean   = colMeans(draws_mat),
    lower  = apply(draws_mat, 2, quantile, probs = probs[1]),
    upper  = apply(draws_mat, 2, quantile, probs = probs[2]),
    width  = apply(draws_mat, 2, function(x) diff(quantile(x, probs)))
  )
}

ci_naive     <- compute_ci(draws_naive)
ci_selective <- compute_ci(draws_selective)
ci_blanket   <- compute_ci(draws_blanket)

## Build comparison table
comparison <- data.frame(
  param_name      = param_names,
  param_type      = param_types,
  der             = as.numeric(der_all),
  flagged         = classification$flagged,
  width_naive     = ci_naive$width,
  width_selective = ci_selective$width,
  width_blanket   = ci_blanket$width,
  ratio_sel_naive = ci_selective$width / ci_naive$width,
  ratio_blk_naive = ci_blanket$width / ci_naive$width,
  stringsAsFactors = FALSE
)

## Report for fixed effects
cat("  Fixed effects:\n")
cat(sprintf("    %-14s  DER       Naive-W   Select-W  Sel/Nai\n", "Parameter"))
cat(sprintf("    %s\n", paste(rep("-", 60), collapse = "")))
for (i in seq_len(P_FIXED)) {
  r <- comparison[i, ]
  cat(sprintf("    %-14s  %7.4f   %8.5f  %8.5f  %7.4f\n",
              covariate_names[i], r$der,
              r$width_naive, r$width_selective, r$ratio_sel_naive))
}

## Blanket over-correction analysis
n_narrowed <- sum(comparison$der < 1.0 & comparison$ratio_blk_naive < 0.999)
cat(sprintf("\n  Blanket correction: %d / %d parameters inappropriately narrowed\n",
            n_narrowed, D_TOTAL))
if (n_narrowed > 0) {
  min_ratio <- min(comparison$ratio_blk_naive[comparison$der < 1.0])
  cat(sprintf("  Worst blanket narrowing: ratio = %.4f (%.1f%% of original width)\n",
              min_ratio, 100 * min_ratio))
}


# =============================================================================
# Section 14: Save results
# =============================================================================

cat("\n--- 14. Saving DER results ---\n")

der_results <- list(

  ## DER values
  der             = der_all,           # named numeric [54]
  der_beta        = der_beta,          # named numeric [3]
  der_theta       = der_theta,         # named numeric [51]
  der_laplace     = der_laplace,       # named numeric [54]

  ## Sandwich ingredients
  H_obs           = H_obs,            # 54 x 54
  H_obs_inv       = H_obs_inv,        # 54 x 54
  J_cluster       = J_cluster,        # 54 x 54
  V_sand          = V_sand,           # 54 x 54

  ## MCMC covariance
  sigma_mcmc      = sigma_mcmc,       # 54 x 54

  ## Posterior means
  beta_hat        = beta_hat,         # named numeric [3]
  theta_hat       = theta_hat,        # numeric [51]
  sigma_theta_hat = sigma_theta_hat,  # scalar

  ## Per-state design quantities
  deff_j          = deff_j,           # numeric [51]
  B_j             = B_j,              # numeric [51]

  ## Classification
  classification  = classification,   # data.frame [54 rows]
  tau_threshold   = TAU_THRESHOLD,
  n_flagged       = n_flagged,

  ## Correction results
  comparison      = comparison,        # data.frame [54 rows]
  scale_factors_sel = scale_factors_sel,
  scale_factors_blk = scale_factors_blk,
  ci_naive        = ci_naive,
  ci_selective    = ci_selective,
  ci_blanket      = ci_blanket,

  ## Corrected draws
  draws_naive     = draws_naive,       # [M x 54]
  draws_selective = draws_selective,    # [M x 54]
  draws_blanket   = draws_blanket,     # [M x 54]
  point_est       = point_est,         # numeric [54]

  ## Diagonal vectors (for downstream analysis)
  diag_V_sand     = diag_Vs,
  diag_MCMC       = diag_mcmc,
  diag_Hinv       = diag_H_inv,

  ## Metadata
  param_names     = param_names,
  param_types     = param_types,
  covariate_names = covariate_names,
  n_obs           = N_OBS,
  n_states        = J_STATES,
  n_psu           = G_psu,
  n_fe            = P_FIXED,
  beta_prior_sd   = BETA_PRIOR_SD,
  n_draws         = M,
  computed_at     = Sys.time()
)

der_results_path <- file.path(PRECOMP_DIR, "der_results.rds")
saveRDS(der_results, der_results_path)
cat(sprintf("  Saved: %s (%.2f MB)\n",
            der_results_path, file.info(der_results_path)$size / 1024^2))


# =============================================================================
# Verification checks
# =============================================================================

cat("\n--- Verification ---\n")

## Check 1: Only beta[2] (poverty_cwc) flagged at tau = 1.2
flagged_names <- classification$param_name[classification$flagged]
cat(sprintf("  Flagged parameters: %s\n", paste(flagged_names, collapse = ", ")))

## Check 2: RE DER all < 1
n_re_gt1 <- sum(der_theta >= 1.0)
cat(sprintf("  RE with DER >= 1: %d [%s]\n",
            n_re_gt1, ifelse(n_re_gt1 == 0, "PASS", "NOTE")))

## Check 3: Selective correction preserves means
mean_shift <- max(abs(colMeans(draws_selective) - colMeans(draws_naive)))
cat(sprintf("  Max mean shift (selective): %.6f [%s]\n",
            mean_shift, ifelse(mean_shift < 0.01, "PASS", "WARNING")))

## Check 4: Unflagged parameters unchanged under selective
unflagged_max_dev <- max(abs(comparison$ratio_sel_naive[!classification$flagged] - 1.0))
cat(sprintf("  Unflagged width ratio deviation: %.2e [%s]\n",
            unflagged_max_dev, ifelse(unflagged_max_dev < 1e-12, "PASS", "WARNING")))


# =============================================================================
# Summary
# =============================================================================

cat("\n==============================================================\n")
cat("  DER ANALYSIS SUMMARY\n")
cat("==============================================================\n")
cat(sprintf("  Parameters: %d total (%d FE + %d RE)\n",
            D_TOTAL, P_FIXED, J_STATES))
cat(sprintf("  Threshold: tau = %.1f\n", TAU_THRESHOLD))
cat(sprintf("  Flagged: %d / %d (%.1f%%)\n",
            n_flagged, D_TOTAL, 100 * n_flagged / D_TOTAL))
cat(sprintf("\n  Key DER values:\n"))
cat(sprintf("    poverty_cwc (within):    DER = %.4f -> CORRECT\n", der_beta[2]))
cat(sprintf("    intercept   (between):   DER = %.4f -> retain\n", der_beta[1]))
cat(sprintf("    tiered_reim (between):   DER = %.4f -> retain\n", der_beta[3]))
cat(sprintf("    RE mean:                 DER = %.4f -> all retain\n",
            mean(der_theta)))
cat(sprintf("\n  Selective correction:\n"))
cat(sprintf("    poverty_cwc width change: +%.1f%%\n",
            (comparison$ratio_sel_naive[2] - 1) * 100))
cat(sprintf("  Blanket correction:\n"))
cat(sprintf("    Parameters inappropriately narrowed: %d / %d\n",
            n_narrowed, D_TOTAL))
cat("==============================================================\n")
cat("  Step 3 Complete (Track A).\n")
cat("==============================================================\n")

}  # end Track A
