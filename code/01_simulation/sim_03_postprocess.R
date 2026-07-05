# =============================================================================
# sim_03_postprocess.R: DER Computation, Cholesky Correction, and Coverage
# =============================================================================
#
# Purpose : Compute the Design Effect Ratio (DER) for each parameter from a
#           fitted Stan model, apply Cholesky-based posterior corrections
#           (blanket and selective), and evaluate frequentist coverage of
#           the resulting credible intervals.
#
#           This script implements the core methodological contribution of
#           the paper: the sandwich-based DER diagnostic and the selective
#           Cholesky correction that targets only survey-sensitive parameters.
#
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Sections 3.1-3.4 (DER Framework), Section 4.3 (Correction Strategies)
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Track   : A (Full Replication)
# Inputs  : sim_02_fit.R (for CmdStanMCMC fit objects)
#           sim_01_dgp.R (for data_list with true parameter values)
# Outputs : Functions for computing DER, applying corrections, evaluating
#           coverage. Called per-replication from the single-rep pipeline.
#
# Mathematical outline:
#
#   The DER for parameter k is defined as (Section 3.1, Eq. 9):
#
#     DER_k = [V_sand]_kk / [Sigma_MCMC]_kk
#
#   where:
#     V_sand = H_obs^{-1} J_cluster H_obs^{-1}   (sandwich variance)
#     Sigma_MCMC = Cov(phi | y)                   (MCMC posterior covariance)
#
#   The sandwich components are:
#     H_obs     : Observed information matrix (negative Hessian of the
#                 weighted log-posterior), d x d where d = p + J.
#     J_cluster : Clustered score outer-product matrix, capturing
#                 within-cluster correlation of score contributions.
#
#   DER > 1 indicates the parameter is survey-sensitive: the model-based
#   posterior underestimates uncertainty relative to the design-based
#   sandwich estimator. DER ~ 1 indicates adequate model-based inference.
#
#   Correction strategies (Section 4.3):
#     - Naive: no correction (baseline).
#     - Blanket Cholesky: correct ALL parameters using V_sand (Williams &
#       Savitsky, 2021). May over-correct parameters with DER ~ 1.
#     - Selective DER-tau: correct only parameters with DER_k > tau.
#       Preserves valid model-based inference for insensitive parameters.
# =============================================================================

library(Matrix)


# =============================================================================
# Section 0: Shared Helper Functions (single canonical definitions)
# =============================================================================
# The sandwich building blocks (build_H_obs_logistic, build_J_cluster,
# build_J_cluster_psu, compute_mcmc_covariance, compute_group_deff,
# compute_shrinkage_factors) and the DER utilities (classify_parameters,
# selective_correct, compute_theoretical_der, compute_der) are defined ONCE
# in code/helpers/ and sourced here. Keeping a single definition of each name
# avoids the risk of divergent local copies; do not re-introduce them here.
# =============================================================================

local({
  helper_files <- c("sandwich_functions.R", "der_functions.R")
  candidate_dirs <- c(
    if (requireNamespace("here", quietly = TRUE)) {
      file.path(here::here(), "code", "helpers")
    },
    file.path(getwd(), "code", "helpers"),
    file.path(getwd(), "..", "helpers")
  )
  helper_dir <- NULL
  for (dir in candidate_dirs) {
    if (all(file.exists(file.path(dir, helper_files)))) {
      helper_dir <- dir
      break
    }
  }
  if (is.null(helper_dir)) {
    stop("Cannot locate code/helpers/. Source sim_03_postprocess.R from ",
         "the replication package root (or open the .Rproj).")
  }
  for (f in helper_files) {
    source(file.path(helper_dir, f))
  }
})


# =============================================================================
# Section 1: Main Sandwich/DER Computation
# =============================================================================

#' Compute sandwich variance and DER from a fitted Stan model
#'
#' Constructs H_obs, J_cluster, and V_sand for the hierarchical logistic
#' model, then computes parameter-specific DER values by comparing V_sand
#' to the MCMC posterior covariance.
#'
#' @param fit A CmdStanMCMC fit object.
#' @param data_list A list with: y (binary outcomes), X (design matrix),
#'   group (cluster indicators), w (survey weights).
#' @param beta_prior_sd Numeric scalar. Prior SD for fixed effects
#'   (default 5, corresponding to beta ~ N(0, 25)).
#' @param strata Optional integer vector (N). Stratum indicator for the
#'   meat; NULL (default) = single stratum. In the balanced simulation the
#'   design has no strata, so NULL is the correct choice.
#' @param center Logical. Center cluster score totals within strata in the
#'   meat (default TRUE; the legacy alternative FALSE reproduces the earlier
#'   estimator).
#' @param df_correct Logical. Apply the C_h/(C_h - 1) finite-cluster
#'   correction in the meat (default TRUE; the legacy alternative FALSE
#'   reproduces the earlier estimator).
#'
#' @return A list with: H_obs, H_obs_inv, J_cluster, V_sand, der,
#'   der_beta, der_theta, sigma_mcmc, beta_hat, theta_hat,
#'   sigma_theta_hat, deff_j, B_j, meat_options.
compute_sandwich_hlr <- function(fit, data_list, beta_prior_sd = 5,
                                 strata = NULL, center = TRUE,
                                 df_correct = TRUE) {

  # -- 1a. Validate inputs ----------------------------------------------------
  stopifnot(
    is.list(data_list),
    all(c("y", "X", "group", "w") %in% names(data_list))
  )

  y     <- as.integer(data_list$y)
  X     <- as.matrix(data_list$X)
  group <- as.integer(data_list$group)
  w     <- as.numeric(data_list$w)

  N <- length(y)
  p <- ncol(X)
  J <- max(group)
  d <- p + J   # total parameter dimension: phi = (beta, theta)

  stopifnot(
    nrow(X) == N, length(group) == N, length(w) == N,
    all(group >= 1L & group <= J),
    all(y %in% c(0L, 1L)),
    all(w > 0)
  )

  # -- 1b. Extract posterior means from Stan fit ------------------------------
  fit_summary <- fit$summary()

  # Fixed effects: beta[1], ..., beta[p]
  beta_hat <- numeric(p)
  for (k in seq_len(p)) {
    idx <- which(fit_summary$variable == paste0("beta[", k, "]"))
    if (length(idx) == 0) {
      stop("Parameter beta[", k, "] not found in fit summary.")
    }
    beta_hat[k] <- fit_summary$mean[idx]
  }

  # Random effects: theta[1], ..., theta[J]
  theta_hat <- numeric(J)
  for (j in seq_len(J)) {
    idx <- which(fit_summary$variable == paste0("theta[", j, "]"))
    if (length(idx) == 0) {
      stop("Parameter theta[", j, "] not found in fit summary.")
    }
    theta_hat[j] <- fit_summary$mean[idx]
  }

  # Random-effect SD
  sigma_idx <- which(fit_summary$variable == "sigma_theta")
  if (length(sigma_idx) == 0) {
    stop("Parameter 'sigma_theta' not found in fit summary.")
  }
  sigma_theta_hat <- fit_summary$mean[sigma_idx]
  stopifnot(sigma_theta_hat > 0)

  # -- 1c. Compute predicted probabilities ------------------------------------
  # Linear predictor: eta_i = X_i' beta_hat + theta_hat[group_i]
  eta <- as.numeric(X %*% beta_hat) + theta_hat[group]
  q   <- plogis(eta)           # predicted probability
  wt  <- q * (1 - q)           # Bernoulli variance function (working weights)

  # -- 1d. Build H_obs (observed information matrix, Section 3.1) -------------
  #
  # H_obs = H_data + H_prior, where:
  #   H_data[a,b] = sum_i w_i * q_i * (1 - q_i) * phi_ia * phi_ib
  # Here phi_i is the d-vector of derivatives of the linear predictor
  # w.r.t. phi = (beta, theta):
  #   d(eta_i)/d(beta_k) = X[i,k]
  #   d(eta_i)/d(theta_j) = 1 if group[i] == j, else 0
  #
  # The block structure is:
  #   H_obs = | H_bb   H_bt |   +   | diag(1/sd_beta^2)  0                |
  #           | H_bt'  H_tt |       | 0                   diag(1/sigma^2)  |
  #
  # where H_bb = X' diag(v) X is the p x p beta-beta block,
  # H_bt is the p x J beta-theta coupling block, and
  # H_tt is a J x J diagonal matrix (theta-theta block).

  H_obs <- build_H_obs_logistic(X, group, w, wt, p, J, N,
                                sigma_theta_hat, beta_prior_sd)

  # -- 1e. Build J_cluster (clustered score outer-product, Section 3.1) -------
  #
  # Default (center = TRUE, df_correct = TRUE; single stratum):
  #   J_cluster = [J/(J-1)] * sum_j (s_j - s_bar)(s_j - s_bar)'
  # Legacy alternative (center = FALSE, df_correct = FALSE):
  #   J_cluster = sum_j s_j s_j'
  # where s_j is the d-vector of cluster-level score sums:
  #   s_j[k] = sum_{i in j} w_i * (y_i - q_i) * X[i,k]   for k = 1..p
  #   s_j[p+j] = sum_{i in j} w_i * (y_i - q_i)           for the theta_j entry
  #   s_j[p+l] = 0 for l != j
  #
  # This matrix captures the within-cluster correlation of score contributions
  # that arises from the complex survey design. See
  # code/helpers/sandwich_functions.R for the canonical implementation.

  r <- w * (y - q)   # weighted score residuals
  J_cluster <- build_J_cluster(X, group, r, p, J, N,
                               cluster = group, strata = strata,
                               center = center, df_correct = df_correct)

  # -- 1f. Assemble sandwich and compute DER ----------------------------------
  #
  # V_sand = H_obs^{-1} J_cluster H_obs^{-1}   (Eq. 6)
  #
  # This is the survey-design-corrected variance estimator. Under complex
  # sampling, V_sand is consistent even when the model does not account
  # for the sampling design (Binder, 1983; Rao & Molina, 2015).

  H_obs_inv <- tryCatch(
    solve(H_obs),
    error = function(e) {
      warning("H_obs is singular or near-singular. Using nearPD fallback. ",
              "Condition number: ", kappa(H_obs, exact = TRUE))
      H_pd <- as.matrix(Matrix::nearPD(H_obs, keepDiag = TRUE)$mat)
      solve(H_pd)
    }
  )

  V_sand <- H_obs_inv %*% J_cluster %*% H_obs_inv
  V_sand <- (V_sand + t(V_sand)) / 2   # enforce symmetry

  # -- 1g. Compute MCMC posterior covariance ----------------------------------
  #
  # Sigma_MCMC is the sample covariance of posterior draws from MCMC.
  # This serves as the denominator in the DER definition (Eq. 9):
  #   DER_k = V_sand[k,k] / Sigma_MCMC[k,k]
  #
  # Using Sigma_MCMC (rather than the Laplace approximation H_obs^{-1})
  # as the denominator ensures that the DER correctly accounts for
  # posterior non-Gaussianity in finite samples.

  sigma_mcmc <- compute_mcmc_covariance(fit, p, J)

  # -- 1h. Compute DER (Eq. 9) -----------------------------------------------
  #
  # DER_k = V_sand[k,k] / Sigma_MCMC[k,k]
  #
  # Interpretation (Section 3.2):
  #   DER_k ~ 1 : parameter k is insensitive to survey design
  #   DER_k > 1 : parameter k is survey-sensitive (under-coverage risk)
  #   DER_k < 1 : parameter k is over-covered by the model
  #
  # The three-tier classification (Section 3.3):
  #   Tier I-a: Within-group FE, DER_k ~ DEFF
  #   Tier I-b: Between-group FE, DER_k ~ DEFF * (1 - B)
  #   Tier II:  Random effects, DER_j = B_j * DEFF_j * kappa_j(J)

  diag_V    <- diag(V_sand)
  diag_mcmc <- diag(sigma_mcmc)

  stopifnot(all(diag_mcmc > 0))
  der <- diag_V / diag_mcmc

  # Name the DER vector
  der_names <- c(paste0("beta[", seq_len(p), "]"),
                 paste0("theta[", seq_len(J), "]"))
  names(der) <- der_names

  der_beta  <- der[seq_len(p)]
  der_theta <- der[(p + 1):d]

  # Laplace-based DER (secondary diagnostic: sandwich vs model
  # information; measures posterior Gaussianity when compared to der)
  der_laplace <- diag_V / diag(H_obs_inv)
  names(der_laplace) <- der_names

  # -- 1i. Per-cluster design effects and shrinkage factors -------------------
  deff_j <- compute_group_deff(group, w, J)
  B_j    <- compute_shrinkage_factors(group, w, wt, sigma_theta_hat, J)

  # -- 1j. Return all results -------------------------------------------------
  list(
    H_obs           = H_obs,
    H_obs_inv       = H_obs_inv,
    J_cluster       = J_cluster,
    V_sand          = V_sand,
    der             = der,
    der_laplace     = der_laplace,
    der_beta        = der_beta,
    der_theta       = der_theta,
    sigma_mcmc      = sigma_mcmc,
    beta_hat        = beta_hat,
    theta_hat       = theta_hat,
    sigma_theta_hat = sigma_theta_hat,
    deff_j          = deff_j,
    B_j             = B_j,
    meat_options    = list(cluster    = "group",
                           stratified = !is.null(strata),
                           center     = center,
                           df_correct = df_correct)
  )
}


# =============================================================================
# Sections 2-5: H_obs / J_cluster / Sigma_MCMC / DEFF / B_j (sourced helpers)
# =============================================================================
# The functions for these quantities (build_H_obs_logistic,
# build_J_cluster, compute_mcmc_covariance, compute_group_deff,
# compute_shrinkage_factors) live in code/helpers/sandwich_functions.R and
# are sourced from Section 0 rather than redefined here. Keeping two copies
# of the same names would be an error hazard. The helpers' build_J_cluster()
# additionally supports the stratified / centered / DF-corrected meat.
# =============================================================================


# =============================================================================
# Section 6: Cholesky Correction Strategies
# =============================================================================
#
# The Cholesky affine transformation (Williams & Savitsky, 2021) replaces
# the MCMC covariance with the sandwich covariance while preserving the
# posterior mean:
#
#   phi_corrected[m] = phi_hat + L_sand * L_mcmc^{-1} * (phi[m] - phi_hat)
#
# where L_sand = chol(V_sand) and L_mcmc = chol(Sigma_MCMC).
#
# Three strategies are implemented:
#
# 1. Naive: no correction. Baseline for comparison.
# 2. Blanket: correct ALL d parameters. This is the standard approach
#    but may over-inflate intervals for insensitive parameters.
# 3. Selective DER-tau: correct only parameters with DER_k > tau.
#    Parameters with DER_k <= tau retain their original MCMC intervals.
#    This is the paper's main innovation (Section 3.4).
# =============================================================================

#' Apply selective Cholesky correction to flagged parameters
#'
#' @param mcmc_draws M x d matrix of MCMC draws.
#' @param sandwich_result List from compute_sandwich_hlr().
#' @param threshold Numeric. DER threshold tau for flagging.
#' @return List with: corrected_draws, flagged_indices, n_corrected.
apply_selective_cholesky <- function(mcmc_draws, sandwich_result,
                                    threshold = 1.2) {

  stopifnot(is.matrix(mcmc_draws), is.list(sandwich_result))

  M <- nrow(mcmc_draws)
  d <- ncol(mcmc_draws)

  der <- sandwich_result$der
  stopifnot(length(der) == d)

  # -- Identify flagged parameters (DER > tau) --------------------------------
  flagged   <- which(der > threshold)
  n_flagged <- length(flagged)

  # If nothing flagged, return original draws unchanged
  if (n_flagged == 0) {
    return(list(
      corrected_draws    = mcmc_draws,
      flagged_indices    = integer(0),
      n_corrected        = 0L,
      fraction_corrected = 0,
      threshold          = threshold,
      der_flagged        = numeric(0)
    ))
  }

  # -- Extract submatrices for flagged parameters -----------------------------
  V_sand_S <- sandwich_result$V_sand[flagged, flagged, drop = FALSE]
  Sigma_S  <- sandwich_result$sigma_mcmc[flagged, flagged, drop = FALSE]

  phi_hat   <- c(sandwich_result$beta_hat, sandwich_result$theta_hat)
  phi_hat_S <- phi_hat[flagged]

  # -- Ensure positive definiteness -------------------------------------------
  V_sand_S <- ensure_pd(V_sand_S, label = "V_sand[S,S]")
  Sigma_S  <- ensure_pd(Sigma_S, label = "Sigma_MCMC[S,S]")

  # -- Compute Cholesky factors -----------------------------------------------
  # R's chol() returns upper triangular; we need lower triangular.
  L_sand <- tryCatch(t(chol(V_sand_S)), error = function(e) {
    warning("Cholesky of V_sand[S,S] failed: ", e$message)
    V_pd <- as.matrix(Matrix::nearPD(V_sand_S, keepDiag = TRUE)$mat)
    t(chol(V_pd))
  })

  L_mcmc <- tryCatch(t(chol(Sigma_S)), error = function(e) {
    warning("Cholesky of Sigma[S,S] failed: ", e$message)
    S_pd <- as.matrix(Matrix::nearPD(Sigma_S, keepDiag = TRUE)$mat)
    t(chol(S_pd))
  })

  # -- Affine transformation: A = L_sand * L_mcmc^{-1} -----------------------
  A <- L_sand %*% solve(L_mcmc)

  # -- Apply transformation to draws ------------------------------------------
  corrected_draws <- mcmc_draws

  # Center, transform, shift back
  centered_S <- sweep(mcmc_draws[, flagged, drop = FALSE], 2, phi_hat_S)
  transformed_S <- t(A %*% t(centered_S))
  corrected_draws[, flagged] <- sweep(transformed_S, 2, -phi_hat_S)

  list(
    corrected_draws    = corrected_draws,
    flagged_indices    = flagged,
    n_corrected        = n_flagged,
    fraction_corrected = n_flagged / d,
    threshold          = threshold,
    der_flagged        = der[flagged]
  )
}


#' Apply blanket Cholesky correction to ALL parameters
#'
#' @param mcmc_draws M x d matrix of MCMC draws.
#' @param sandwich_result List from compute_sandwich_hlr().
#' @return List with corrected_draws and metadata.
apply_blanket_cholesky <- function(mcmc_draws, sandwich_result) {
  # Blanket = selective with threshold = -Inf (flag everything)
  result <- apply_selective_cholesky(mcmc_draws, sandwich_result,
                                    threshold = -Inf)
  result$threshold <- NA_real_
  result
}


# =============================================================================
# Section 7: Credible Interval and Coverage Computation
# =============================================================================

#' Compute quantile-based credible intervals
#'
#' @param draws M x d matrix of posterior draws.
#' @param level Numeric. Credible level (default 0.90 for 90% CI).
#' @return Data.frame with param, lower, upper, width, median, mean.
compute_credible_intervals <- function(draws, level = 0.90) {
  stopifnot(is.matrix(draws), level > 0, level < 1)

  alpha <- (1 - level) / 2
  probs <- c(alpha, 0.5, 1 - alpha)

  d <- ncol(draws)
  param_names <- colnames(draws)
  if (is.null(param_names)) {
    param_names <- paste0("param[", seq_len(d), "]")
  }

  quants <- apply(draws, 2, quantile, probs = probs)
  means  <- colMeans(draws)

  data.frame(
    param  = param_names,
    lower  = quants[1, ],
    median = quants[2, ],
    upper  = quants[3, ],
    width  = quants[3, ] - quants[1, ],
    mean   = means,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


#' Compute coverage indicators for one replication
#'
#' Checks whether the true parameter values fall within the credible
#' intervals. Returns a logical vector.
#'
#' @param true_values Numeric vector of length d. True parameter values.
#' @param ci_df Data.frame from compute_credible_intervals().
#' @return Named logical vector of length d.
compute_coverage <- function(true_values, ci_df) {
  stopifnot(length(true_values) == nrow(ci_df))

  covered <- (true_values >= ci_df$lower) & (true_values <= ci_df$upper)
  names(covered) <- ci_df$param
  covered
}


# =============================================================================
# Section 8: Utility: Ensure Positive Definiteness
# =============================================================================

#' Ensure a matrix is positive definite
#'
#' Checks eigenvalues and applies Matrix::nearPD if needed. This is a
#' defensive measure against numerical issues in the sandwich computation.
#'
#' @param mat Symmetric numeric matrix.
#' @param label Character. Label for warning messages.
#' @param tol Numeric. Minimum acceptable eigenvalue ratio (default 1e-10).
#' @return A positive definite matrix.
ensure_pd <- function(mat, label = "matrix", tol = 1e-10) {
  stopifnot(is.matrix(mat), nrow(mat) == ncol(mat))
  mat <- (mat + t(mat)) / 2

  eig <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig)
  max_eig <- max(eig)

  if (min_eig <= 0 || (max_eig > 0 && min_eig / max_eig < tol)) {
    warning(label, " is not sufficiently positive definite. ",
            "Min eigenvalue = ", format(min_eig, digits = 4),
            ". Applying nearPD correction.")
    mat <- as.matrix(Matrix::nearPD(mat, keepDiag = TRUE)$mat)
  }

  mat
}


# =============================================================================
# Section 9: DER Classification and Theoretical Predictions
# =============================================================================
# classify_parameters() and compute_theoretical_der() are sourced from
# code/helpers/der_functions.R (see Section 0). The classification follows
# the paper's type-based taxonomy (Tier I-a within-FE, I-b between-FE,
# II random effects, III reserved for hyperparameters). This is the single
# taxonomy used across the package; no alternative tier scheme is defined
# in this file.
# =============================================================================

#' Classify parameters by type: within-group FE, between-group FE, RE
#'
#' In our simulation design:
#'   - beta[1] (intercept): between-group FE
#'   - beta[2] (x_within):  within-group FE (group-mean centered)
#'   - beta[3] (z_between): between-group FE
#'   - theta[1..J]:         random effects
#'
#' @param p Integer. Number of fixed effects.
#' @param J Integer. Number of clusters.
#'
#' @return Character vector of length d = p + J with parameter types.
classify_param_types <- function(p, J) {
  d <- p + J
  types <- character(d)

  types[1] <- "fe_between"    # intercept
  types[2] <- "fe_within"     # x_within (group-mean centered)
  types[3] <- "fe_between"    # z_between
  types[(p + 1):d] <- "re"

  types
}


# =============================================================================
# Section 10: Current Runner Location
# =============================================================================
#
# The active v2 single-replication pipeline is defined in sim_04_run_single.R
# and driven by sim_05_run_batch.R. This file provides post-processing helpers
# used by that runner.


# =============================================================================
# Section 11: Self-Validation
# =============================================================================

if (interactive()) {
  cat("=== sim_03_postprocess.R: Self-validation ===\n\n")

  set.seed(123)
  N_test <- 100
  J_test <- 5
  p_test <- 3

  X_test <- cbind(1, rnorm(N_test),
                  rep(rnorm(J_test), each = N_test / J_test))
  group_test <- rep(1:J_test, each = N_test / J_test)
  w_test     <- runif(N_test, 0.5, 2.0)
  wt_test    <- runif(N_test, 0.1, 0.25)

  # Test H_obs construction
  H_test <- build_H_obs_logistic(X_test, group_test, w_test, wt_test,
                                 p_test, J_test, N_test,
                                 sigma_theta = 0.5, beta_prior_sd = 5)
  d_test <- p_test + J_test
  stopifnot(nrow(H_test) == d_test, ncol(H_test) == d_test)
  stopifnot(max(abs(H_test - t(H_test))) < 1e-12)
  cat("  H_obs construction: OK\n")

  # Test J_cluster construction
  r_test <- rnorm(N_test)
  Jmat_test <- build_J_cluster(X_test, group_test, r_test,
                               p_test, J_test, N_test)
  stopifnot(nrow(Jmat_test) == d_test, ncol(Jmat_test) == d_test)
  eig_vals <- eigen(Jmat_test, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(all(eig_vals >= -1e-10))
  cat("  J_cluster construction: OK (PSD verified)\n")

  # Test group DEFF
  deff_test <- compute_group_deff(group_test, w_test, J_test)
  stopifnot(length(deff_test) == J_test, all(deff_test >= 1))
  cat("  Group DEFF: OK\n")

  # Test shrinkage factors
  B_test <- compute_shrinkage_factors(group_test, w_test, wt_test,
                                      sigma_theta = 0.5, J_test)
  stopifnot(all(B_test > 0), all(B_test < 1))
  cat("  Shrinkage factors: OK\n")

  cat("\nAll self-validation checks passed.\n")
}
