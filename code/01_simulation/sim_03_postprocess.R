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
#           coverage. Called per-replication from sim_04_run.R.
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
#'
#' @return A list with: H_obs, H_obs_inv, J_cluster, V_sand, der,
#'   der_beta, der_theta, sigma_mcmc, beta_hat, theta_hat,
#'   sigma_theta_hat, deff_j, B_j.
compute_sandwich_hlr <- function(fit, data_list, beta_prior_sd = 5) {

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
  # J_cluster = sum_j s_j s_j'
  # where s_j is the d-vector of cluster-level score sums:
  #   s_j[k] = sum_{i in j} w_i * (y_i - q_i) * X[i,k]   for k = 1..p
  #   s_j[p+j] = sum_{i in j} w_i * (y_i - q_i)           for the theta_j entry
  #   s_j[p+l] = 0 for l != j
  #
  # This matrix captures the within-cluster correlation of score contributions
  # that arises from the complex survey design.

  r <- w * (y - q)   # weighted score residuals
  J_cluster <- build_J_cluster(X, group, r, p, J, N)

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
    der_beta        = der_beta,
    der_theta       = der_theta,
    sigma_mcmc      = sigma_mcmc,
    beta_hat        = beta_hat,
    theta_hat       = theta_hat,
    sigma_theta_hat = sigma_theta_hat,
    deff_j          = deff_j,
    B_j             = B_j
  )
}


# =============================================================================
# Section 2: H_obs Construction (Block-Efficient)
# =============================================================================
#
# H_obs is the observed information matrix (negative Hessian) of the
# weighted log-posterior. For computational efficiency, we construct it
# in block form rather than forming the full N x d design matrix.
#
# Block structure for phi = (beta_1..beta_p, theta_1..theta_J):
#
#   H_data = | X' diag(v) X       X' diag(v) Z  |
#            | Z' diag(v) X       diag(h_jj)    |
#
# where v_i = w_i * q_i * (1 - q_i) and Z is the N x J cluster indicator.
# Since Z is sparse (each row has exactly one 1), we compute the blocks
# using group-level summations.
# =============================================================================

#' Build the observed information matrix for hierarchical logistic regression
#'
#' @param X N x p design matrix.
#' @param group N-vector of cluster indices (1..J).
#' @param w N-vector of survey weights.
#' @param wt N-vector of working weights q*(1-q).
#' @param p Number of fixed effects.
#' @param J Number of clusters.
#' @param N Number of observations.
#' @param sigma_theta Posterior mean of sigma_theta.
#' @param beta_prior_sd Prior SD for beta (default 5).
#'
#' @return d x d observed information matrix (d = p + J).
build_H_obs_logistic <- function(X, group, w, wt, p, J, N,
                                 sigma_theta, beta_prior_sd = 5) {
  d <- p + J
  H <- matrix(0, d, d)

  # Combined weight: survey weight * working weight
  v <- w * wt   # length N

  # -- Beta-beta block (p x p): H_bb = X' diag(v) X --------------------------
  H_bb <- crossprod(X, X * v)
  H[1:p, 1:p] <- H_bb

  # -- Beta-theta block (p x J) and theta-theta diagonal ----------------------
  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    v_j   <- v[idx_j]

    # Beta-theta coupling: sum_i v_i * X[i,] for i in cluster j
    if (length(idx_j) == 1) {
      bt_j <- v_j * X[idx_j, ]
    } else {
      bt_j <- colSums(X[idx_j, , drop = FALSE] * v_j)
    }
    H[1:p, p + j] <- bt_j
    H[p + j, 1:p] <- bt_j

    # Theta-theta diagonal: sum_i v_i for i in cluster j
    H[p + j, p + j] <- sum(v_j)
  }

  # -- Add prior contributions ------------------------------------------------

  # Beta prior: beta ~ N(0, beta_prior_sd^2)
  # => precision contribution = 1/beta_prior_sd^2 added to diagonal
  if (is.finite(beta_prior_sd)) {
    beta_prior_prec <- 1 / beta_prior_sd^2
    for (k in seq_len(p)) {
      H[k, k] <- H[k, k] + beta_prior_prec
    }
  }

  # Theta prior: theta_j ~ N(0, sigma_theta^2)
  # => precision contribution = 1/sigma_theta^2 added to diagonal
  theta_prior_prec <- 1 / sigma_theta^2
  for (j in seq_len(J)) {
    H[p + j, p + j] <- H[p + j, p + j] + theta_prior_prec
  }

  # Verify symmetry (should hold by construction)
  stopifnot(max(abs(H - t(H))) < 1e-10)

  H
}


# =============================================================================
# Section 3: J_cluster Construction (Cluster-Level Score Sums)
# =============================================================================
#
# J_cluster = sum_j s_j s_j' is a positive semi-definite matrix that
# captures within-cluster correlation of the score function. It is the
# "meat" of the sandwich variance estimator.
#
# Each cluster score vector s_j has dimension d = p + J, with:
#   - Entries 1..p: sum of weighted score residuals times covariates
#   - Entry p+j: sum of weighted score residuals (for theta_j)
#   - All other entries: zero (no cross-cluster contributions)
# =============================================================================

#' Build the clustered score outer-product matrix
#'
#' @param X N x p design matrix.
#' @param group N-vector of cluster indices.
#' @param r N-vector of weighted score residuals: r_i = w_i * (y_i - q_i).
#' @param p Number of fixed effects.
#' @param J Number of clusters.
#' @param N Number of observations.
#'
#' @return d x d clustered score outer-product matrix (PSD).
build_J_cluster <- function(X, group, r, p, J, N) {
  d    <- p + J
  Jmat <- matrix(0, d, d)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)

    # Construct the d-dimensional score sum for cluster j
    s_j <- numeric(d)

    # Fixed-effect score components: sum_{i in j} r_i * X[i,]
    if (length(idx_j) == 1) {
      s_j[1:p] <- r[idx_j] * X[idx_j, ]
    } else {
      s_j[1:p] <- colSums(X[idx_j, , drop = FALSE] * r[idx_j])
    }

    # Random-effect score component: sum_{i in j} r_i
    s_j[p + j] <- sum(r[idx_j])

    # Accumulate outer product: J_cluster += s_j s_j'
    Jmat <- Jmat + tcrossprod(s_j)
  }

  # Symmetrize (enforces numerical symmetry)
  Jmat <- (Jmat + t(Jmat)) / 2

  Jmat
}


# =============================================================================
# Section 4: MCMC Posterior Covariance
# =============================================================================

#' Compute the posterior covariance matrix from MCMC draws
#'
#' Extracts draws for beta and theta, merges chains, and computes the
#' sample covariance. This matrix forms the denominator of the DER.
#'
#' @param fit CmdStanMCMC object.
#' @param p Number of fixed effects.
#' @param J Number of clusters.
#'
#' @return d x d posterior covariance matrix (d = p + J).
compute_mcmc_covariance <- function(fit, p, J) {
  d <- p + J

  beta_names  <- paste0("beta[", seq_len(p), "]")
  theta_names <- paste0("theta[", seq_len(J), "]")
  all_names   <- c(beta_names, theta_names)

  # Extract all post-warmup draws, merged across chains
  draws <- fit$draws(variables = all_names, format = "matrix")
  stopifnot(is.matrix(draws), ncol(draws) == d)

  sigma_mcmc <- cov(draws)
  sigma_mcmc <- (sigma_mcmc + t(sigma_mcmc)) / 2   # enforce symmetry

  sigma_mcmc
}


# =============================================================================
# Section 5: Per-Cluster Design Effects and Shrinkage Factors
# =============================================================================
#
# These quantities appear in the theoretical DER decomposition (Theorem 2):
#   DER_theta_j = B_j * DEFF_j * kappa_j(J)
#
# DEFF_j is the Kish design effect for cluster j (Eq. 3):
#   DEFF_j = n_j * sum(w_i^2) / (sum w_i)^2     for i in cluster j
#
# B_j is the shrinkage factor (Eq. 10):
#   B_j = sigma^2_theta / (sigma^2_theta + V_tilde_j)
# where V_tilde_j = 1 / sum_{i in j} w_i * q_i * (1 - q_i)
# is the effective conditional variance from cluster j data.
#
# B_j controls how much the random effect theta_j is pulled toward the
# grand mean. Parameters with more shrinkage (B_j closer to 1) are more
# insensitive to the survey design.
# =============================================================================

#' Compute per-cluster Kish design effect
#'
#' @param group N-vector of cluster indices.
#' @param w N-vector of survey weights.
#' @param J Number of clusters.
#' @return Numeric vector of length J with DEFF_j values.
compute_group_deff <- function(group, w, J) {
  deff_j <- numeric(J)
  for (j in seq_len(J)) {
    idx_j  <- which(group == j)
    w_j    <- w[idx_j]
    n_j    <- length(idx_j)
    deff_j[j] <- n_j * sum(w_j^2) / (sum(w_j))^2
  }
  deff_j
}


#' Compute per-cluster shrinkage factors B_j
#'
#' @param group N-vector of cluster indices.
#' @param w N-vector of survey weights.
#' @param wt N-vector of working weights q*(1-q).
#' @param sigma_theta Posterior mean of sigma_theta.
#' @param J Number of clusters.
#' @return Numeric vector of length J with B_j in (0, 1).
compute_shrinkage_factors <- function(group, w, wt, sigma_theta, J) {
  sigma2 <- sigma_theta^2
  B_j    <- numeric(J)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    # Effective information from cluster j: sum_{i in j} w_i * q_i * (1-q_i)
    info_j    <- sum(w[idx_j] * wt[idx_j])
    V_tilde_j <- 1 / info_j
    B_j[j]    <- sigma2 / (sigma2 + V_tilde_j)
  }

  stopifnot(all(B_j > 0), all(B_j < 1))
  B_j
}


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
# The three-tier classification (Section 3.3, Table 2) assigns each
# parameter to one of three tiers based on its DER value relative to
# the threshold tau:
#
#   Tier I (survey-sensitive):    DER > tau
#     -> Apply Cholesky correction to inflate posterior variance
#   Tier II (survey-neutral):     1/tau < DER <= tau
#     -> No correction needed; standard posterior is adequate
#   Tier III (survey-insulated):  DER <= 1/tau
#     -> No correction needed; shrinkage absorbs design effects
# =============================================================================

#' Classify parameters into three tiers based on DER
#'
#' @param der Named numeric vector of DER values.
#' @param tau Numeric scalar. Classification threshold (default: 1.2).
#'
#' @return A character vector with elements "sensitive", "neutral", or
#'   "insulated", corresponding to Tiers I, II, III.
classify_parameters <- function(der, tau = 1.2) {
  stopifnot(is.numeric(der), length(der) > 0, tau > 1)

  classification <- character(length(der))
  names(classification) <- names(der)

  classification[der > tau]                <- "sensitive"     # Tier I
  classification[der > 1/tau & der <= tau] <- "neutral"       # Tier II
  classification[der <= 1/tau]             <- "insulated"     # Tier III

  classification
}


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


#' Compute theoretical DER predictions from the decomposition
#'
#' The closed-form DER decomposition (Theorem 1, Section 3.2) predicts:
#'   Within-group FE:  DER approx DEFF
#'   Between-group FE: DER approx DEFF * (1 - B_bar)
#'   Random effects:   DER_j = B_j * DEFF_j * kappa_j(J)
#'
#' @param B_j Numeric vector of length J. Shrinkage factors.
#' @param deff_j Numeric vector of length J. Per-cluster DEFF.
#' @param J Integer. Number of clusters.
#' @param sigma_theta Numeric scalar. Random-effect SD.
#'
#' @return A list with theoretical DER predictions.
compute_theoretical_der <- function(B_j, deff_j, J, sigma_theta) {

  stopifnot(
    length(B_j) == J,
    length(deff_j) == J,
    J >= 2L,
    sigma_theta > 0
  )

  B_bar    <- mean(B_j)
  DEFF_bar <- mean(deff_j)

  kappa_j    <- (J - 1) * (1 - B_j) / (J * (1 - B_j) + B_j)
  der_re     <- B_j * deff_j * kappa_j

  der_fe_within  <- DEFF_bar
  der_fe_between <- DEFF_bar * (1 - B_bar)
  der_mu         <- DEFF_bar * (1 - B_bar)

  list(
    der_re_theory         = der_re,
    der_fe_within_theory  = der_fe_within,
    der_fe_between_theory = der_fe_between,
    der_mu_theory         = der_mu,
    B_bar                 = B_bar,
    DEFF_bar              = DEFF_bar
  )
}


# =============================================================================
# Section 10: Single-Replication Pipeline
# =============================================================================
#
# This function orchestrates the complete post-fitting pipeline for one
# replication: sandwich/DER computation, four correction strategies,
# credible intervals, and coverage evaluation. It is called from the
# parallel runner (sim_04_run.R).
# =============================================================================

#' Run the complete post-processing pipeline for one replication
#'
#' @param scenario Single-row data.frame or named list with scenario metadata.
#' @param rep_id Integer. Replication index.
#' @param stan_model Pre-compiled CmdStanModel object.
#' @param mcmc_settings List of MCMC parameters.
#' @param output_dir Character or NULL. If non-NULL, save result as .rds.
#' @param thresholds Numeric vector. DER thresholds (default c(1.2, 1.5)).
#' @param ci_level Numeric. Credible interval level (default 0.90).
#' @param beta_true Numeric vector. True fixed-effect coefficients.
#' @param beta_prior_sd Numeric. Prior SD for beta (default 5).
#' @param quiet Logical. Suppress Stan output.
#'
#' @return A list with: scenario_id, rep_id, der, coverage, ci_widths,
#'   convergence, timing, and other diagnostics.
run_single_rep <- function(scenario, rep_id, stan_model,
                           mcmc_settings = NULL, output_dir = NULL,
                           thresholds = c(1.2, 1.5), ci_level = 0.90,
                           beta_true = NULL, beta_prior_sd = 5,
                           quiet = TRUE) {

  t0     <- proc.time()["elapsed"]
  timing <- list()

  # -- Extract scenario parameters --------------------------------------------
  if (is.data.frame(scenario)) {
    scenario <- as.list(scenario[1, , drop = TRUE])
  }
  stopifnot(all(c("scenario_id", "J", "n_j", "icc", "cv_w",
                   "informative", "sigma_theta") %in% names(scenario)))

  scenario_id <- scenario$scenario_id
  J   <- as.integer(scenario$J)
  n_j <- as.integer(scenario$n_j)
  p   <- 3L
  d   <- p + J

  if (is.null(beta_true)) {
    beta_true <- if (exists("SIM_PARAMS")) SIM_PARAMS$beta_true else c(-0.5, 0.5, 0.3)
  }
  if (is.null(mcmc_settings)) {
    mcmc_settings <- if (exists("SIM_PARAMS")) SIM_PARAMS$mcmc else list(
      chains = 4L, warmup = 1000L, sampling = 1500L,
      adapt_delta = 0.90, max_treedepth = 12L
    )
  }

  # Reproducible seed: deterministic mapping from (scenario, rep_id) to seed.
  # Uses get_rep_seed() from sim_00_config.R, which computes:
  #   seed = base_seed + scenario_number * 10000 + rep_id
  # The scenario_number is the sorted row index in the 54-scenario grid.
  grid_for_seed <- build_scenario_grid()
  scenario_number <- which(grid_for_seed$scenario_id == scenario_id)
  if (length(scenario_number) == 0) {
    # Fallback: use a hash-based seed if scenario not found in grid
    scenario_number <- as.integer(abs(sum(utf8ToInt(scenario_id))) %% 54L) + 1L
  }
  seed <- get_rep_seed(scenario_number, rep_id)

  # Parameter types for coverage reporting (Section 3.3):
  #   beta_0 = intercept (between-group, confounded with theta)
  #   beta_1 = x_within  (within-group, orthogonal to theta)
  #   beta_2 = z_between (between-group, confounded with theta)
  #   theta_1..theta_J = random effects
  param_types <- c("fe_between", "fe_within", "fe_between", rep("re", J))

  # Null result skeleton for early-exit cases
  null_result <- list(
    scenario_id = scenario_id, rep_id = rep_id,
    convergence_failed = FALSE, fitting_failed = FALSE,
    der = NULL, der_beta = NULL, der_theta = NULL,
    B_j = NULL, deff_j = NULL, beta_hat = NULL,
    sigma_theta_hat = NULL, coverage = NULL, ci_widths = NULL,
    n_corrected = NULL, convergence = NULL, timing = NULL,
    data_diagnostics = NULL, param_types = param_types, seed = seed
  )

  # -- Step 1: Generate data --------------------------------------------------
  t1 <- proc.time()["elapsed"]
  data_gen <- tryCatch(
    generate_survey_data(
      J = J, n_j = n_j, icc = scenario$icc, beta = beta_true,
      cv_w = scenario$cv_w, informative = scenario$informative, seed = seed
    ),
    error = function(e) { warning("DGP failed: ", e$message); NULL }
  )
  timing$data_gen <- as.numeric(proc.time()["elapsed"] - t1)

  if (is.null(data_gen)) {
    null_result$fitting_failed <- TRUE
    null_result$timing <- timing
    return(null_result)
  }

  data_diag <- list(
    prevalence           = mean(data_gen$y),
    weight_cv_empirical  = data_gen$weight_cv_empirical,
    weight_deff_empirical = data_gen$weight_deff_empirical,
    theta_sd_empirical   = sd(data_gen$theta_true)
  )

  # -- Step 2: Fit Stan model -------------------------------------------------
  t2 <- proc.time()["elapsed"]
  fit_result <- fit_hlr(data_gen, stan_model, mcmc_settings,
                        seed = seed + 1000000L, quiet = quiet)
  timing$fitting <- as.numeric(proc.time()["elapsed"] - t2)

  if (!fit_result$succeeded) {
    null_result$fitting_failed <- TRUE
    null_result$timing <- timing
    null_result$data_diagnostics <- data_diag
    return(null_result)
  }

  # -- Step 3: Check convergence ----------------------------------------------
  if (!fit_result$convergence$passed) {
    null_result$convergence_failed <- TRUE
    null_result$convergence <- fit_result$convergence
    null_result$timing <- timing
    null_result$data_diagnostics <- data_diag
    return(null_result)
  }

  # -- Step 4: Compute sandwich variance and DER ------------------------------
  t4 <- proc.time()["elapsed"]
  sandwich_data <- list(
    y = data_gen$y, X = data_gen$X,
    group = data_gen$group, w = data_gen$weights
  )

  sandwich <- tryCatch(
    compute_sandwich_hlr(fit_result$fit, sandwich_data, beta_prior_sd),
    error = function(e) {
      warning("Sandwich failed: ", e$message); NULL
    }
  )
  timing$sandwich <- as.numeric(proc.time()["elapsed"] - t4)

  if (is.null(sandwich)) {
    null_result$fitting_failed <- TRUE
    null_result$convergence <- fit_result$convergence
    null_result$timing <- timing
    null_result$data_diagnostics <- data_diag
    return(null_result)
  }

  # -- Step 5: Apply correction strategies ------------------------------------
  t5 <- proc.time()["elapsed"]

  # Extract MCMC draws
  beta_names  <- paste0("beta[", seq_len(p), "]")
  theta_names <- paste0("theta[", seq_len(J), "]")
  all_param_names <- c(beta_names, theta_names)

  mcmc_draws <- fit_result$fit$draws(
    variables = all_param_names, format = "matrix"
  )
  stopifnot(ncol(mcmc_draws) == d)

  # True values for coverage evaluation
  true_values <- c(beta_true, data_gen$theta_true)
  names(true_values) <- all_param_names

  # Strategy 1: Naive (no correction)
  ci_naive  <- compute_credible_intervals(mcmc_draws, level = ci_level)
  cov_naive <- compute_coverage(true_values, ci_naive)

  # Strategy 2: Blanket Cholesky
  blanket_result <- tryCatch(
    apply_blanket_cholesky(mcmc_draws, sandwich),
    error = function(e) { warning("Blanket failed: ", e$message); NULL }
  )
  if (!is.null(blanket_result)) {
    ci_blanket  <- compute_credible_intervals(
      blanket_result$corrected_draws, level = ci_level
    )
    cov_blanket <- compute_coverage(true_values, ci_blanket)
  } else {
    ci_blanket  <- ci_naive
    cov_blanket <- cov_naive
  }

  # Strategies 3+: Selective DER-tau
  ci_selective  <- list()
  cov_selective <- list()
  sel_results   <- list()

  for (tau in thresholds) {
    tau_label <- paste0("DER-", tau)

    sel <- tryCatch(
      apply_selective_cholesky(mcmc_draws, sandwich, threshold = tau),
      error = function(e) {
        warning("Selective (tau=", tau, ") failed: ", e$message); NULL
      }
    )

    if (!is.null(sel)) {
      sel_results[[tau_label]] <- sel
      ci_selective[[tau_label]] <- compute_credible_intervals(
        sel$corrected_draws, level = ci_level
      )
      cov_selective[[tau_label]] <- compute_coverage(
        true_values, ci_selective[[tau_label]]
      )
    } else {
      ci_selective[[tau_label]]  <- ci_naive
      cov_selective[[tau_label]] <- cov_naive
    }
  }

  timing$correction <- as.numeric(proc.time()["elapsed"] - t5)

  # -- Step 6: Assemble results -----------------------------------------------
  coverage  <- list(naive = cov_naive, blanket = cov_blanket)
  ci_widths <- list(naive = ci_naive$width, blanket = ci_blanket$width)
  n_corrected <- c(naive = 0L,
                   blanket = ifelse(!is.null(blanket_result),
                                    blanket_result$n_corrected, 0L))

  for (tau_label in names(cov_selective)) {
    coverage[[tau_label]]  <- cov_selective[[tau_label]]
    ci_widths[[tau_label]] <- ci_selective[[tau_label]]$width
    n_corrected[tau_label] <- sel_results[[tau_label]]$n_corrected
  }

  timing$total <- as.numeric(proc.time()["elapsed"] - t0)

  result <- list(
    scenario_id        = scenario_id,
    rep_id             = rep_id,
    convergence_failed = FALSE,
    fitting_failed     = FALSE,
    der                = sandwich$der,
    der_beta           = sandwich$der_beta,
    der_theta          = sandwich$der_theta,
    B_j                = sandwich$B_j,
    deff_j             = sandwich$deff_j,
    beta_hat           = sandwich$beta_hat,
    sigma_theta_hat    = sandwich$sigma_theta_hat,
    coverage           = coverage,
    ci_widths          = ci_widths,
    n_corrected        = n_corrected,
    convergence        = fit_result$convergence,
    timing             = timing,
    data_diagnostics   = data_diag,
    param_types        = param_types,
    seed               = seed
  )

  # -- Step 7: Save to disk (if requested) ------------------------------------
  if (!is.null(output_dir)) {
    scenario_dir <- file.path(output_dir, scenario_id)
    if (!dir.exists(scenario_dir)) {
      dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
    }
    saveRDS(result, file.path(scenario_dir,
                              sprintf("rep_%04d.rds", rep_id)))
  }

  result
}


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
