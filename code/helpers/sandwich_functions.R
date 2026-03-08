# =============================================================================
# sandwich_functions.R: Sandwich Variance Estimation Functions
# =============================================================================
#
# Purpose : Build the observed information matrix (H_obs), clustered score
#           outer-product matrix (J_cluster), and sandwich variance estimator
#           (V_sand) for hierarchical logistic regression models fitted with
#           survey weights.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Contents:
#   1. build_H_obs_logistic()        -- Observed information (bread)
#   2. build_J_cluster()             -- Clustered meat matrix (group-level)
#   3. build_J_cluster_psu()         -- Clustered meat matrix (PSU-level)
#   4. compute_sandwich()            -- Full sandwich variance estimator
#   5. compute_mcmc_covariance()     -- Posterior covariance from MCMC draws
#   6. compute_group_deff()          -- Per-group Kish design effect
#   7. compute_shrinkage_factors()   -- Per-group shrinkage factors B_j
#
# Key equations (Paper Section 2.2-2.3):
#
#   H_obs = -d^2 ell_w / d phi d phi^T  (observed information, Eq. 5)
#         = H_data + H_prior
#         = sum_i w_i q_i(1-q_i) phi_i phi_i^T + diag(prior precisions)
#
#   J_cluster = sum_j s_j s_j^T         (clustered meat, Eq. 6)
#         where s_j = sum_{i in j} w_i (y_i - q_i) * phi_i
#
#   V_sand = H_obs^{-1} J_cluster H_obs^{-1}   (sandwich, Eq. 7)
#
# Parameter ordering: phi = (beta_1, ..., beta_p, theta_1, ..., theta_J)
# with dimension d = p + J.
#
# =============================================================================

library(Matrix)


###############################################################################
## Section 1 : Observed Information Matrix (H_obs)
###############################################################################

#' Build the observed information matrix for hierarchical logistic regression
#'
#' Constructs H_obs using efficient block structure (Paper Eq. 5):
#'
#'   H_obs = | H_bb   H_bt |   +   | diag(1/sigma_beta^2)   0                |
#'           | H_bt^T H_tt |       | 0                       diag(1/sigma_theta^2) |
#'
#' where:
#'   H_bb[a,b] = sum_i v_i * X[i,a] * X[i,b]       (p x p)
#'   H_bt[a,j] = sum_{i in j} v_i * X[i,a]          (p x J)
#'   H_tt[j,j] = sum_{i in j} v_i                   (J x J diagonal)
#'
#' and v_i = w_i * q_i * (1 - q_i) is the combined survey-working weight.
#'
#' @param X             Numeric matrix (N x p). Design matrix for fixed effects.
#' @param group         Integer vector (N). Group indices (1, ..., J).
#' @param w             Numeric vector (N). Normalized survey weights.
#' @param wt            Numeric vector (N). Working weights: q_i * (1 - q_i).
#' @param p             Integer. Number of fixed effects.
#' @param J             Integer. Number of groups.
#' @param N             Integer. Number of observations.
#' @param sigma_theta   Numeric scalar. Posterior mean of sigma_theta.
#' @param beta_prior_sd Numeric scalar. Prior SD for beta (default 5).
#'   Set to Inf for flat prior.
#'
#' @return Numeric matrix (d x d) where d = p + J.
build_H_obs_logistic <- function(X, group, w, wt, p, J, N,
                                 sigma_theta, beta_prior_sd = 5) {
  d <- p + J
  H <- matrix(0, d, d)

  ## Combined weight: v_i = w_i * q_i * (1 - q_i)
  v <- w * wt

  ## --- Beta-beta block (p x p): H_bb = X^T diag(v) X ---
  ## This is the standard weighted Fisher information for the FE.
  H_bb <- crossprod(X, X * v)
  H[1:p, 1:p] <- H_bb

  ## --- Beta-theta and theta-theta blocks ---
  ## For each group j:
  ##   H_bt[, j] = sum_{i in j} v_i * X[i, ]    (p-vector)
  ##   H_tt[j, j] = sum_{i in j} v_i             (scalar, diagonal)
  ## Cross-group theta entries are zero (groups are independent conditional
  ## on beta and sigma_theta).
  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    v_j   <- v[idx_j]

    ## Beta-theta coupling: reflects the correlation between FE and RE
    ## for group j. Strong coupling indicates the FE parameter absorbs
    ## information from group j (relevant for DER decomposition).
    if (length(idx_j) == 1L) {
      bt_j <- v_j * X[idx_j, ]
    } else {
      bt_j <- colSums(X[idx_j, , drop = FALSE] * v_j)
    }
    H[1:p, p + j] <- bt_j
    H[p + j, 1:p] <- bt_j

    ## Theta-theta diagonal: effective sample size for group j
    H[p + j, p + j] <- sum(v_j)
  }

  ## --- Prior contributions ---

  ## Beta prior: N(0, beta_prior_sd^2) => precision = 1/beta_prior_sd^2
  if (is.finite(beta_prior_sd)) {
    beta_prior_prec <- 1 / beta_prior_sd^2
    for (k in seq_len(p)) {
      H[k, k] <- H[k, k] + beta_prior_prec
    }
  }

  ## Theta prior: N(0, sigma_theta^2) => precision = 1/sigma_theta^2
  ## This is the shrinkage contribution from the hierarchical model.
  theta_prior_prec <- 1 / sigma_theta^2
  for (j in seq_len(J)) {
    H[p + j, p + j] <- H[p + j, p + j] + theta_prior_prec
  }

  ## Verify symmetry (should hold by construction)
  stopifnot(max(abs(H - t(H))) < 1e-10)

  H
}


###############################################################################
## Section 2 : Clustered Meat Matrix (Group-Level)
###############################################################################

#' Build the clustered score outer-product matrix (group-level clustering)
#'
#' Computes J_cluster = sum_j s_j s_j^T (Paper Eq. 6), where s_j is the
#' d-dimensional score sum for group j:
#'
#'   s_j = (sum_{i in j} r_i * X[i, ],  0, ..., 0, sum_{i in j} r_i, 0, ..., 0)
#'           [--- beta block, length p ---]  [--- theta block, entry j only ---]
#'
#' where r_i = w_i * (y_i - q_i) is the weighted score residual.
#'
#' The block-sparse structure of s_j reflects the fact that observation i
#' only contributes to the score for its own group's random effect theta_j.
#'
#' @param X     Numeric matrix (N x p). Design matrix.
#' @param group Integer vector (N). Group indices (1, ..., J).
#' @param r     Numeric vector (N). Weighted score residuals:
#'   r_i = w_i * (y_i - q_i).
#' @param p     Integer. Number of fixed effects.
#' @param J     Integer. Number of groups.
#' @param N     Integer. Number of observations.
#'
#' @return Numeric matrix (d x d) where d = p + J. Positive semi-definite
#'   by construction (sum of rank-1 outer products).
build_J_cluster <- function(X, group, r, p, J, N) {
  d <- p + J
  Jmat <- matrix(0, d, d)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)

    ## Score sum for group j: d-vector
    s_j <- numeric(d)

    ## Beta components: sum_{i in j} r_i * X[i, ]
    if (length(idx_j) == 1L) {
      s_j[1:p] <- r[idx_j] * X[idx_j, ]
    } else {
      s_j[1:p] <- colSums(X[idx_j, , drop = FALSE] * r[idx_j])
    }

    ## Theta_j component: sum_{i in j} r_i
    s_j[p + j] <- sum(r[idx_j])

    ## Outer product contribution: s_j s_j^T
    Jmat <- Jmat + tcrossprod(s_j)
  }

  ## Symmetrize (should be symmetric by construction; enforce numerically)
  Jmat <- (Jmat + t(Jmat)) / 2

  Jmat
}


###############################################################################
## Section 3 : Clustered Meat Matrix (PSU-Level)
###############################################################################

#' Build PSU-level clustered score outer-product matrix
#'
#' A sensitivity analysis variant of build_J_cluster() that uses primary
#' sampling units (PSUs) as the clustering level rather than the model's
#' group structure (e.g., states). This is relevant when the survey design's
#' PSUs do not align with the model's random-effect grouping variable.
#'
#' For each PSU g (which may span multiple states), the score sum is:
#'   s_g[1:p] = sum_{i in PSU g} r_i * X[i, ]
#'   s_g[p+j] = sum_{i in PSU g, state=j} r_i   for each state j in PSU g
#'
#' @param X        Numeric matrix (N x p). Design matrix.
#' @param state    Integer vector (N). State indices (1, ..., J_states).
#' @param psu      Integer vector (N). PSU indices (1, ..., G).
#' @param r        Numeric vector (N). Weighted score residuals.
#' @param p        Integer. Number of fixed effects.
#' @param J_states Integer. Number of states (random-effect groups).
#' @param N        Integer. Number of observations.
#'
#' @return Numeric matrix (d x d) where d = p + J_states.
build_J_cluster_psu <- function(X, state, psu, r, p, J_states, N) {

  d <- p + J_states
  G <- max(psu)

  stopifnot(
    nrow(X) == N,
    length(state) == N,
    length(psu)   == N,
    length(r)     == N,
    all(state >= 1L & state <= J_states),
    all(psu   >= 1L & psu   <= G)
  )

  Jmat <- matrix(0.0, d, d)

  for (g in seq_len(G)) {
    idx_g <- which(psu == g)
    if (length(idx_g) == 0L) next

    s_g <- numeric(d)

    ## Beta components: sum_{i in PSU g} r_i * X[i, ]
    if (length(idx_g) == 1L) {
      s_g[1:p] <- r[idx_g] * X[idx_g, ]
    } else {
      s_g[1:p] <- colSums(X[idx_g, , drop = FALSE] * r[idx_g])
    }

    ## Theta components: for each state that PSU g overlaps,
    ## accumulate r_i into the corresponding theta entry
    for (j_state in unique(state[idx_g])) {
      idx_gj <- idx_g[state[idx_g] == j_state]
      s_g[p + j_state] <- s_g[p + j_state] + sum(r[idx_gj])
    }

    ## Outer-product contribution
    Jmat <- Jmat + tcrossprod(s_g)
  }

  ## Enforce exact symmetry
  Jmat <- (Jmat + t(Jmat)) / 2.0

  Jmat
}


###############################################################################
## Section 4 : Full Sandwich Variance Estimator
###############################################################################

#' Compute the full sandwich variance estimator
#'
#' Assembles the sandwich variance V_sand = H_obs^{-1} J_cluster H_obs^{-1}
#' (Paper Eq. 7). Includes a nearPD fallback for near-singular H_obs.
#'
#' @param H_obs     Numeric matrix (d x d). Observed information matrix.
#' @param J_cluster Numeric matrix (d x d). Clustered score outer-product.
#'
#' @return A list with components:
#'   \describe{
#'     \item{V_sand}{Numeric matrix (d x d). Sandwich variance estimator.}
#'     \item{H_obs_inv}{Numeric matrix (d x d). Inverse of H_obs.}
#'   }
compute_sandwich <- function(H_obs, J_cluster) {

  stopifnot(
    is.matrix(H_obs),
    is.matrix(J_cluster),
    nrow(H_obs) == ncol(H_obs),
    nrow(J_cluster) == ncol(J_cluster),
    nrow(H_obs) == nrow(J_cluster)
  )

  d <- nrow(H_obs)

  ## Invert H_obs (with nearPD fallback for numerical stability)
  H_obs_inv <- tryCatch(
    solve(H_obs),
    error = function(e) {
      warning("H_obs is singular or near-singular. ",
              "Using nearPD fallback. ",
              "Condition number: ",
              tryCatch(kappa(H_obs, exact = TRUE),
                       error = function(e2) NA_real_))
      H_pd <- as.matrix(Matrix::nearPD(H_obs, keepDiag = TRUE)$mat)
      solve(H_pd)
    }
  )

  ## Sandwich: V_sand = H_inv J H_inv
  V_sand <- H_obs_inv %*% J_cluster %*% H_obs_inv

  ## Enforce symmetry
  V_sand <- (V_sand + t(V_sand)) / 2

  list(
    V_sand    = V_sand,
    H_obs_inv = H_obs_inv
  )
}


###############################################################################
## Section 5 : MCMC Posterior Covariance
###############################################################################

#' Compute posterior covariance matrix from MCMC draws
#'
#' Extracts posterior draws for all parameters (beta and theta), merges
#' across chains, and computes the sample covariance matrix Sigma_MCMC.
#' This serves as the denominator in the DER computation (Paper Eq. 8).
#'
#' @param fit A cmdstanr fit object with a draws() method.
#' @param p   Integer. Number of fixed effects.
#' @param J   Integer. Number of groups (random effects).
#'
#' @return Numeric matrix (d x d) where d = p + J. The posterior covariance
#'   matrix Sigma_MCMC = Cov(phi | y).
compute_mcmc_covariance <- function(fit, p, J) {
  d <- p + J

  ## Parameter names in Stan
  beta_names  <- paste0("beta[", seq_len(p), "]")
  theta_names <- paste0("theta[", seq_len(J), "]")
  all_names   <- c(beta_names, theta_names)

  ## Extract draws as matrix (iterations x parameters), merging all chains
  draws <- fit$draws(variables = all_names, format = "matrix")

  stopifnot(
    is.matrix(draws),
    ncol(draws) == d
  )

  ## Compute sample covariance
  sigma_mcmc <- cov(draws)

  ## Enforce exact symmetry
  sigma_mcmc <- (sigma_mcmc + t(sigma_mcmc)) / 2

  sigma_mcmc
}


###############################################################################
## Section 6 : Per-Group Design Effect (DEFF)
###############################################################################

#' Compute per-group Kish design effect (DEFF)
#'
#' For each group j, computes the Kish (1965) design effect:
#'
#'   DEFF_j = n_j * sum_{i in j} w_i^2 / (sum_{i in j} w_i)^2
#'
#' This measures the variance inflation due to unequal weighting within
#' each group. DEFF = 1 for equal weights; DEFF > 1 for unequal weights.
#' The DEFF enters the DER decomposition (Paper Theorem 2, Eq. 12).
#'
#' @param group Integer vector (N). Group indices (1, ..., J).
#' @param w     Numeric vector (N). Survey weights.
#' @param J     Integer. Number of groups.
#'
#' @return Numeric vector of length J with per-group DEFF values (all >= 1).
compute_group_deff <- function(group, w, J) {
  deff_j <- numeric(J)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    w_j   <- w[idx_j]
    n_j   <- length(idx_j)

    ## Kish DEFF: n * sum(w^2) / (sum(w))^2
    deff_j[j] <- n_j * sum(w_j^2) / (sum(w_j))^2
  }

  deff_j
}


###############################################################################
## Section 7 : Shrinkage Factors
###############################################################################

#' Compute per-group shrinkage factors B_j
#'
#' The shrinkage factor B_j measures how much the hierarchical model
#' pulls group j's random effect toward the grand mean (Paper Eq. 9):
#'
#'   B_j = sigma_theta^2 / (sigma_theta^2 + V_tilde_j)
#'
#' where V_tilde_j = 1 / sum_{i in j} w_i * q_i * (1 - q_i) is the
#' effective conditional variance from group j's data.
#'
#' B_j close to 1 means strong shrinkage (small/noisy groups).
#' B_j close to 0 means weak shrinkage (large/precise groups).
#'
#' The shrinkage factor appears in the DER decomposition:
#'   - For RE: DER_j = B_j * DEFF_j * kappa_j (Theorem 2)
#'   - For between-FE: DER_beta = DEFF * (1 - B_bar) (Theorem 1)
#'
#' @param group        Integer vector (N). Group indices (1, ..., J).
#' @param w            Numeric vector (N). Survey weights.
#' @param wt           Numeric vector (N). Working weights: q_i * (1 - q_i).
#' @param sigma_theta  Numeric scalar. Posterior mean of sigma_theta.
#' @param J            Integer. Number of groups.
#'
#' @return Numeric vector of length J with shrinkage factors in (0, 1).
compute_shrinkage_factors <- function(group, w, wt, sigma_theta, J) {
  sigma2 <- sigma_theta^2
  B_j <- numeric(J)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)

    ## Effective information from group j
    info_j <- sum(w[idx_j] * wt[idx_j])

    ## Conditional variance: inverse of effective information
    V_tilde_j <- 1 / info_j

    ## Shrinkage factor (Paper Eq. 9)
    B_j[j] <- sigma2 / (sigma2 + V_tilde_j)
  }

  stopifnot(all(B_j > 0), all(B_j < 1))

  B_j
}


###############################################################################
## Section 8 : Self-Test
###############################################################################

if (interactive()) {
  cat("=== sandwich_functions.R: Self-test ===\n\n")

  set.seed(123)
  N_test <- 100
  J_test <- 5
  p_test <- 3

  ## Synthetic data
  X_test     <- cbind(1, rnorm(N_test), rep(rnorm(J_test), each = N_test / J_test))
  group_test <- rep(1:J_test, each = N_test / J_test)
  w_test     <- runif(N_test, 0.5, 2.0)
  wt_test    <- runif(N_test, 0.1, 0.25)

  ## Test H_obs construction
  H_test <- build_H_obs_logistic(X_test, group_test, w_test, wt_test,
                                 p_test, J_test, N_test,
                                 sigma_theta = 0.5, beta_prior_sd = 5)
  d_test <- p_test + J_test
  stopifnot(
    nrow(H_test) == d_test,
    ncol(H_test) == d_test,
    max(abs(H_test - t(H_test))) < 1e-12
  )
  cat("  build_H_obs_logistic: OK (", d_test, "x", d_test, ")\n")

  ## Test J_cluster construction
  r_test <- rnorm(N_test)
  Jmat_test <- build_J_cluster(X_test, group_test, r_test,
                               p_test, J_test, N_test)
  stopifnot(
    nrow(Jmat_test) == d_test,
    ncol(Jmat_test) == d_test,
    max(abs(Jmat_test - t(Jmat_test))) < 1e-12
  )
  eigen_vals <- eigen(Jmat_test, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(all(eigen_vals >= -1e-10))
  cat("  build_J_cluster: OK (PSD verified)\n")

  ## Test compute_sandwich
  sand_test <- compute_sandwich(H_test, Jmat_test)
  stopifnot(
    nrow(sand_test$V_sand) == d_test,
    ncol(sand_test$V_sand) == d_test
  )
  cat("  compute_sandwich: OK\n")

  ## Test compute_group_deff
  deff_test <- compute_group_deff(group_test, w_test, J_test)
  stopifnot(
    length(deff_test) == J_test,
    all(deff_test >= 1)
  )
  cat("  compute_group_deff: OK (all >= 1)\n")

  ## Test compute_shrinkage_factors
  B_test <- compute_shrinkage_factors(group_test, w_test, wt_test,
                                      sigma_theta = 0.5, J_test)
  stopifnot(
    length(B_test) == J_test,
    all(B_test > 0),
    all(B_test < 1)
  )
  cat("  compute_shrinkage_factors: OK (all in (0,1))\n")

  cat("\nAll self-tests passed.\n")
}
