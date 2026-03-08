# =============================================================================
# der_functions.R: Core DER Computation and Classification Functions
# =============================================================================
#
# Purpose : Compute Design Effect Ratios (DER) from sandwich variance and
#           posterior covariance matrices. Implements the three-tier
#           classification scheme and selective Cholesky correction.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Contents:
#   1. compute_der()          -- Main DER computation from matrices
#   2. classify_parameters()  -- Three-tier DER classification
#   3. selective_correct()    -- Selective Cholesky correction (Algorithm 2)
#   4. compute_theoretical_der() -- Closed-form DER predictions (Theorems 1-3)
#
# Key equations (Paper Section 2.3-3.2):
#   DER_p = diag(V_sand)[p] / diag(Sigma_MCMC)[p]
#   V_sand = H_obs^{-1} J_cluster H_obs^{-1}   (sandwich variance)
#   Sigma_MCMC = Cov(phi | y)                    (posterior covariance from MCMC)
#
# =============================================================================

library(Matrix)


###############################################################################
## Section 1 : Main DER Computation
###############################################################################

#' Compute Design Effect Ratios (DER) from variance matrices
#'
#' Computes the parameter-specific Design Effect Ratio (DER), defined as the
#' ratio of the sandwich (design-adjusted) variance to the model-based
#' posterior variance (Paper Eq. 8):
#'
#'   DER_p = diag(V_sand)[p] / diag(Sigma_MCMC)[p]
#'
#' A DER of 1.0 indicates that the model-based variance already accounts for
#' the survey design. DER > 1 indicates survey-induced under-coverage;
#' DER < 1 indicates the hierarchical model provides more shrinkage than
#' the design effect would suggest.
#'
#' @param V_sand    Numeric matrix (d x d). Sandwich variance estimator.
#'   Computed as H_obs^{-1} %*% J_cluster %*% H_obs^{-1}.
#' @param sigma_mcmc Numeric matrix (d x d). Posterior covariance from MCMC
#'   draws, i.e., Cov(phi | y) estimated from the merged posterior chains.
#' @param H_obs_inv Numeric matrix (d x d). Optional. Inverse of the observed
#'   information matrix. If provided, also computes the Laplace-based DER
#'   (V_sand / H_obs_inv) as a secondary diagnostic.
#' @param param_names Character vector of length d. Optional. Names for each
#'   parameter (e.g., "beta[1]", "theta[1]", ...).
#'
#' @return A list with components:
#'   \describe{
#'     \item{der}{Named numeric vector of length d. Primary DER values
#'       (V_sand / Sigma_MCMC).}
#'     \item{der_laplace}{Named numeric vector of length d or NULL. Laplace-
#'       based DER values (V_sand / H_obs_inv), only if H_obs_inv is provided.}
#'     \item{diag_V}{Numeric vector of length d. Diagonal of V_sand.}
#'     \item{diag_mcmc}{Numeric vector of length d. Diagonal of Sigma_MCMC.}
#'   }
#'
#' @details
#' The DER definition uses the MCMC posterior covariance as the denominator
#' rather than the Laplace approximation (H_obs_inv). This ensures consistency
#' with the Cholesky correction framework (Algorithm 2), which transforms
#' MCMC draws from Sigma_MCMC to V_sand. For the normal model, the Laplace
#' and MCMC denominators coincide exactly.
compute_der <- function(V_sand, sigma_mcmc, H_obs_inv = NULL,
                        param_names = NULL) {

  ## --- Input validation ---
  stopifnot(
    is.matrix(V_sand),
    is.matrix(sigma_mcmc),
    nrow(V_sand) == ncol(V_sand),
    nrow(sigma_mcmc) == ncol(sigma_mcmc),
    nrow(V_sand) == nrow(sigma_mcmc)
  )

  d <- nrow(V_sand)

  ## Extract diagonals
  diag_V    <- diag(V_sand)
  diag_mcmc <- diag(sigma_mcmc)

  ## Validate positive diagonals
  if (any(diag_V <= 0)) {
    warning(sprintf("V_sand has %d non-positive diagonal entries.",
                    sum(diag_V <= 0)))
  }
  stopifnot("sigma_mcmc diagonal must be all positive" = all(diag_mcmc > 0))

  ## Primary DER: V_sand / Sigma_MCMC (Paper Eq. 8)
  der <- diag_V / diag_mcmc

  ## Laplace-based DER (secondary, for cross-validation)
  der_laplace <- NULL
  if (!is.null(H_obs_inv)) {
    stopifnot(
      is.matrix(H_obs_inv),
      nrow(H_obs_inv) == d,
      ncol(H_obs_inv) == d
    )
    diag_H <- diag(H_obs_inv)
    stopifnot(all(diag_H > 0))
    der_laplace <- diag_V / diag_H
  }

  ## Apply names
  if (!is.null(param_names)) {
    stopifnot(length(param_names) == d)
    names(der) <- param_names
    if (!is.null(der_laplace)) names(der_laplace) <- param_names
  }

  list(
    der         = der,
    der_laplace = der_laplace,
    diag_V      = diag_V,
    diag_mcmc   = diag_mcmc
  )
}


###############################################################################
## Section 2 : Three-Tier DER Classification
###############################################################################

#' Classify parameters into three DER-based tiers
#'
#' Implements the three-tier classification scheme (Paper Section 3.2,
#' Table 1):
#'
#'   Tier I   (Survey-Robust):     DER < lo
#'     -> No correction needed; hierarchical shrinkage absorbs design effect.
#'   Tier II  (Survey-Sensitive):  lo <= DER <= hi
#'     -> Borderline; monitor but correction may not improve coverage.
#'   Tier III (Survey-Inflated):   DER > hi
#'     -> Correction needed; sandwich variance exceeds model-based variance.
#'
#' The default thresholds (lo = 0.80, hi = 1.20) correspond to a 20%
#' tolerance band around DER = 1.0. The paper's simulation study (Section 4)
#' validates that these thresholds achieve good coverage properties.
#'
#' @param der Numeric vector. DER values (named or unnamed).
#' @param lo  Numeric scalar. Lower threshold (default 0.80).
#' @param hi  Numeric scalar. Upper threshold (default 1.20).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{parameter}{Character. Parameter name (from names of der, or
#'       "param_1", "param_2", etc. if unnamed).}
#'     \item{der}{Numeric. DER value.}
#'     \item{tier}{Character. Tier label: "Tier I", "Tier II", or "Tier III".}
#'     \item{tier_label}{Character. Full descriptive label.}
#'     \item{action}{Character. Recommended action.}
#'   }
classify_parameters <- function(der, lo = 0.80, hi = 1.20) {

  stopifnot(
    is.numeric(der),
    length(der) >= 1L,
    is.numeric(lo), is.numeric(hi),
    lo > 0, hi > lo
  )

  ## Parameter names
  if (is.null(names(der))) {
    param_names <- paste0("param_", seq_along(der))
  } else {
    param_names <- names(der)
  }

  ## Classification
  tier <- character(length(der))
  tier_label <- character(length(der))
  action <- character(length(der))

  ## Tier I: Survey-Robust (DER < lo)
  idx_i <- der < lo
  tier[idx_i]       <- "Tier I"
  tier_label[idx_i] <- "Survey-Robust"
  action[idx_i]     <- "No correction needed"

  ## Tier II: Survey-Sensitive (lo <= DER <= hi)
  idx_ii <- der >= lo & der <= hi
  tier[idx_ii]       <- "Tier II"
  tier_label[idx_ii] <- "Survey-Sensitive"
  action[idx_ii]     <- "Monitor; correction optional"

  ## Tier III: Survey-Inflated (DER > hi)
  idx_iii <- der > hi
  tier[idx_iii]       <- "Tier III"
  tier_label[idx_iii] <- "Survey-Inflated"
  action[idx_iii]     <- "Correction recommended"

  data.frame(
    parameter  = param_names,
    der        = as.numeric(der),
    tier       = tier,
    tier_label = tier_label,
    action     = action,
    stringsAsFactors = FALSE
  )
}


###############################################################################
## Section 3 : Selective Cholesky Correction (Algorithm 2)
###############################################################################

#' Apply selective Cholesky correction to posterior draws
#'
#' Implements Algorithm 2 from the paper (Section 3.3). For parameters
#' classified as Tier III (DER > tau), the posterior draws are transformed
#' so that their marginal variance matches the sandwich variance. Parameters
#' in Tier I and II are left unchanged.
#'
#' The correction proceeds as follows:
#'   1. Partition parameters into S (to correct, DER > tau) and S^c (to keep).
#'   2. Compute the target covariance for the corrected subset by inflating
#'      Sigma_MCMC[S,S] to match V_sand[S,S], while preserving the
#'      correlation structure.
#'   3. Apply the affine transformation:
#'        phi_corrected[S] = mu[S] + L_target %*% L_mcmc_inv %*% (phi[S] - mu[S])
#'      where L_target and L_mcmc are Cholesky factors of the respective
#'      covariance submatrices.
#'
#' @param draws Numeric matrix (n_draws x d). Posterior draws, one row per
#'   iteration, columns corresponding to the d model parameters.
#' @param der   Numeric vector of length d. DER values for each parameter.
#' @param V_sand Numeric matrix (d x d). Sandwich variance estimator.
#' @param sigma_mcmc Numeric matrix (d x d). MCMC posterior covariance.
#' @param tau   Numeric scalar. DER threshold above which correction is
#'   applied (default 1.20, matching the Tier III cutoff).
#'
#' @return A list with components:
#'   \describe{
#'     \item{draws_corrected}{Numeric matrix (n_draws x d). Corrected
#'       posterior draws.}
#'     \item{corrected_idx}{Integer vector. Column indices of corrected
#'       parameters.}
#'     \item{n_corrected}{Integer. Number of parameters corrected.}
#'   }
#'
#' @details
#' The selective approach is critical: blanket correction (applying sandwich
#' inflation to ALL parameters) causes catastrophic over-inflation for
#' Tier I parameters (random effects), collapsing their coverage from ~90%
#' to ~20% (Paper Section 4, Figure 3).
selective_correct <- function(draws, der, V_sand, sigma_mcmc, tau = 1.20) {

  ## --- Input validation ---
  stopifnot(
    is.matrix(draws),
    is.numeric(der),
    ncol(draws) == length(der),
    is.matrix(V_sand),
    is.matrix(sigma_mcmc),
    nrow(V_sand) == length(der),
    nrow(sigma_mcmc) == length(der),
    is.numeric(tau), tau > 0
  )

  d <- length(der)
  n_draws <- nrow(draws)

  ## Identify parameters to correct (Tier III: DER > tau)
  correct_idx <- which(der > tau)
  n_corrected <- length(correct_idx)

  ## If no parameters need correction, return original draws
  if (n_corrected == 0L) {
    return(list(
      draws_corrected = draws,
      corrected_idx   = integer(0),
      n_corrected     = 0L
    ))
  }

  ## --- Extract submatrices for the correction subset S ---
  S <- correct_idx
  Sigma_S <- sigma_mcmc[S, S, drop = FALSE]
  V_S     <- V_sand[S, S, drop = FALSE]

  ## Posterior means (column means of draws)
  mu <- colMeans(draws)

  ## --- Cholesky decomposition ---
  ## L_mcmc: Cholesky factor of Sigma_MCMC[S,S]
  ## L_target: Cholesky factor of V_sand[S,S]
  ##
  ## For numerical stability, project to nearest PD if needed.
  L_mcmc <- tryCatch(
    chol(Sigma_S),
    error = function(e) {
      chol(as.matrix(Matrix::nearPD(Sigma_S, keepDiag = TRUE)$mat))
    }
  )

  L_target <- tryCatch(
    chol(V_S),
    error = function(e) {
      chol(as.matrix(Matrix::nearPD(V_S, keepDiag = TRUE)$mat))
    }
  )

  ## --- Affine transformation ---
  ## For each draw: phi_new[S] = mu[S] + L_target^T %*% (L_mcmc^T)^{-1}
  ##                              %*% (phi[S] - mu[S])
  ##
  ## Using upper Cholesky: R = chol(Sigma) means Sigma = R^T R
  ## We need: L_target^T %*% (L_mcmc^T)^{-1} = t(L_target) %*% solve(t(L_mcmc))
  ##        = (solve(L_mcmc) %*% L_target)^T
  ##
  ## More efficiently: transform = solve(L_mcmc, L_target)
  ## Then: phi_new[S] = mu[S] + (phi[S] - mu[S]) %*% transform

  transform_mat <- solve(L_mcmc, L_target)

  ## Center, transform, and re-center
  draws_corrected <- draws
  centered_S <- draws[, S, drop = FALSE] - matrix(mu[S], n_draws,
                                                    n_corrected,
                                                    byrow = TRUE)
  transformed_S <- centered_S %*% transform_mat
  draws_corrected[, S] <- transformed_S + matrix(mu[S], n_draws,
                                                   n_corrected,
                                                   byrow = TRUE)

  list(
    draws_corrected = draws_corrected,
    corrected_idx   = correct_idx,
    n_corrected     = n_corrected
  )
}


###############################################################################
## Section 4 : Theoretical DER Predictions (Theorems 1-3)
###############################################################################

#' Compute theoretical DER predictions from closed-form formulas
#'
#' Evaluates the DER decomposition formulas derived in the paper (Section 3.1,
#' Theorems 1-3). These provide closed-form predictions that can be compared
#' against the brute-force DER computed via Algorithm 1.
#'
#' Simple formulas (Theorems 1-2, exact for random-intercept model):
#'   - Within-cluster FE:  DER_beta_w = DEFF_bar
#'   - Between-cluster FE: DER_beta_b = DEFF_bar * (1 - B_bar)
#'   - Random effects:     DER_j = B_j * DEFF_j * kappa_j(J)
#'     where kappa_j = (J-1)(1-B_j) / (J(1-B_j) + B_j)
#'
#' General formula (Theorem 3, multi-covariate model):
#'   Computes DER via block matrix inversion, accounting for the beta-theta
#'   coupling matrix C.
#'
#' @param B_j     Numeric vector of length J. Per-group shrinkage factors.
#' @param deff_j  Numeric vector of length J. Per-group Kish design effects.
#' @param J       Integer. Number of groups.
#' @param sigma_theta Numeric scalar. Random-effect SD.
#' @param H_obs   Optional (p+J) x (p+J) observed information matrix.
#'   Required for the general formula.
#' @param J_cluster Optional (p+J) x (p+J) clustered score matrix.
#'   Required for the general formula.
#' @param p       Optional integer. Number of fixed effects. Required with
#'   H_obs and J_cluster.
#'
#' @return A list with components:
#'   \describe{
#'     \item{der_re_simple}{Numeric vector of length J. Simple formula DER
#'       for random effects.}
#'     \item{der_fe_within_simple}{Numeric scalar. Simple formula DER for
#'       within-group FE.}
#'     \item{der_fe_between_simple}{Numeric scalar. Simple formula DER for
#'       between-group FE.}
#'     \item{der_re_general}{Numeric vector of length J. General formula DER
#'       for RE (if H_obs/J_cluster provided).}
#'     \item{der_fe_general}{Numeric vector of length p. General formula DER
#'       for FE (if H_obs/J_cluster provided).}
#'   }
compute_theoretical_der <- function(B_j, deff_j, J, sigma_theta,
                                    H_obs = NULL, J_cluster = NULL,
                                    p = NULL) {

  stopifnot(
    length(B_j) == J,
    length(deff_j) == J,
    J >= 2L,
    sigma_theta > 0
  )

  result <- list()
  B_bar    <- mean(B_j)
  DEFF_bar <- mean(deff_j)

  ## =================================================================
  ## Part A: Simple formulas (Theorems 1-2)
  ## =================================================================

  ## Random effects: DER_j = B_j * DEFF_j * kappa_j
  ## kappa_j accounts for the finite-group correction
  ## (Paper Theorem 2, Eq. 12)
  kappa_j <- (J - 1) * (1 - B_j) / (J * (1 - B_j) + B_j)
  result$der_re_simple <- B_j * deff_j * kappa_j

  ## Within-cluster FE: DER = DEFF (Paper Theorem 1, Eq. 10)
  ## Robust to covariates because group-mean centering eliminates
  ## the beta-theta coupling.
  result$der_fe_within_simple <- DEFF_bar

  ## Between-cluster FE: DER = DEFF * (1 - B_bar)
  ## (Paper Theorem 1, Eq. 11)
  result$der_fe_between_simple <- DEFF_bar * (1 - B_bar)

  ## Overall mean: DER_mu = DEFF * (1 - B_bar)
  result$der_mu_simple <- DEFF_bar * (1 - B_bar)

  ## =================================================================
  ## Part B: General formula (Theorem 3)
  ## =================================================================

  if (!is.null(H_obs) && !is.null(J_cluster) && !is.null(p)) {
    d <- p + J
    stopifnot(
      nrow(H_obs) == d, ncol(H_obs) == d,
      nrow(J_cluster) == d, ncol(J_cluster) == d
    )

    ## Full sandwich computation
    H_inv <- tryCatch(solve(H_obs), error = function(e) {
      as.matrix(solve(Matrix::nearPD(H_obs, keepDiag = TRUE)$mat))
    })

    V_sand_gen <- H_inv %*% J_cluster %*% H_inv
    V_sand_gen <- (V_sand_gen + t(V_sand_gen)) / 2

    ## General DER (Laplace-based for theoretical comparison)
    result$der_fe_general <- diag(V_sand_gen)[1:p] / diag(H_inv)[1:p]
    result$der_re_general <- diag(V_sand_gen)[(p + 1):d] / diag(H_inv)[(p + 1):d]
  }

  result
}


###############################################################################
## Section 5 : Self-Test
###############################################################################

if (interactive()) {
  cat("=== der_functions.R: Self-test ===\n\n")

  ## Test compute_der with synthetic matrices
  set.seed(42)
  d_test <- 10
  A <- matrix(rnorm(d_test^2), d_test, d_test)
  V_test <- crossprod(A)
  B <- matrix(rnorm(d_test^2), d_test, d_test)
  S_test <- crossprod(B)

  result_test <- compute_der(V_test, S_test)
  stopifnot(
    length(result_test$der) == d_test,
    all(result_test$der > 0),
    is.null(result_test$der_laplace)
  )
  cat("  compute_der: OK\n")

  ## Test classify_parameters
  der_test <- c(0.5, 0.9, 1.0, 1.1, 1.5, 2.0)
  names(der_test) <- paste0("param_", 1:6)
  class_test <- classify_parameters(der_test)
  stopifnot(
    nrow(class_test) == 6,
    class_test$tier[1] == "Tier I",
    class_test$tier[3] == "Tier II",
    class_test$tier[6] == "Tier III"
  )
  cat("  classify_parameters: OK\n")

  ## Test selective_correct
  draws_test <- matrix(rnorm(1000 * d_test), 1000, d_test)
  der_sc <- c(rep(0.5, 5), rep(2.0, 5))  # 5 Tier I + 5 Tier III
  sc_result <- selective_correct(draws_test, der_sc, V_test, S_test, tau = 1.20)
  stopifnot(
    ncol(sc_result$draws_corrected) == d_test,
    nrow(sc_result$draws_corrected) == 1000,
    sc_result$n_corrected == 5
  )
  cat("  selective_correct: OK\n")

  ## Test compute_theoretical_der
  J_t <- 20
  B_t <- runif(J_t, 0.3, 0.9)
  deff_t <- runif(J_t, 1.0, 3.0)
  theory_test <- compute_theoretical_der(B_t, deff_t, J_t, sigma_theta = 0.5)
  stopifnot(
    length(theory_test$der_re_simple) == J_t,
    is.numeric(theory_test$der_fe_within_simple),
    is.numeric(theory_test$der_fe_between_simple)
  )
  cat("  compute_theoretical_der: OK\n")

  cat("\nAll self-tests passed.\n")
}
