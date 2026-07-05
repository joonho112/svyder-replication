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
#   2. classify_parameters()  -- Tier classification (paper taxonomy:
#                                I-a / I-b / II by parameter type, with
#                                Tier III reserved for hyperparameters)
#   3. selective_correct()    -- Selective correction of flagged parameters
#                                (block Cholesky or marginal rescale)
#   4. compute_theoretical_der() -- Closed-form DER predictions (Theorems 1-2
#      + the general matrix formula)
#
# This file is the SINGLE canonical home of these functions; do not
# re-define these names elsewhere in the package. The taxonomy implemented
# here is the one used by the application and the manuscript: tiers are
# assigned by PARAMETER TYPE (I-a within-cluster FE, I-b between-cluster FE,
# II random effects, III hyperparameters), and the correction decision is
# DER > tau.
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
## Section 2 : Tier Classification (Paper Taxonomy)
###############################################################################

#' Classify parameters into the paper's tier taxonomy
#'
#' Implements the tier scheme used by the application and the manuscript
#' (Paper Section 3.3 / Section 5). Tiers are assigned by PARAMETER TYPE:
#'
#'   Tier I-a (Survey-dominated):        within-cluster fixed effects
#'     -> identified from within-cluster contrasts; inherit the design
#'        effect (DER ~ DEFF); primary correction candidates.
#'   Tier I-b (Protected, between):      between-cluster fixed effects
#'     -> shielded by confounding with cluster means (DER ~ DEFF(1-B)).
#'   Tier II  (Protected, random):       random effects
#'     -> shielded by hierarchical shrinkage.
#'   Tier III (Hyperparameters):         variance components etc.
#'     -> reserved; the sandwich target for hyperparameters is outside the
#'        (beta, theta) parameter block treated by the DER framework.
#'
#' The correction decision is orthogonal to the tier: a parameter is
#' flagged for correction if and only if DER > tau (default 1.2).
#'
#' This type-based taxonomy is the only classifier used in the package.
#'
#' @param der         Numeric vector of DER values (named or unnamed).
#' @param param_types Character vector, same length as der. Each element
#'   one of "fe_within", "fe_between", "re", "hyper".
#' @param tau         Numeric scalar. Correction threshold (default 1.2).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{param_name}{Character. Parameter name (from names of der, or
#'       "param_1", "param_2", ... if unnamed).}
#'     \item{param_type}{Character. As supplied.}
#'     \item{der}{Numeric. DER value.}
#'     \item{tier}{Character. "I-a", "I-b", "II", or "III".}
#'     \item{tier_label}{Character. Descriptive label.}
#'     \item{flagged}{Logical. DER > tau.}
#'     \item{action}{Character. "CORRECT" or "retain".}
#'   }
classify_parameters <- function(der, param_types, tau = 1.2) {

  stopifnot(
    is.numeric(der),
    length(der) >= 1L,
    is.character(param_types),
    length(param_types) == length(der),
    all(param_types %in% c("fe_within", "fe_between", "re", "hyper")),
    is.numeric(tau), length(tau) == 1L, tau > 0
  )

  ## Parameter names
  if (is.null(names(der))) {
    param_names <- paste0("param_", seq_along(der))
  } else {
    param_names <- names(der)
  }

  tier_map <- c(
    fe_within  = "I-a",
    fe_between = "I-b",
    re         = "II",
    hyper      = "III"
  )
  label_map <- c(
    fe_within  = "Survey-dominated",
    fe_between = "Protected (between)",
    re         = "Protected (random effects)",
    hyper      = "Hyperparameter (reserved)"
  )

  flagged <- as.numeric(der) > tau

  data.frame(
    param_name = param_names,
    param_type = param_types,
    der        = as.numeric(der),
    tier       = unname(tier_map[param_types]),
    tier_label = unname(label_map[param_types]),
    flagged    = flagged,
    action     = ifelse(flagged, "CORRECT", "retain"),
    stringsAsFactors = FALSE
  )
}


###############################################################################
## Section 3 : Selective Correction (Block Cholesky / Marginal)
###############################################################################

#' Apply selective correction to posterior draws of flagged parameters
#'
#' Rescales the posterior draws of parameters flagged by DER > tau so that
#' their covariance matches the sandwich variance target, leaving unflagged
#' parameters bitwise unchanged. The correction preserves posterior means.
#'
#' Two methods are provided (this mirrors svyder::der_correct(), whose
#' block path is the canonical reference implementation):
#'
#'   "block_cholesky" (default): joint Cholesky rescaling of the flagged
#'     block. Writing Sigma_F = L1 L1^T and V_F = L2 L2^T for the MCMC and
#'     sandwich covariance of the flagged block F, each draw is mapped
#'     through
#'         phi* = phi_hat + L2 %*% solve(L1) %*% (phi - phi_hat),
#'     so the corrected block has covariance exactly V_F (marginal
#'     variances AND off-diagonal structure). This is the paper's
#'     Algorithm 1 (CCC) correction.
#'
#'   "marginal": per-parameter sqrt(DER_k) rescaling. Matches marginal
#'     variances only; the off-diagonal structure of the flagged block is
#'     NOT matched. It is retained, explicitly named, as a documented
#'     alternative to the block method.
#'
#' The two methods coincide exactly when exactly one parameter is flagged
#' (|F| = 1); with |F| > 1 (e.g., the design-PSU-level target in the NSECE
#' application) the block method is required for a valid joint correction.
#'
#' @param draws Numeric matrix (n_draws x d). Posterior draws.
#' @param der   Numeric vector of length d. DER values for each parameter.
#' @param V_sand Numeric matrix (d x d). Sandwich variance estimator.
#' @param sigma_mcmc Numeric matrix (d x d). MCMC posterior covariance.
#' @param tau   Numeric scalar. DER threshold above which correction is
#'   applied (default 1.20).
#' @param method "block_cholesky" (default) or "marginal".
#' @param point_est Optional numeric vector of length d. Centering point
#'   (posterior means). Defaults to colMeans(draws).
#'
#' @return A list with components:
#'   \describe{
#'     \item{draws_corrected}{Numeric matrix (n_draws x d).}
#'     \item{corrected_idx}{Integer vector. Indices of corrected parameters.}
#'     \item{n_corrected}{Integer. Number of parameters corrected.}
#'     \item{method}{Character. Method actually applied.}
#'     \item{scale_factors}{Numeric vector of length d. Marginal SD ratios
#'       sqrt(diag(V_sand)/diag(Sigma_MCMC)) for flagged parameters, 1
#'       elsewhere. The block method also matches these marginals exactly.}
#'     \item{threshold}{Numeric. The tau used.}
#'   }
#'
#' @details
#' The selective approach is critical: blanket correction (applying sandwich
#' inflation to ALL parameters) causes catastrophic over-inflation for
#' protected parameters (random effects), collapsing their coverage from
#' ~90% to ~20% (Paper Section 4, Figure 3).
selective_correct <- function(draws, der, V_sand, sigma_mcmc, tau = 1.20,
                              method = c("block_cholesky", "marginal"),
                              point_est = NULL) {

  method <- match.arg(method)

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

  ## Centering point: posterior means
  if (is.null(point_est)) {
    point_est <- colMeans(draws)
  }
  stopifnot(length(point_est) == d)

  diag_V    <- diag(V_sand)
  diag_mcmc <- diag(sigma_mcmc)

  ## Identify parameters to correct (DER > tau)
  correct_idx <- which(der > tau)
  n_corrected <- length(correct_idx)

  scale_factors   <- rep(1.0, d)
  draws_corrected <- draws

  if (n_corrected > 0L) {

    S <- correct_idx

    ## Marginal SD ratios (reported for both methods; the block method
    ## also matches these marginals exactly)
    scale_factors[S] <- sqrt(diag_V[S] / diag_mcmc[S])

    if (method == "marginal" || n_corrected == 1L) {
      ## Per-parameter rescale (identical to block method when |S| = 1)
      for (i in S) {
        draws_corrected[, i] <- point_est[i] +
          scale_factors[i] * (draws[, i] - point_est[i])
      }
    } else {
      ## Block Cholesky: phi* = mu + L2 %*% solve(L1) %*% (phi - mu)
      Sigma_S <- sigma_mcmc[S, S, drop = FALSE]
      V_S     <- V_sand[S, S, drop = FALSE]

      L1 <- tryCatch(
        t(chol(Sigma_S)),
        error = function(e) {
          stop("MCMC covariance of the flagged block is not positive ",
               "definite; cannot apply block correction.", call. = FALSE)
        }
      )
      L2 <- tryCatch(
        t(chol(V_S)),
        error = function(e) {
          if (!requireNamespace("Matrix", quietly = TRUE)) {
            stop("Sandwich covariance of the flagged block is not positive ",
                 "definite and the Matrix package is unavailable for a ",
                 "nearPD fallback.", call. = FALSE)
          }
          warning("Sandwich covariance of the flagged block is not positive ",
                  "definite (this can occur when the number of clusters is ",
                  "smaller than the flagged block); using the nearest ",
                  "positive-definite matrix.", call. = FALSE)
          t(chol(as.matrix(Matrix::nearPD(V_S)$mat)))
        }
      )

      A <- L2 %*% solve(L1)
      centered <- sweep(draws[, S, drop = FALSE], 2L, point_est[S])
      draws_corrected[, S] <- sweep(centered %*% t(A), 2L,
                                    point_est[S], FUN = "+")
    }
  }

  list(
    draws_corrected = draws_corrected,
    corrected_idx   = correct_idx,
    n_corrected     = n_corrected,
    method          = method,
    scale_factors   = scale_factors,
    threshold       = tau
  )
}


###############################################################################
## Section 4 : Theoretical DER Predictions (Theorems 1-2 + general formula)
###############################################################################

#' Compute theoretical DER predictions from closed-form formulas
#'
#' Evaluates the DER decomposition formulas derived in the paper (Section 3,
#' Theorems 1-2, plus the general matrix formula). These provide closed-form
#' predictions that can be compared against the brute-force DER computed via
#' Algorithm 1.
#'
#' Simple formulas (Theorems 1-2; exact for the BALANCED random-intercept
#' model with COMMON within-group design effect -- evaluated at per-group
#' (B_j, DEFF_j) they are the balanced/common-DEFF approximation, not exact;
#' the ensemble-exact unbalanced form is supplement S-B.3, eq.
#' (S-B: unbalanced-der)):
#'   - Within-cluster FE:  DER_beta_w = DEFF_bar
#'   - Between-cluster FE: DER_beta_b = DEFF_bar * (1 - B_bar)
#'   - Random effects:     DER_j ~= B_j * DEFF_j * kappa_j(J)
#'     where kappa_j = (J-1)(1-B_j) / (J(1-B_j) + B_j)
#'
#' General formula (main paper eq. (2.6) / supplement S-B.1, multi-covariate
#' model): computes DER via block matrix inversion, accounting for the
#' beta-theta coupling matrix C. This is the exact object under
#' heterogeneity; the simple formulas above are its homogeneous-case
#' specialization.
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
#'     \item{der_re_theory}{Numeric vector of length J. Simple formula DER
#'       for random effects.}
#'     \item{der_fe_within_theory}{Numeric scalar. Simple formula DER for
#'       within-group FE.}
#'     \item{der_fe_between_theory}{Numeric scalar. Simple formula DER for
#'       between-group FE.}
#'     \item{der_mu_theory}{Numeric scalar. Simple formula DER for the
#'       overall mean.}
#'     \item{B_bar}{Numeric scalar. Mean shrinkage factor.}
#'     \item{DEFF_bar}{Numeric scalar. Mean per-group design effect.}
#'     \item{der_re_general}{Numeric vector of length J. General formula DER
#'       for RE (only if H_obs/J_cluster provided).}
#'     \item{der_fe_general}{Numeric vector of length p. General formula DER
#'       for FE (only if H_obs/J_cluster provided).}
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
  ## Part A: Simple formulas (Theorems 1-2; balanced/common-DEFF
  ## approximation when evaluated at heterogeneous per-group inputs)
  ## =================================================================
  ## Field names use the *_theory convention consumed by the per-
  ## replication outputs (run_single_rep) and the Figure 2 script.

  ## Random effects: DER_j ~= B_j * DEFF_j * kappa_j
  ## kappa_j accounts for the finite-group correction (Paper Theorem 2).
  ## Exact only under common (B, DEFF); with heterogeneous per-group
  ## inputs this is the balanced/common-DEFF approximation (the exact
  ## ensemble form is supplement S-B.3).
  kappa_j <- (J - 1) * (1 - B_j) / (J * (1 - B_j) + B_j)
  result$der_re_theory <- B_j * deff_j * kappa_j

  ## Within-cluster FE: DER = DEFF (Paper Theorem 1, Eq. 10)
  ## Robust to covariates because group-mean centering eliminates
  ## the beta-theta coupling.
  result$der_fe_within_theory <- DEFF_bar

  ## Between-cluster FE: DER = DEFF * (1 - B_bar)
  ## (Paper Theorem 1, Eq. 11)
  result$der_fe_between_theory <- DEFF_bar * (1 - B_bar)

  ## Overall mean: DER_mu = DEFF * (1 - B_bar)
  result$der_mu_theory <- DEFF_bar * (1 - B_bar)

  ## Summary factors
  result$B_bar    <- B_bar
  result$DEFF_bar <- DEFF_bar

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

  ## Test classify_parameters (paper taxonomy: type-based tiers)
  der_test <- c(2.6, 0.03, 0.03, 0.9, 1.5, 0.5)
  names(der_test) <- paste0("param_", 1:6)
  types_test <- c("fe_within", "fe_between", "fe_between", "re", "re", "hyper")
  class_test <- classify_parameters(der_test, types_test, tau = 1.2)
  stopifnot(
    nrow(class_test) == 6,
    class_test$tier[1] == "I-a",
    class_test$tier[2] == "I-b",
    class_test$tier[4] == "II",
    class_test$tier[6] == "III",
    identical(class_test$flagged, der_test > 1.2),
    class_test$action[1] == "CORRECT",
    class_test$action[2] == "retain"
  )
  cat("  classify_parameters: OK\n")

  ## Test selective_correct (block Cholesky): corrected block covariance
  ## must equal V_sand[S,S]; unflagged columns must be untouched
  draws_test <- matrix(rnorm(5000 * d_test), 5000, d_test)
  der_sc <- c(rep(0.5, 5), rep(2.0, 5))  # last 5 flagged
  sc_result <- selective_correct(draws_test, der_sc, V_test, S_test,
                                 tau = 1.20, method = "block_cholesky")
  stopifnot(
    ncol(sc_result$draws_corrected) == d_test,
    nrow(sc_result$draws_corrected) == 5000,
    sc_result$n_corrected == 5,
    sc_result$method == "block_cholesky",
    identical(sc_result$draws_corrected[, 1:5], draws_test[, 1:5])
  )
  S_idx <- sc_result$corrected_idx
  emp_cov <- cov(sc_result$draws_corrected[, S_idx])
  ## Corrected covariance = A %*% Sigma_emp %*% t(A), which converges to
  ## V_sand[S,S] as the empirical covariance converges to sigma_mcmc[S,S].
  ## Verify the exact algebraic identity instead of the sampling limit:
  mu_S <- colMeans(draws_test)[S_idx]
  cen  <- sweep(draws_test[, S_idx], 2L, mu_S)
  L1 <- t(chol(S_test[S_idx, S_idx]))
  L2 <- t(chol(V_test[S_idx, S_idx]))
  A  <- L2 %*% solve(L1)
  target_cov <- A %*% cov(cen) %*% t(A)
  stopifnot(max(abs(emp_cov - target_cov)) < 1e-8)
  cat("  selective_correct (block_cholesky): OK\n")

  ## Marginal method: matches sqrt(DER) rescale per parameter
  sc_marg <- selective_correct(draws_test, der_sc, V_test, S_test,
                               tau = 1.20, method = "marginal")
  k1 <- S_idx[1]
  sf <- sqrt(diag(V_test)[k1] / diag(S_test)[k1])
  manual <- colMeans(draws_test)[k1] +
    sf * (draws_test[, k1] - colMeans(draws_test)[k1])
  stopifnot(max(abs(sc_marg$draws_corrected[, k1] - manual)) < 1e-10)
  cat("  selective_correct (marginal): OK\n")

  ## |S| = 1: block and marginal coincide exactly
  der_one <- c(rep(0.5, 9), 2.0)
  sc_b1 <- selective_correct(draws_test, der_one, V_test, S_test,
                             tau = 1.20, method = "block_cholesky")
  sc_m1 <- selective_correct(draws_test, der_one, V_test, S_test,
                             tau = 1.20, method = "marginal")
  stopifnot(identical(sc_b1$draws_corrected, sc_m1$draws_corrected))
  cat("  selective_correct (|S| = 1 coincidence): OK\n")

  ## Test compute_theoretical_der
  J_t <- 20
  B_t <- runif(J_t, 0.3, 0.9)
  deff_t <- runif(J_t, 1.0, 3.0)
  theory_test <- compute_theoretical_der(B_t, deff_t, J_t, sigma_theta = 0.5)
  stopifnot(
    length(theory_test$der_re_theory) == J_t,
    is.numeric(theory_test$der_fe_within_theory),
    is.numeric(theory_test$der_fe_between_theory),
    is.numeric(theory_test$der_mu_theory),
    abs(theory_test$der_fe_within_theory - mean(deff_t)) < 1e-12
  )
  cat("  compute_theoretical_der: OK\n")

  cat("\nAll self-tests passed.\n")
}
