# =============================================================================
# sim_01_dgp.R: Data Generating Process
# =============================================================================
#
# Purpose : Implement the data generating process (DGP) for the simulation
#           study. Generates clustered binary survey data from a hierarchical
#           logistic regression model with either non-informative or
#           informative survey weights.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Section 4.1 (Data Generating Process)
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Track   : A (Full Replication)
# Inputs  : sim_00_config.R (for icc_to_sigma2, SIM_PARAMS)
# Outputs : generate_survey_data() function returning a list suitable for
#           Stan model fitting
#
# Model (Section 4.1, Eq. 15-16):
#   y_ij | p_ij ~ Bernoulli(p_ij)
#   logit(p_ij) = beta_0 + beta_1 * x_wc_ij + beta_2 * z_bc_j + theta_j
#   theta_j ~ N(0, sigma^2_theta)
#
# Covariates:
#   x_wc_ij : within-cluster covariate (group-mean centered, so purely
#             within-cluster variation). This ensures beta_1 is a
#             within-group parameter unconfounded with theta_j.
#   z_bc_j  : between-cluster covariate (constant within clusters).
#             beta_2 is a between-group parameter confounded with theta_j.
#
# Weights (Section 4.1, Eq. 17-18):
#   Non-informative: log(w_ij) ~ N(0, sigma^2_w)
#     where sigma^2_w = log(1 + CV_w^2) so that CV(w) = CV_w.
#   Informative:     log(w_ij) = alpha * theta_j + epsilon_ij
#     where alpha and sigma_eps are calibrated so that:
#       (a) CV(w) ~ CV_w (matching the target design effect)
#       (b) Cor(theta_j, mean_j(log w)) ~ 0.3 (moderate informativeness)
# =============================================================================


# =============================================================================
# Section 1: Main Data Generation Function
# =============================================================================

#' Generate survey data from a hierarchical logistic model
#'
#' This is the core DGP function. It generates all components of one simulated
#' dataset: random effects, covariates, binary outcomes, and survey weights.
#'
#' @param J Integer. Number of clusters (PSUs).
#' @param n_j Integer. Number of observations per cluster (balanced).
#' @param icc Numeric in (0,1). Target latent-scale ICC.
#' @param beta Numeric vector of length 3. True coefficients:
#'   c(intercept, within-cluster effect, between-cluster effect).
#' @param cv_w Numeric > 0. Target coefficient of variation of survey weights.
#' @param informative Logical. If TRUE, weights depend on random effects.
#' @param seed Integer or NULL. Random seed for reproducibility.
#'
#' @return A list with all data needed for model fitting and evaluation:
#'   y, X, group, weights, theta_true, sigma2_theta, sigma_theta,
#'   weight_cv_empirical, weight_deff_empirical, and scenario metadata.
generate_survey_data <- function(J, n_j, icc, beta, cv_w,
                                 informative, seed = NULL) {

  # -- Input validation -------------------------------------------------------
  stopifnot(
    is.numeric(J), length(J) == 1, J >= 2, J == as.integer(J),
    is.numeric(n_j), length(n_j) == 1, n_j >= 2, n_j == as.integer(n_j),
    is.numeric(icc), length(icc) == 1, icc > 0, icc < 1,
    is.numeric(beta), length(beta) == 3,
    is.numeric(cv_w), length(cv_w) == 1, cv_w > 0,
    is.logical(informative), length(informative) == 1
  )

  J   <- as.integer(J)
  n_j <- as.integer(n_j)
  N   <- J * n_j

  if (!is.null(seed)) set.seed(seed)

  # -- Step 1: Compute random-effect variance from ICC (Eq. 8) ----------------
  # sigma^2_theta = ICC * pi^2 / (3 * (1 - ICC))
  sigma2_theta <- icc_to_sigma2(icc)
  sigma_theta  <- sqrt(sigma2_theta)

  # -- Step 2: Generate between-cluster covariate -----------------------------
  # z_bc_j ~ N(0, 1), constant within each cluster j.
  # This covariate is confounded with the random intercept theta_j,
  # making beta_2 a "between-group" parameter.
  z_j <- rnorm(J, mean = 0, sd = 1)

  # -- Step 3: Generate cluster random effects --------------------------------
  # theta_j ~ N(0, sigma^2_theta)
  # These represent unobserved cluster-level heterogeneity.
  theta_j <- rnorm(J, mean = 0, sd = sigma_theta)

  # -- Step 4: Expand cluster-level quantities to observation level -----------
  group_idx <- rep(seq_len(J), each = n_j)
  z_ij      <- z_j[group_idx]
  theta_ij  <- theta_j[group_idx]

  # -- Step 5: Generate within-cluster covariate (group-mean centered) --------
  # Raw covariate: x_raw_ij ~ N(0, 1)
  # Group-mean centering: x_wc_ij = x_raw_ij - mean_j(x_raw)
  # After centering, x_wc has zero between-cluster variation, ensuring
  # that beta_1 is orthogonal to theta_j (a "within-group" parameter).
  # This is critical for the DER tier classification (Section 3.3):
  #   - beta_1 (within-group) should have DER ~ DEFF (Tier I-a)
  #   - beta_0, beta_2 (between-group) should have DER < DEFF (Tier I-b)
  x_raw <- rnorm(N, mean = 0, sd = 1)
  group_means <- tapply(x_raw, group_idx, mean)
  x_ij <- x_raw - group_means[group_idx]

  # -- Step 6: Assemble design matrix -----------------------------------------
  # X = [1, x_wc, z_bc] : N x 3 matrix
  X <- cbind(
    intercept = rep(1, N),
    x_within  = x_ij,
    z_between = z_ij
  )

  # -- Step 7: Generate binary outcomes (Eq. 15-16) ---------------------------
  # eta_ij = X_ij %*% beta + theta_j
  # p_ij = logit^{-1}(eta_ij)
  # y_ij ~ Bernoulli(p_ij)
  eta <- X %*% beta + theta_ij
  p   <- plogis(as.numeric(eta))
  y   <- rbinom(N, size = 1, prob = p)

  # -- Step 8: Generate survey weights (Eq. 17-18) ----------------------------
  wt_result <- generate_survey_weights(
    J           = J,
    n_j         = n_j,
    cv_w        = cv_w,
    informative = informative,
    theta_j     = theta_j,
    seed        = NULL   # seed already set upstream; continue the RNG stream
  )

  # -- Step 9: Assemble and return the complete dataset -----------------------
  list(
    y                    = y,
    X                    = X,
    group                = group_idx,
    weights              = wt_result$weights,
    theta_true           = theta_j,
    z_j                  = z_j,
    sigma2_theta         = sigma2_theta,
    sigma_theta          = sigma_theta,
    weight_cv_empirical  = wt_result$cv_empirical,
    weight_deff_empirical = wt_result$deff_empirical,
    J                    = J,
    n_j                  = n_j,
    N                    = N,
    icc                  = icc,
    beta                 = beta,
    cv_w                 = cv_w,
    informative          = informative,
    prevalence           = mean(y)
  )
}


# =============================================================================
# Section 2: Survey Weight Generation
# =============================================================================
#
# Survey weights control the design effect (DEFF), which is the key driver
# of DER. We generate lognormal weights with a target CV, producing
# DEFF ~ 1 + CV^2 (the Kish approximation, Section 2.1, Eq. 3).
#
# Two regimes:
#
# (A) Non-informative weights (Eq. 17):
#     log(w_ij) ~ N(0, sigma^2_w), independent of (y, theta).
#     The weights affect efficiency but not consistency.
#     sigma^2_w = log(1 + CV_w^2) from the lognormal CV formula.
#
# (B) Informative weights (Eq. 18):
#     log(w_ij) = alpha * theta_j + epsilon_ij, epsilon_ij ~ N(0, sigma^2_eps)
#     The weights depend on the random effects, creating confounding:
#     units in clusters with large theta_j (high outcome probability)
#     have systematically different selection probabilities.
#
#     We calibrate alpha and sigma_eps so that:
#       - Total log-variance = alpha^2 * sigma^2_theta + sigma^2_eps
#                            = log(1 + CV_w^2)  [matching target CV]
#       - Signal fraction = alpha^2 * sigma^2_theta / total_var = 0.30
#         [producing moderate informativeness, Cor(theta, mean_j(log w)) ~ 0.3]
#
# After generation, weights are normalized within each cluster so that
# sum(w_ij) = n_j per cluster. This is consistent with pseudo-likelihood
# scaling conventions (Pfeffermann et al., 1998).
# =============================================================================

#' Generate survey weights (non-informative or informative)
#'
#' @param J Integer. Number of clusters.
#' @param n_j Integer. Cluster size.
#' @param cv_w Numeric > 0. Target CV of weights.
#' @param informative Logical. Whether weights depend on random effects.
#' @param theta_j Numeric vector of length J. True random effects.
#' @param seed Integer or NULL. Random seed (usually NULL).
#'
#' @return A list with:
#'   weights: Numeric vector of length N (normalized within clusters),
#'   cv_empirical: Empirical CV of raw (pre-normalization) weights,
#'   deff_empirical: Empirical Kish DEFF = 1 + CV^2.
generate_survey_weights <- function(J, n_j, cv_w, informative,
                                    theta_j, seed = NULL) {

  stopifnot(
    is.numeric(J), length(J) == 1, J >= 2,
    is.numeric(n_j), length(n_j) == 1, n_j >= 2,
    is.numeric(cv_w), length(cv_w) == 1, cv_w > 0,
    is.logical(informative), length(informative) == 1,
    is.numeric(theta_j), length(theta_j) == J
  )

  if (!is.null(seed)) set.seed(seed)

  N <- as.integer(J) * as.integer(n_j)
  group_idx <- rep(seq_len(J), each = n_j)

  # Total log-variance needed for the target CV.
  # For lognormal X with log(X) ~ N(mu, sigma^2):
  #   CV(X) = sqrt(exp(sigma^2) - 1)
  # Inverting: sigma^2 = log(1 + CV^2)
  log_var_total <- log(1 + cv_w^2)

  if (!informative) {
    # -- Non-informative weights (Eq. 17) ------------------------------------
    # log(w_ij) ~ N(0, sigma^2_w), independent of theta
    sigma_w <- sqrt(log_var_total)
    log_w   <- rnorm(N, mean = 0, sd = sigma_w)

  } else {
    # -- Informative weights (Eq. 18) ----------------------------------------
    # Decompose total log-variance into signal (from theta) and noise:
    #   Var(log w) = alpha^2 * Var(theta) + Var(epsilon)
    #              = alpha^2 * sigma^2_theta + sigma^2_eps
    #              = log_var_total
    #
    # Strategy: set signal_fraction = 0.30
    #   => alpha^2 * sigma^2_theta = 0.30 * log_var_total
    #   => sigma^2_eps = 0.70 * log_var_total
    #
    # For n_j = 50, averaging over epsilon within clusters yields:
    #   Cor(theta_j, mean_j(log w)) ~ alpha * sd(theta) /
    #     sqrt(alpha^2 * sigma^2_theta + sigma^2_eps / n_j)
    # With signal_fraction = 0.30 and n_j = 50, this correlation is
    # approximately 0.3-0.5, representing moderate informativeness.

    sigma2_theta_emp <- var(theta_j)
    # Guard against degenerate theta (extremely unlikely but safe)
    if (sigma2_theta_emp < 1e-10) {
      sigma2_theta_emp <- icc_to_sigma2(0.05)
    }

    signal_fraction <- 0.30
    alpha2     <- signal_fraction * log_var_total / sigma2_theta_emp
    alpha      <- sqrt(alpha2)
    sigma_eps2 <- (1 - signal_fraction) * log_var_total
    sigma_eps  <- sqrt(sigma_eps2)

    epsilon <- rnorm(N, mean = 0, sd = sigma_eps)
    log_w   <- alpha * theta_j[group_idx] + epsilon
  }

  # -- Exponentiate to get raw weights ----------------------------------------
  w_raw <- exp(log_w)

  # -- Compute empirical diagnostics before normalization ---------------------
  cv_empirical   <- sd(w_raw) / mean(w_raw)
  deff_empirical <- 1 + cv_empirical^2

  # -- Normalize within each cluster ------------------------------------------
  # Rescale so that sum(w_ij) = n_j for each cluster j.
  # This is the standard pseudo-likelihood normalization convention.
  w_normalized <- normalize_weights_by_group(w_raw, group_idx, n_j)

  list(
    weights        = w_normalized,
    cv_empirical   = cv_empirical,
    deff_empirical = deff_empirical
  )
}


# =============================================================================
# Section 3: Weight Normalization
# =============================================================================

#' Normalize weights within each cluster
#'
#' Rescales weights so that sum(w_ij) = n_j for each cluster j.
#' This normalization is standard in pseudo-likelihood estimation and
#' ensures that each cluster contributes proportionally to its sample size.
#'
#' @param w Numeric vector of length N. Raw weights.
#' @param group Integer vector of length N. Cluster membership (1 to J).
#' @param n_j Integer. Target sum per cluster.
#' @return Numeric vector of length N. Normalized weights.
normalize_weights_by_group <- function(w, group, n_j) {
  stopifnot(length(w) == length(group))

  w_norm <- numeric(length(w))
  groups <- sort(unique(group))

  for (j in groups) {
    idx   <- which(group == j)
    w_sum <- sum(w[idx])
    if (w_sum < .Machine$double.eps) {
      stop("Cluster ", j, " has zero total weight. ",
           "This indicates a degenerate weight generation.")
    }
    w_norm[idx] <- w[idx] * n_j / w_sum
  }

  return(w_norm)
}


# =============================================================================
# Section 4: Diagnostic Helpers
# =============================================================================

#' Compute empirical diagnostics for a generated dataset
#'
#' Checks that the generated data match the intended design: weight CV,
#' DEFF, prevalence, random-effect variability, and informativeness
#' correlation.
#'
#' @param data_list List returned by generate_survey_data().
#' @return A named list of diagnostic quantities.
diagnose_generated_data <- function(data_list) {
  stopifnot(is.list(data_list), "y" %in% names(data_list))

  w   <- data_list$weights
  y   <- data_list$y
  grp <- data_list$group

  # Weight diagnostics
  wt_cv   <- sd(w) / mean(w)
  wt_deff <- 1 + wt_cv^2

  # Outcome diagnostics
  prevalence    <- mean(y)
  prev_by_group <- tapply(y, grp, mean)

  # Random effect recovery
  theta_sd_empirical <- sd(data_list$theta_true)

  # Informativeness check: Cor(theta_j, mean_j(w))
  # If weights are non-informative, this should be near zero.
  # If informative, this should be substantially positive.
  group_mean_w <- tapply(w, grp, mean)
  cor_theta_w  <- cor(data_list$theta_true, group_mean_w)

  list(
    weight_cv       = wt_cv,
    weight_deff     = wt_deff,
    weight_min      = min(w),
    weight_max      = max(w),
    weight_ratio    = max(w) / min(w),
    prevalence      = prevalence,
    prevalence_range = range(prev_by_group),
    theta_sd        = theta_sd_empirical,
    cor_theta_w     = cor_theta_w,
    N               = data_list$N,
    J               = data_list$J
  )
}


#' Print a formatted diagnostics report
#'
#' @param data_list List returned by generate_survey_data().
print_data_diagnostics <- function(data_list) {
  dx <- diagnose_generated_data(data_list)

  cat("----------------------------------------------------------------\n")
  cat("Generated Data Diagnostics\n")
  cat("----------------------------------------------------------------\n")
  cat(sprintf("N = %d, J = %d, n_j = %d\n",
              dx$N, dx$J, data_list$n_j))
  cat(sprintf("Prevalence:         %.3f\n", dx$prevalence))
  cat(sprintf("Prevalence range:   [%.3f, %.3f]\n",
              dx$prevalence_range[1], dx$prevalence_range[2]))
  cat(sprintf("theta SD:           %.4f (target: %.4f)\n",
              dx$theta_sd, data_list$sigma_theta))
  cat(sprintf("Weight CV:          %.4f (target: %.4f)\n",
              dx$weight_cv, data_list$cv_w))
  cat(sprintf("Weight DEFF:        %.4f (target: %.4f)\n",
              dx$weight_deff, 1 + data_list$cv_w^2))
  cat(sprintf("Weight range:       [%.4f, %.4f] (ratio: %.2f)\n",
              dx$weight_min, dx$weight_max, dx$weight_ratio))
  cat(sprintf("Cor(theta, mean w): %.4f (%s)\n",
              dx$cor_theta_w,
              ifelse(data_list$informative, "informative", "non-informative")))
  cat("----------------------------------------------------------------\n")
}


# =============================================================================
# Section 5: Self-Validation
# =============================================================================

if (interactive()) {
  # Source configuration if not already loaded
  if (!exists("icc_to_sigma2")) {
    source(file.path(dirname(sys.frame(1)$ofile), "sim_00_config.R"))
  }

  cat("\n=== Testing generate_survey_data ===\n")

  # Test 1: Non-informative weights, small cluster count
  cat("\nTest 1: Non-informative, J=20, cv_w=0.3, icc=0.05\n")
  d1 <- generate_survey_data(
    J = 20, n_j = 50, icc = 0.05,
    beta = c(-0.5, 0.5, 0.3),
    cv_w = 0.3, informative = FALSE, seed = 42
  )
  print_data_diagnostics(d1)

  stopifnot(length(d1$y) == 1000)
  stopifnot(nrow(d1$X) == 1000, ncol(d1$X) == 3)
  stopifnot(all(d1$y %in% c(0L, 1L)))
  stopifnot(all(d1$weights > 0))

  # Verify within-cluster weight normalization: sum(w) = n_j per cluster
  wt_sums <- tapply(d1$weights, d1$group, sum)
  stopifnot(all(abs(wt_sums - 50) < 1e-10))

  # Verify group-mean centering of x_within: mean = 0 per cluster
  x_means <- tapply(d1$X[, "x_within"], d1$group, mean)
  stopifnot(all(abs(x_means) < 1e-10))

  # Test 2: Informative weights, high CV and ICC
  cat("\nTest 2: Informative, J=50, cv_w=1.5, icc=0.30\n")
  d2 <- generate_survey_data(
    J = 50, n_j = 50, icc = 0.30,
    beta = c(-0.5, 0.5, 0.3),
    cv_w = 1.5, informative = TRUE, seed = 123
  )
  print_data_diagnostics(d2)
  stopifnot(length(d2$y) == 2500)

  # Test 3: Reproducibility
  cat("\nTest 3: Reproducibility check\n")
  d3a <- generate_survey_data(
    J = 20, n_j = 50, icc = 0.15,
    beta = c(-0.5, 0.5, 0.3),
    cv_w = 1.0, informative = FALSE, seed = 999
  )
  d3b <- generate_survey_data(
    J = 20, n_j = 50, icc = 0.15,
    beta = c(-0.5, 0.5, 0.3),
    cv_w = 1.0, informative = FALSE, seed = 999
  )
  stopifnot(identical(d3a$y, d3b$y))
  stopifnot(identical(d3a$weights, d3b$weights))
  cat("  Reproducibility confirmed.\n")

  cat("\nAll self-validation checks passed.\n")
}
