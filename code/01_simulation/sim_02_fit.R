# =============================================================================
# sim_02_fit.R: Stan Model Fitting and Convergence Diagnostics
# =============================================================================
#
# Purpose : Fit the weighted hierarchical logistic regression Stan model
#           (hlr_weighted.stan) to a single simulated dataset using cmdstanr.
#           Includes Stan data preparation, MCMC sampling, convergence
#           checking, and posterior draw extraction.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Section 4.2 (Model Fitting)
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Track   : A (Full Replication)
# Inputs  : stan/hlr_weighted.stan (Stan model file)
#           sim_00_config.R (for SIM_PARAMS, MCMC settings)
#           sim_01_dgp.R (for generate_survey_data output format)
# Outputs : fit_hlr() function returning a CmdStanMCMC fit object with
#           convergence diagnostics
#
# Stan model overview (hlr_weighted.stan):
#   The model implements a hierarchical logistic regression with:
#     - Non-centered parameterization (NCP) for random effects:
#         theta_j = sigma_theta * eta_raw_j, eta_raw_j ~ N(0,1)
#     - Pseudo-maximum likelihood via survey weight scaling:
#         target += w_i * bernoulli_logit_lpmf(y_i | linpred_i)
#     - Score residuals computed in generated quantities block:
#         score_residual_i = w_i * (y_i - q_i)
#       These are needed for sandwich variance estimation (sim_03).
#
#   Priors:
#     beta ~ N(0, 5)        [weakly informative on logit scale]
#     sigma_theta ~ N+(0, 2) [half-normal, allows moderate clustering]
# =============================================================================


# =============================================================================
# Section 1: Model Fitting
# =============================================================================

#' Fit the weighted hierarchical logistic regression via cmdstanr
#'
#' Prepares the Stan data list from the DGP output, runs MCMC sampling,
#' checks convergence, and returns a structured result. The function is
#' designed to be called within the single-replication pipeline
#' (sim_04_run.R).
#'
#' @param data_list List. Output from generate_survey_data().
#' @param stan_model A CmdStanModel object (pre-compiled).
#' @param mcmc_settings List with elements: chains, warmup, sampling,
#'   adapt_delta, max_treedepth. Defaults to SIM_PARAMS$mcmc.
#' @param seed Integer or NULL. MCMC seed for reproducibility.
#' @param quiet Logical. If TRUE, suppress Stan sampling output.
#'
#' @return A list with components:
#'   fit: CmdStanMCMC object (or NULL if fitting failed),
#'   convergence: list from check_convergence() (or NULL),
#'   timing: wall-clock seconds,
#'   error: error message or NULL,
#'   succeeded: logical.
fit_hlr <- function(data_list, stan_model, mcmc_settings = NULL,
                    seed = NULL, quiet = TRUE) {

  # -- Input validation -------------------------------------------------------
  stopifnot(is.list(data_list))
  required_fields <- c("y", "X", "group", "weights", "J", "N")
  missing_fields  <- setdiff(required_fields, names(data_list))
  if (length(missing_fields) > 0) {
    stop("data_list is missing required fields: ",
         paste(missing_fields, collapse = ", "))
  }

  if (!inherits(stan_model, "CmdStanModel")) {
    stop("stan_model must be a CmdStanModel object. ",
         "Compile it with cmdstanr::cmdstan_model().")
  }

  # Default MCMC settings from global configuration
  if (is.null(mcmc_settings)) {
    if (exists("SIM_PARAMS")) {
      mcmc_settings <- SIM_PARAMS$mcmc
    } else {
      mcmc_settings <- list(
        chains        = 4L,
        warmup        = 1000L,
        sampling      = 1500L,
        adapt_delta   = 0.90,
        max_treedepth = 12L
      )
    }
  }

  # -- Prepare Stan data ------------------------------------------------------
  stan_data <- prepare_stan_data(data_list)

  # -- Run MCMC sampling ------------------------------------------------------
  fit       <- NULL
  error_msg <- NULL
  t_start   <- proc.time()["elapsed"]

  tryCatch({
    fit <- stan_model$sample(
      data            = stan_data,
      chains          = mcmc_settings$chains,
      parallel_chains = as.integer(Sys.getenv(
        "PARALLEL_CHAINS",
        min(mcmc_settings$chains, parallel::detectCores())
      )),
      iter_warmup     = mcmc_settings$warmup,
      iter_sampling   = mcmc_settings$sampling,
      adapt_delta     = mcmc_settings$adapt_delta,
      max_treedepth   = mcmc_settings$max_treedepth,
      seed            = seed,
      refresh         = ifelse(quiet, 0, 200),
      show_messages   = !quiet,
      show_exceptions = !quiet
    )
  }, error = function(e) {
    error_msg <<- conditionMessage(e)
  })

  t_elapsed <- as.numeric(proc.time()["elapsed"] - t_start)

  # -- Check convergence ------------------------------------------------------
  convergence <- NULL
  if (!is.null(fit)) {
    convergence <- tryCatch(
      check_convergence(fit),
      error = function(e) {
        list(
          passed        = FALSE,
          rhat_max      = NA_real_,
          ess_bulk_min  = NA_real_,
          ess_tail_min  = NA_real_,
          n_divergent   = NA_integer_,
          pct_divergent = NA_real_,
          details       = paste("Convergence check failed:",
                                conditionMessage(e))
        )
      }
    )
  }

  list(
    fit         = fit,
    convergence = convergence,
    timing      = t_elapsed,
    error       = error_msg,
    succeeded   = !is.null(fit)
  )
}


# =============================================================================
# Section 2: Stan Data Preparation
# =============================================================================
#
# The Stan model (hlr_weighted.stan) expects a specific data format:
#   N     : total observations
#   J     : number of clusters
#   p     : number of fixed-effect predictors (= 3)
#   y     : integer array of binary outcomes
#   X     : N x p design matrix
#   group : integer array of cluster indicators (1..J)
#   w     : positive real vector of normalized survey weights
# =============================================================================

#' Prepare the Stan data list from DGP output
#'
#' Validates dimensions and types, then formats the data as required by
#' hlr_weighted.stan. The pseudo-likelihood approach scales each
#' observation's log-likelihood contribution by its survey weight w_i.
#'
#' @param data_list List from generate_survey_data().
#' @return Named list formatted for Stan.
prepare_stan_data <- function(data_list) {

  N <- data_list$N
  J <- data_list$J
  p <- ncol(data_list$X)

  # Dimension validation
  stopifnot(
    length(data_list$y) == N,
    nrow(data_list$X) == N,
    length(data_list$group) == N,
    length(data_list$weights) == N,
    all(data_list$y %in% c(0L, 1L)),
    all(data_list$group >= 1L & data_list$group <= J),
    all(data_list$weights > 0)
  )

  list(
    N     = N,
    J     = J,
    p     = p,
    y     = as.integer(data_list$y),
    X     = data_list$X,
    group = as.integer(data_list$group),
    w     = as.numeric(data_list$weights)
  )
}


# =============================================================================
# Section 3: Convergence Diagnostics
# =============================================================================
#
# We evaluate four standard MCMC diagnostics on the key parameters
# (beta and sigma_theta, not the J individual theta_j):
#
#   1. Rhat (split-R-hat, Vehtari et al. 2021): measures between-chain vs
#      within-chain variability. Rhat <= 1.01 is ideal; we use 1.10 as a
#      permissive threshold for the simulation.
#
#   2. Effective sample size (bulk): measures the number of independent draws.
#      ESS >= 100 per chain is a minimal requirement (Vehtari et al., 2021).
#
#   3. Effective sample size (tail): ESS in the tails of the posterior,
#      important for credible interval estimation.
#
#   4. Divergent transitions: indicate regions where the Hamiltonian dynamics
#      are unreliable. Zero is ideal; we tolerate up to 5%.
#
# Replications that fail convergence are excluded from coverage analysis
# to ensure valid inference.
# =============================================================================

#' Check MCMC convergence diagnostics
#'
#' @param fit CmdStanMCMC object.
#' @param rhat_threshold Numeric. Maximum acceptable Rhat (default: 1.10).
#' @param min_ess Integer. Minimum acceptable ESS (default: 100).
#' @param max_div_pct Numeric. Maximum acceptable divergent fraction
#'   (default: 0.05 = 5%).
#'
#' @return A list with: passed (logical), rhat_max, ess_bulk_min,
#'   ess_tail_min, n_divergent, pct_divergent, n_max_treedepth, details.
check_convergence <- function(fit, rhat_threshold = 1.10,
                              min_ess = 100, max_div_pct = 0.05) {

  if (!inherits(fit, "CmdStanMCMC")) {
    stop("fit must be a CmdStanMCMC object.")
  }

  # Focus on key parameters (beta, sigma_theta), not all theta[j]
  key_pars <- c("beta", "sigma_theta")
  summ     <- fit$summary(variables = key_pars)

  # Rhat
  rhat_vals <- summ$rhat
  rhat_max  <- max(rhat_vals, na.rm = TRUE)

  # Effective sample size (bulk and tail)
  ess_bulk_min <- min(summ$ess_bulk, na.rm = TRUE)
  ess_tail_min <- min(summ$ess_tail, na.rm = TRUE)

  # Divergent transitions
  diag_summary <- fit$diagnostic_summary(quiet = TRUE)
  n_divergent  <- sum(diag_summary$num_divergent)
  n_chains     <- length(diag_summary$num_divergent)

  # Total post-warmup iterations
  total_iter <- tryCatch({
    metadata <- fit$metadata()
    metadata$iter_sampling * n_chains
  }, error = function(e) {
    nrow(fit$draws(format = "matrix"))
  })

  pct_divergent   <- n_divergent / total_iter
  n_max_treedepth <- sum(diag_summary$num_max_treedepth)

  # -- Evaluate all criteria --------------------------------------------------
  issues <- character(0)

  rhat_ok <- rhat_max <= rhat_threshold
  if (!rhat_ok) {
    issues <- c(issues, sprintf("Rhat max = %.4f (threshold: %.2f)",
                                rhat_max, rhat_threshold))
  }

  ess_bulk_ok <- ess_bulk_min >= min_ess
  if (!ess_bulk_ok) {
    issues <- c(issues, sprintf("ESS bulk min = %.0f (threshold: %d)",
                                ess_bulk_min, min_ess))
  }

  ess_tail_ok <- ess_tail_min >= min_ess
  if (!ess_tail_ok) {
    issues <- c(issues, sprintf("ESS tail min = %.0f (threshold: %d)",
                                ess_tail_min, min_ess))
  }

  div_ok <- pct_divergent <= max_div_pct
  if (!div_ok) {
    issues <- c(issues,
                sprintf("Divergent transitions: %d / %d (%.1f%%)",
                        n_divergent, total_iter, pct_divergent * 100))
  }

  passed  <- rhat_ok && ess_bulk_ok && ess_tail_ok && div_ok
  details <- if (passed) {
    "All convergence diagnostics passed."
  } else {
    paste("CONVERGENCE ISSUES:", paste(issues, collapse = "; "))
  }

  list(
    passed          = passed,
    rhat_max        = rhat_max,
    ess_bulk_min    = ess_bulk_min,
    ess_tail_min    = ess_tail_min,
    n_divergent     = n_divergent,
    pct_divergent   = pct_divergent,
    n_max_treedepth = n_max_treedepth,
    details         = details
  )
}


# =============================================================================
# Section 4: Posterior Draw Extraction
# =============================================================================
#
# These helper functions extract quantities from the fitted model that are
# needed by the post-processing step (sim_03_postprocess.R):
#   - Parameter summaries (mean, quantiles, coverage)
#   - Full posterior draws as a matrix (for Cholesky correction)
#   - Score residuals (for sandwich variance estimation)
# =============================================================================

#' Extract posterior summaries for key parameters
#'
#' @param fit CmdStanMCMC object.
#' @param data_list List from generate_survey_data() (for true values).
#' @return A data.frame with parameter-level summaries and coverage indicators.
extract_parameter_summaries <- function(fit, data_list) {

  if (!inherits(fit, "CmdStanMCMC")) {
    stop("fit must be a CmdStanMCMC object.")
  }

  par_names  <- c("beta[1]", "beta[2]", "beta[3]", "sigma_theta")
  summ       <- fit$summary(variables = par_names)

  true_vals  <- c(data_list$beta, data_list$sigma_theta)
  par_labels <- c("beta_0 (intercept)", "beta_1 (x_within)",
                  "beta_2 (z_between)", "sigma_theta")

  result <- data.frame(
    parameter  = par_labels,
    true_value = true_vals,
    mean       = summ$mean,
    median     = summ$median,
    sd         = summ$sd,
    q025       = summ$q5,
    q975       = summ$q95,
    rhat       = summ$rhat,
    ess_bulk   = summ$ess_bulk,
    ess_tail   = summ$ess_tail,
    stringsAsFactors = FALSE
  )

  result$bias    <- result$mean - result$true_value
  result$covered <- (result$true_value >= result$q025) &
                    (result$true_value <= result$q975)

  return(result)
}


#' Extract posterior draws as a matrix for key parameters
#'
#' Returns a matrix with all post-warmup draws merged across chains.
#' Row count = chains * iter_sampling, column count = number of parameters.
#'
#' @param fit CmdStanMCMC object.
#' @param variables Character vector of parameter names to extract.
#' @return Numeric matrix [iterations x parameters].
extract_draws_matrix <- function(fit,
                                 variables = c("beta", "sigma_theta")) {
  if (!inherits(fit, "CmdStanMCMC")) {
    stop("fit must be a CmdStanMCMC object.")
  }
  fit$draws(variables = variables, format = "matrix")
}


#' Extract score residuals from the Stan model
#'
#' The score residuals s_i = w_i * (y_i - q_i) are computed in the
#' generated quantities block of hlr_weighted.stan. They are needed
#' for constructing the clustered meat matrix J_cluster in the
#' sandwich variance estimator. See Section 2.3 of Lee (2026).
#'
#' @param fit CmdStanMCMC object.
#' @return A matrix [iterations x N] of score residuals.
extract_score_residuals <- function(fit) {
  if (!inherits(fit, "CmdStanMCMC")) {
    stop("fit must be a CmdStanMCMC object.")
  }
  fit$draws(variables = "score_residual", format = "matrix")
}


# =============================================================================
# Section 5: Model Compilation Helper
# =============================================================================

#' Compile the Stan model from source
#'
#' @param stan_file Character. Path to the .stan file.
#' @param quiet Logical. Suppress compilation messages.
#' @return A CmdStanModel object.
compile_stan_model <- function(stan_file, quiet = FALSE) {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required. Install with:\n",
         "  install.packages('cmdstanr', repos = c(",
         "    'https://stan-dev.r-universe.dev',",
         "    'https://cloud.r-project.org'))")
  }

  if (!file.exists(stan_file)) {
    stop("Stan file not found: ", stan_file)
  }

  if (!quiet) {
    cat(sprintf("Compiling Stan model: %s\n", stan_file))
  }

  model <- cmdstanr::cmdstan_model(
    stan_file   = stan_file,
    cpp_options = list(stan_threads = FALSE),
    quiet       = quiet
  )

  if (!quiet) {
    cat("Compilation successful.\n")
  }

  return(model)
}


# =============================================================================
# Section 6: Printing Helpers
# =============================================================================

#' Print convergence diagnostics summary
#'
#' @param conv_result List from check_convergence().
print_convergence <- function(conv_result) {
  cat("----------------------------------------------------------------\n")
  cat("Convergence Diagnostics\n")
  cat("----------------------------------------------------------------\n")
  cat(sprintf("Status:            %s\n",
              ifelse(conv_result$passed, "PASSED", "*** FAILED ***")))
  cat(sprintf("Rhat max:          %.4f\n", conv_result$rhat_max))
  cat(sprintf("ESS bulk min:      %.0f\n", conv_result$ess_bulk_min))
  cat(sprintf("ESS tail min:      %.0f\n", conv_result$ess_tail_min))
  cat(sprintf("Divergences:       %d (%.2f%%)\n",
              conv_result$n_divergent, conv_result$pct_divergent * 100))
  cat(sprintf("Max treedepth:     %d\n", conv_result$n_max_treedepth))
  cat(sprintf("Details:           %s\n", conv_result$details))
  cat("----------------------------------------------------------------\n")
}


#' Print model fit summary
#'
#' @param fit_result List from fit_hlr().
print_fit_summary <- function(fit_result) {
  cat("================================================================\n")
  cat("Model Fit Summary\n")
  cat("================================================================\n")
  cat(sprintf("Succeeded:  %s\n", fit_result$succeeded))
  cat(sprintf("Timing:     %.1f seconds\n", fit_result$timing))
  if (!is.null(fit_result$error)) {
    cat(sprintf("Error:      %s\n", fit_result$error))
  }
  if (!is.null(fit_result$convergence)) {
    print_convergence(fit_result$convergence)
  }
  cat("================================================================\n")
}


# =============================================================================
# Section 7: Self-Validation (structural tests only -- no Stan fitting)
# =============================================================================

if (interactive()) {
  cat("\n=== Testing sim_02_fit.R helper functions ===\n")

  # Test prepare_stan_data with a small synthetic dataset
  if (exists("generate_survey_data")) {
    cat("\nTest: prepare_stan_data\n")
    test_data <- generate_survey_data(
      J = 10, n_j = 20, icc = 0.10,
      beta = c(-0.5, 0.5, 0.3),
      cv_w = 0.3, informative = FALSE, seed = 42
    )
    stan_data <- prepare_stan_data(test_data)

    stopifnot(stan_data$N == 200)
    stopifnot(stan_data$J == 10)
    stopifnot(stan_data$p == 3)
    stopifnot(length(stan_data$y) == 200)
    stopifnot(nrow(stan_data$X) == 200)
    stopifnot(ncol(stan_data$X) == 3)
    stopifnot(all(stan_data$group >= 1 & stan_data$group <= 10))
    stopifnot(all(stan_data$w > 0))
    cat("  prepare_stan_data: PASSED\n")
  } else {
    cat("  Skipping prepare_stan_data test.\n")
    cat("  Source sim_01_dgp.R first.\n")
  }

  cat("\nAll available self-validation checks passed.\n")
}
