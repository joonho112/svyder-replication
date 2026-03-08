# =============================================================================
# app_02_model_fit.R: Hierarchical Logistic Regression via cmdstanr
# =============================================================================
#
# Purpose : Compile and fit the weighted hierarchical logistic regression
#           (hlr_weighted.stan) to the NSECE 2019 data prepared in Step 1.
#           Run convergence diagnostics (Go/No-Go gate), extract posterior
#           draws, and save model outputs for DER analysis.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Section 5 (Application: NSECE 2019)
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
# Track   : A/B
# Inputs  :
#   Track A: data/precomputed/application/stan_data.rds (from app_01)
#            stan/hlr_weighted.stan
#   Track B: data/precomputed/application/model_summary.rds
# Outputs :
#   data/precomputed/application/model_fit.rds     (CmdStanMCMC object)
#   data/precomputed/application/model_summary.rds  (draws + diagnostics)
#
# Model specification (Equation 16 in the paper):
#   y_i | beta, theta_{j[i]} ~ Bernoulli(logit^{-1}(X_i beta + theta_{j[i]}))
#   theta_j ~ N(0, sigma_theta^2)
#
# The model uses:
#   - Non-centered parameterization (NCP) for random effects
#   - Pseudo-likelihood weighting with normalized survey weights
#   - 4 chains x 2000 warmup + 2000 sampling iterations
#   - adapt_delta = 0.95, max_treedepth = 14
#
# Convergence criteria (Go/No-Go gate):
#   - Divergent transitions: 0
#   - Rhat (max): < 1.01
#   - ESS bulk (min): > 400
#   - ESS tail (min): > 200
# =============================================================================


# =============================================================================
# Track selection
# =============================================================================

USE_PRECOMPUTED <- TRUE


# =============================================================================
# Setup
# =============================================================================

cat("==============================================================\n")
cat("  Application Step 2: Stan Model Fitting\n")
cat("==============================================================\n\n")

## Detect project root
if (requireNamespace("here", quietly = TRUE)) {
  PROJECT_ROOT <- here::here()
} else {
  PROJECT_ROOT <- getwd()
}
cat(sprintf("[Setup] Project root: %s\n\n", PROJECT_ROOT))

## Paths
PRECOMP_DIR  <- file.path(PROJECT_ROOT, "data", "precomputed", "application")
STAN_DIR     <- file.path(PROJECT_ROOT, "stan")


# =============================================================================
# Track B: Load pre-computed results
# =============================================================================

if (USE_PRECOMPUTED) {

  cat("--- Track B: Loading pre-computed model results ---\n\n")

  summary_path <- file.path(PRECOMP_DIR, "model_summary.rds")

  if (!file.exists(summary_path)) {
    stop(
      "Pre-computed model summary not found:\n  ", summary_path,
      "\n\nTo use Track B, ensure this file is present.",
      "\nSet USE_PRECOMPUTED <- FALSE for Track A (requires Stan)."
    )
  }

  model_summary <- readRDS(summary_path)

  ## Validate essential fields
  expected_fields <- c("summary", "draws_beta", "draws_sigma",
                       "draws_theta", "diagnostics", "mcmc_config")
  missing_fields  <- setdiff(expected_fields, names(model_summary))
  if (length(missing_fields) > 0) {
    stop("model_summary.rds missing fields: ",
         paste(missing_fields, collapse = ", "))
  }

  ## Print summary
  cat(sprintf("  draws_beta:  %d x %d\n",
              nrow(model_summary$draws_beta), ncol(model_summary$draws_beta)))
  cat(sprintf("  draws_theta: %d x %d\n",
              nrow(model_summary$draws_theta), ncol(model_summary$draws_theta)))
  cat(sprintf("  draws_sigma: length %d\n", length(model_summary$draws_sigma)))
  cat(sprintf("  Diagnostics:\n"))
  cat(sprintf("    Divergent transitions: %d\n",
              model_summary$diagnostics$n_divergent))
  cat(sprintf("    Rhat max: %.4f\n",
              model_summary$diagnostics$rhat_max))
  cat(sprintf("    ESS bulk min: %.0f\n",
              model_summary$diagnostics$ess_bulk_min))

  ## Print posterior summary for key parameters
  summ_df <- model_summary$summary
  for (par in c("beta[1]", "beta[2]", "beta[3]", "sigma_theta")) {
    row <- summ_df[summ_df$variable == par, , drop = FALSE]
    if (nrow(row) > 0) {
      cat(sprintf("    %-14s mean = %6.3f, sd = %.3f\n",
                  par, row$mean, row$sd))
    }
  }

  cat("\n  [Track B] Pre-computed model results loaded successfully.\n")
  cat("==============================================================\n")
  cat("  Step 2 Complete (Track B).\n")
  cat("==============================================================\n")

} else {

# =============================================================================
# Track A: Full model fitting pipeline
# =============================================================================

cat("--- Track A: Full Stan model fitting pipeline ---\n\n")
cat("  NOTE: This step requires cmdstanr and takes approximately\n")
cat("  20-60 minutes depending on hardware.\n\n")


# -----------------------------------------------------------------------------
# Section 1: Load Stan data
# -----------------------------------------------------------------------------

cat("--- 1. Loading Stan data ---\n")

stan_data_path <- file.path(PRECOMP_DIR, "stan_data.rds")
if (!file.exists(stan_data_path)) {
  stop("Stan data not found:\n  ", stan_data_path,
       "\nRun app_01_data_prep.R first.")
}

stan_data <- readRDS(stan_data_path)

expected_fields <- c("N", "J", "p", "y", "X", "group", "w")
missing_fields  <- setdiff(expected_fields, names(stan_data))
if (length(missing_fields) > 0) {
  stop("Stan data missing fields: ", paste(missing_fields, collapse = ", "))
}

cat(sprintf("  N = %d, J = %d, p = %d\n",
            stan_data$N, stan_data$J, stan_data$p))
cat(sprintf("  Outcome prevalence: %.3f\n", mean(stan_data$y)))


# -----------------------------------------------------------------------------
# Section 2: Compile Stan model
# -----------------------------------------------------------------------------

cat("\n--- 2. Compiling Stan model ---\n")

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop(
    "Package 'cmdstanr' is required but not installed.\n",
    "Install with:\n",
    "  install.packages('cmdstanr',\n",
    "    repos = c('https://stan-dev.r-universe.dev',\n",
    "              'https://cloud.r-project.org'))\n",
    "Then run cmdstanr::install_cmdstan() to install CmdStan."
  )
}

library(cmdstanr)

stan_file <- file.path(STAN_DIR, "hlr_weighted.stan")
if (!file.exists(stan_file)) {
  stop("Stan model file not found: ", stan_file)
}

cat(sprintf("  Stan file: %s\n", stan_file))

stan_model <- cmdstan_model(
  stan_file       = stan_file,
  force_recompile = FALSE,
  quiet           = FALSE
)
cat("  [PASS] Compilation successful.\n")


# -----------------------------------------------------------------------------
# Section 3: MCMC configuration
# -----------------------------------------------------------------------------
# Configuration tuned for real survey data with moderate random-effect
# structure. Higher adapt_delta and max_treedepth than typical defaults
# to ensure reliable exploration of the posterior geometry.

mcmc_config <- list(
  chains        = 4L,        # standard number of chains
  warmup        = 2000L,     # extended warmup for real data
  sampling      = 2000L,     # 2000 post-warmup iterations per chain
  adapt_delta   = 0.95,      # elevated to prevent divergences
  max_treedepth = 14L,       # elevated for complex geometry
  seed          = 20260307L, # reproducibility seed
  refresh       = 500L       # progress reporting frequency
)


# -----------------------------------------------------------------------------
# Section 4: Fit model
# -----------------------------------------------------------------------------

cat("\n--- 3. Fitting Stan model ---\n")
cat(sprintf("  Configuration: %d chains x %d warmup + %d sampling\n",
            mcmc_config$chains, mcmc_config$warmup, mcmc_config$sampling))
cat(sprintf("  adapt_delta = %.2f, max_treedepth = %d\n",
            mcmc_config$adapt_delta, mcmc_config$max_treedepth))

parallel_chains <- min(mcmc_config$chains, parallel::detectCores())
cat(sprintf("  parallel_chains = %d (detected cores: %d)\n\n",
            parallel_chains, parallel::detectCores()))

t_start <- proc.time()["elapsed"]

fit <- stan_model$sample(
  data            = stan_data,
  chains          = mcmc_config$chains,
  parallel_chains = parallel_chains,
  iter_warmup     = mcmc_config$warmup,
  iter_sampling   = mcmc_config$sampling,
  adapt_delta     = mcmc_config$adapt_delta,
  max_treedepth   = mcmc_config$max_treedepth,
  seed            = mcmc_config$seed,
  refresh         = mcmc_config$refresh,
  show_messages   = TRUE,
  show_exceptions = TRUE
)

t_elapsed <- as.numeric(proc.time()["elapsed"] - t_start)
cat(sprintf("\n  Sampling completed in %.1f seconds (%.1f minutes).\n",
            t_elapsed, t_elapsed / 60))


# -----------------------------------------------------------------------------
# Section 5: Convergence diagnostics (Go/No-Go gate)
# -----------------------------------------------------------------------------
# The Go/No-Go gate enforces strict convergence criteria before proceeding
# to DER analysis. Failed convergence would invalidate the sandwich variance
# calculations that underpin the entire DER framework.

cat("\n--- 4. Convergence diagnostics ---\n")

## Extract diagnostic quantities
diag_summ       <- fit$diagnostic_summary(quiet = TRUE)
n_divergent     <- sum(diag_summ$num_divergent)
n_max_treedepth <- sum(diag_summ$num_max_treedepth)

## Summary for key parameters (exclude score_residual for speed)
key_pars <- c("beta", "sigma_theta", "theta")
summ_key <- fit$summary(variables = key_pars, .cores = 1L)

rhat_max     <- max(summ_key$rhat, na.rm = TRUE)
ess_bulk_min <- min(summ_key$ess_bulk, na.rm = TRUE)
ess_tail_min <- min(summ_key$ess_tail, na.rm = TRUE)

diagnostics <- list(
  n_divergent       = as.integer(n_divergent),
  n_max_treedepth   = as.integer(n_max_treedepth),
  rhat_max          = rhat_max,
  ess_bulk_min      = ess_bulk_min,
  ess_tail_min      = ess_tail_min,
  rhat_by_param     = setNames(summ_key$rhat, summ_key$variable),
  ess_bulk_by_param = setNames(summ_key$ess_bulk, summ_key$variable),
  ess_tail_by_param = setNames(summ_key$ess_tail, summ_key$variable)
)

## Report
pass_div      <- n_divergent  == 0L
pass_rhat     <- rhat_max     <  1.01
pass_ess_bulk <- ess_bulk_min >  400
pass_ess_tail <- ess_tail_min >  200

cat(sprintf("  Divergent transitions: %d  [%s]\n",
            n_divergent, ifelse(pass_div, "PASS", "FAIL")))
cat(sprintf("  Max treedepth hits:    %d  [%s]\n",
            n_max_treedepth,
            ifelse(n_max_treedepth == 0L, "PASS", "WARNING")))
cat(sprintf("  Rhat max:              %.4f  [%s]\n",
            rhat_max, ifelse(pass_rhat, "PASS", "FAIL")))
cat(sprintf("  ESS bulk min:          %.0f  [%s]\n",
            ess_bulk_min, ifelse(pass_ess_bulk, "PASS", "FAIL")))
cat(sprintf("  ESS tail min:          %.0f  [%s]\n",
            ess_tail_min, ifelse(pass_ess_tail, "PASS", "FAIL")))

## Enforce Go/No-Go gate
if (!pass_div)      stop("FAILED: divergent transitions > 0")
if (!pass_rhat)     stop("FAILED: Rhat max >= 1.01")
if (!pass_ess_bulk) stop("FAILED: ESS bulk min <= 400")
if (!pass_ess_tail) stop("FAILED: ESS tail min <= 200")

cat("  Go/No-Go: ALL CHECKS PASSED\n")


# -----------------------------------------------------------------------------
# Section 6: Extract posterior draws
# -----------------------------------------------------------------------------
# Draws are extracted in matrix format for downstream DER computation.
# The DER analysis requires the full joint posterior covariance, which
# is estimated from the MCMC draws as Sigma_MCMC = Cov(phi^(1), ..., phi^(M))
# where phi = (beta, theta) is the (p + J)-dimensional parameter vector
# (Definition 1, denominator).

cat("\n--- 5. Extracting posterior draws ---\n")

## Fixed effects: matrix [draws x p]
draws_beta <- fit$draws(variables = "beta", format = "matrix")
stopifnot(ncol(draws_beta) == stan_data$p)

## sigma_theta: vector [draws]
draws_sigma_mat <- fit$draws(variables = "sigma_theta", format = "matrix")
draws_sigma     <- as.numeric(draws_sigma_mat[, 1L])

## Random effects: matrix [draws x J]
draws_theta <- fit$draws(variables = "theta", format = "matrix")
stopifnot(ncol(draws_theta) == stan_data$J)

cat(sprintf("  draws_beta:  %d x %d\n", nrow(draws_beta), ncol(draws_beta)))
cat(sprintf("  draws_sigma: length %d\n", length(draws_sigma)))
cat(sprintf("  draws_theta: %d x %d\n", nrow(draws_theta), ncol(draws_theta)))


# -----------------------------------------------------------------------------
# Section 7: Posterior summary
# -----------------------------------------------------------------------------

cat("\n--- 6. Posterior summary ---\n")

## Print key parameter summaries
for (par in c("beta[1]", "beta[2]", "beta[3]", "sigma_theta")) {
  row <- summ_key[summ_key$variable == par, , drop = FALSE]
  if (nrow(row) > 0) {
    q5  <- if ("q5"  %in% names(row)) row$q5  else NA
    q95 <- if ("q95" %in% names(row)) row$q95 else NA
    cat(sprintf("  %-14s  mean = %6.3f, sd = %.3f, 90%% CI [%.3f, %.3f]\n",
                par, row$mean, row$sd, q5, q95))
  }
}

## Random effect summary
theta_means <- colMeans(draws_theta)
cat(sprintf("\n  Random effects (theta):\n"))
cat(sprintf("    Mean of means: %.4f\n", mean(theta_means)))
cat(sprintf("    SD of means:   %.4f\n", sd(theta_means)))
cat(sprintf("    Sum of means:  %.4f (should be ~0)\n", sum(theta_means)))


# -----------------------------------------------------------------------------
# Section 8: Save outputs
# -----------------------------------------------------------------------------

cat("\n--- 7. Saving outputs ---\n")

## Save CmdStanMCMC fit object
fit_path <- file.path(PRECOMP_DIR, "model_fit.rds")
fit$save_object(file = fit_path)
cat(sprintf("  Saved fit: %s (%.1f MB)\n",
            fit_path, file.info(fit_path)$size / 1024^2))

## Assemble and save summary list
model_summary <- list(
  summary     = summ_key,          # data.frame with variable/mean/sd/...
  diagnostics = diagnostics,        # convergence diagnostics
  draws_beta  = draws_beta,         # matrix [M x p]
  draws_sigma = draws_sigma,        # vector [M]
  draws_theta = draws_theta,        # matrix [M x J]
  mcmc_config = mcmc_config,        # MCMC configuration
  timing      = t_elapsed           # elapsed seconds
)

summary_path <- file.path(PRECOMP_DIR, "model_summary.rds")
saveRDS(model_summary, summary_path)
cat(sprintf("  Saved summary: %s (%.1f MB)\n",
            summary_path, file.info(summary_path)$size / 1024^2))


# =============================================================================
# Final report
# =============================================================================

cat(sprintf("\n  Total elapsed time: %.1f seconds (%.1f minutes)\n",
            t_elapsed, t_elapsed / 60))
cat("\n==============================================================\n")
cat("  Step 2 Complete (Track A).\n")
cat("==============================================================\n")

}  # end Track A
