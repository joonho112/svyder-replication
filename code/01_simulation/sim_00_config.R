# =============================================================================
# sim_00_config.R: Simulation Design Configuration
# =============================================================================
#
# Purpose : Define the 54-scenario factorial design for the Monte Carlo
#           simulation study. This script sets up the complete experimental
#           grid, global parameters, and helper functions that all subsequent
#           simulation scripts depend upon.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Section 4 (Simulation Study), Table 2 (Factorial Design)
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Track   : A (Full Replication)
# Inputs  : None (self-contained configuration)
# Outputs : Objects loaded into the R environment:
#             - SIM_PARAMS        : Global simulation parameters (list)
#             - build_scenario_grid() : Function returning 54-row data.frame
#             - icc_to_sigma2()   : ICC-to-variance conversion
#             - make_scenario_id(): Scenario naming utility
#             - get_rep_seed()    : Deterministic seed generator
#
# Usage   : source("code/01_simulation/sim_00_config.R")
# =============================================================================


# =============================================================================
# Section 1: Global Simulation Parameters
# =============================================================================
#
# These parameters are shared across all 54 scenarios. They define the fixed
# aspects of the simulation that do not vary across the factorial grid.
#
# Rationale for specific choices:
#   - n_j = 50 : Moderate cluster size, typical of many survey designs.
#     Large enough for within-cluster estimation stability, small enough
#     for shrinkage effects to be visible.
#   - beta_true = (-0.5, 0.5, 0.3) : Intercept giving ~38% prevalence,
#     moderate within-group effect, moderate between-group effect.
#     These are on the logit scale.
#   - R = 200 : Sufficient for coverage estimation with MC SE < 2%
#     (since SE of coverage ~ sqrt(0.9 * 0.1 / 200) ~ 0.021).
# =============================================================================

SIM_PARAMS <- list(

  # -- Cluster structure (balanced design) ------------------------------------
  # All clusters have the same size n_j. This simplification allows clean
  # separation of weight-driven design effects from size-driven effects.
  n_j = 50L,

  # -- True regression coefficients -------------------------------------------
  # The hierarchical logistic model is:
  #   logit(p_ij) = beta_0 + beta_1 * x_wc_ij + beta_2 * z_bc_j + theta_j
  #
  # Parameter roles (Section 4.1):
  #   beta_0 = -0.5 : intercept (between-group parameter)
  #   beta_1 =  0.5 : within-cluster covariate effect (group-mean centered)
  #   beta_2 =  0.3 : between-cluster covariate effect
  #
  # beta_0 and beta_2 are "between-group" parameters because they are
  # confounded with the cluster-level random effect theta_j. beta_1 is a
  # "within-group" parameter because group-mean centering removes the
  # between-cluster component.
  beta_true = c(-0.5, 0.5, 0.3),

  # -- Monte Carlo replications -----------------------------------------------
  R       = 200L,                   # main Monte Carlo replications
  R_quick = 50L,                    # reduced count for debugging/testing

  # -- MCMC settings for cmdstanr ---------------------------------------------
  # These settings balance computational cost against inferential quality.
  #   - 4 chains : Standard for convergence diagnostics (Rhat, ESS).
  #   - 1000 warmup : Sufficient for adaptation in most scenarios.
  #   - 1500 sampling : Yields 6000 total post-warmup draws (4 * 1500).
  #   - adapt_delta = 0.90 : Slightly above default (0.80) to reduce
  #     divergent transitions in high-ICC scenarios.
  #   - max_treedepth = 12 : Above default (10) for robustness.
  mcmc = list(
    chains        = 4L,
    warmup        = 1000L,
    sampling      = 1500L,
    adapt_delta   = 0.90,
    max_treedepth = 12L
  ),

  # -- Stan model path (relative to replication package root) -----------------
  stan_model_path = "stan/hlr_weighted.stan",

  # -- DER classification threshold -------------------------------------------
  # Parameters with DER > tau are classified as "survey-sensitive" and
  # receive the Cholesky-based variance correction. See Section 3.3.
  tau = 1.2,

  # -- Credible interval level ------------------------------------------------
  ci_level = 0.90,

  # -- Seed management --------------------------------------------------------
  # Seeds are deterministic: seed = base_seed + scenario_number * 10000 + rep_id
  # This ensures full reproducibility and non-overlapping RNG streams.
  base_seed = 20260308L
)


# =============================================================================
# Section 2: ICC-to-Variance Conversion
# =============================================================================
#
# For a logistic model, the level-1 variance on the latent scale is fixed at
#   Var(logistic) = pi^2 / 3 ~ 3.29
# by the logistic distribution. The intraclass correlation is therefore:
#
#   ICC = sigma^2_theta / (sigma^2_theta + pi^2/3)         ... (Eq. 7)
#
# Inverting this gives the random-effect variance:
#
#   sigma^2_theta = ICC * pi^2 / (3 * (1 - ICC))           ... (Eq. 8)
#
# Reference: Goldstein, Browne & Rasbash (2002, JRSS-A) for the latent
# variable interpretation of ICC in binary response models.
# =============================================================================

#' Convert latent-scale ICC to random-effect variance
#'
#' @param icc Numeric scalar in (0, 1). Target intraclass correlation.
#' @return Numeric scalar. sigma^2_theta = ICC * pi^2 / (3 * (1 - ICC)).
icc_to_sigma2 <- function(icc) {
  stopifnot(is.numeric(icc), all(icc > 0), all(icc < 1))
  icc * pi^2 / (3 * (1 - icc))
}

#' Convert latent-scale ICC to random-effect standard deviation
#'
#' @param icc Numeric scalar in (0, 1).
#' @return Numeric scalar. sigma_theta = sqrt(sigma^2_theta).
icc_to_sigma <- function(icc) {
  sqrt(icc_to_sigma2(icc))
}


# =============================================================================
# Section 3: Factorial Design Grid
# =============================================================================
#
# The simulation crosses four factors (Section 4.1, Table 2):
#
#   Factor 1 -- Number of clusters J in {20, 50, 100}
#     Rationale: Spans small (J=20, typical of many state-level surveys)
#     to moderate (J=100) cluster counts. DER behavior depends on J through
#     the shrinkage factor B_j (Theorem 2).
#
#   Factor 2 -- Weight CV (cv_w) in {0.3, 1.0, 2.0}
#     Rationale: Controls the Kish design effect via DEFF ~ 1 + CV_w^2
#     (Section 2.1). CV = 0.3 gives mild weighting (DEFF ~ 1.09),
#     CV = 1.0 gives moderate (DEFF ~ 2.0), CV = 2.0 gives severe
#     (DEFF ~ 5.0). These span the range observed in practice
#     (Kish, 1965; Potthoff et al., 1992).
#
#   Factor 3 -- Latent ICC in {0.05, 0.15, 0.30}
#     Rationale: Determines the balance between within- and between-cluster
#     variance on the latent scale. Low ICC (0.05) means weak clustering;
#     high ICC (0.30) means strong clustering. The ICC affects DER through
#     the shrinkage factor: stronger clustering => more shrinkage => lower
#     DER for random effects (Theorem 2).
#
#   Factor 4 -- Weight informativeness in {non-informative, informative}
#     Rationale: Tests whether the DER framework correctly identifies
#     additional bias from informative sampling. Under informative weights,
#     log(w_ij) depends on the random effect theta_j, creating confounding.
#
# Total: 3 x 3 x 3 x 2 = 54 scenarios.
# Each scenario is replicated R = 200 times.
# =============================================================================

#' Build the 54-scenario factorial grid
#'
#' @return A data.frame with 54 rows and columns:
#'   scenario_id, J, n_j, N, cv_w, icc, informative,
#'   sigma2_theta, sigma_theta, approx_deff
build_scenario_grid <- function() {

  factors <- expand.grid(
    J           = c(20L, 50L, 100L),
    cv_w        = c(0.3, 1.0, 2.0),
    icc         = c(0.05, 0.15, 0.30),
    informative = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  # -- Derived quantities -----------------------------------------------------

  # Cluster size is fixed across scenarios
  factors$n_j <- SIM_PARAMS$n_j

  # Total sample size: N = J * n_j
  factors$N <- factors$J * factors$n_j

  # Random-effect variance from ICC (Section 4.1, Eq. 8)
  factors$sigma2_theta <- vapply(
    factors$icc, icc_to_sigma2, numeric(1)
  )
  factors$sigma_theta <- sqrt(factors$sigma2_theta)

  # Approximate Kish design effect: DEFF ~ 1 + CV_w^2 (Section 2.1, Eq. 3)
  # This is the theoretical DEFF under lognormal weights with the given CV.
  factors$approx_deff <- 1 + factors$cv_w^2

  # -- Scenario identifiers ---------------------------------------------------
  factors$scenario_id <- make_scenario_id(
    J           = factors$J,
    cv_w        = factors$cv_w,
    icc         = factors$icc,
    informative = factors$informative
  )

  # -- Column ordering for readability ----------------------------------------
  col_order <- c(
    "scenario_id", "J", "n_j", "N", "cv_w", "icc",
    "informative", "sigma2_theta", "sigma_theta", "approx_deff"
  )
  factors <- factors[, col_order]

  # Sort by scenario_id for deterministic ordering
  factors <- factors[order(factors$scenario_id), ]
  rownames(factors) <- NULL

  return(factors)
}


# =============================================================================
# Section 4: Scenario Naming Utilities
# =============================================================================
#
# Each scenario receives a unique, human-readable identifier that encodes
# all four design factors. This facilitates file management, logging, and
# result retrieval.
#
# Format: "J020_CV030_ICC005_NI"
#   - J: cluster count, zero-padded to 3 digits
#   - CV: cv_w * 100, zero-padded to 3 digits
#   - ICC: icc * 100, zero-padded to 3 digits
#   - NI = non-informative weights, IN = informative weights
# =============================================================================

#' Create scenario ID string(s)
#'
#' @param J Integer vector. Number of clusters.
#' @param cv_w Numeric vector. Weight coefficient of variation.
#' @param icc Numeric vector. Latent-scale ICC.
#' @param informative Logical vector. Weight informativeness.
#' @return Character vector of scenario IDs.
make_scenario_id <- function(J, cv_w, icc, informative) {
  stopifnot(
    length(J) == length(cv_w),
    length(J) == length(icc),
    length(J) == length(informative)
  )

  j_str   <- sprintf("J%03d", J)
  cv_str  <- sprintf("CV%03d", round(cv_w * 100))
  icc_str <- sprintf("ICC%03d", round(icc * 100))
  inf_str <- ifelse(informative, "IN", "NI")

  paste(j_str, cv_str, icc_str, inf_str, sep = "_")
}


#' Parse a scenario ID back to its design factors
#'
#' @param scenario_id Character scalar. E.g., "J020_CV030_ICC005_NI".
#' @return A named list with J, cv_w, icc, informative.
parse_scenario_id <- function(scenario_id) {
  stopifnot(is.character(scenario_id), length(scenario_id) == 1)

  parts <- strsplit(scenario_id, "_")[[1]]
  if (length(parts) != 4) {
    stop("Invalid scenario_id format: ", scenario_id,
         ". Expected format: J020_CV030_ICC005_NI")
  }

  J           <- as.integer(sub("^J", "", parts[1]))
  cv_w        <- as.numeric(sub("^CV", "", parts[2])) / 100
  icc         <- as.numeric(sub("^ICC", "", parts[3])) / 100
  informative <- switch(parts[4],
                        "NI" = FALSE,
                        "IN" = TRUE,
                        stop("Invalid informativeness code: ", parts[4]))

  list(J = J, cv_w = cv_w, icc = icc, informative = informative)
}


#' Look up a scenario row from the grid by its ID
#'
#' @param scenario_id Character scalar.
#' @param grid Data.frame (default: builds the grid fresh).
#' @return Single-row data.frame matching the scenario.
get_scenario <- function(scenario_id, grid = build_scenario_grid()) {
  idx <- which(grid$scenario_id == scenario_id)
  if (length(idx) == 0) {
    stop("Scenario not found: ", scenario_id)
  }
  grid[idx, , drop = FALSE]
}


# =============================================================================
# Section 5: Seed Management
# =============================================================================
# Deterministic seed assignment ensures full reproducibility.
# Each (scenario, replication) pair receives a unique seed:
#   seed = base_seed + scenario_number * 10000 + rep_id
#
# The scenario_number is the row index in the sorted scenario grid (1..54).
# This guarantees non-overlapping RNG streams across all runs.
# =============================================================================

#' Compute a deterministic seed for a given scenario and replication
#'
#' @param scenario_number Integer. Row index in the scenario grid (1..54).
#' @param rep_id Integer. Replication number (1..R).
#' @param base_seed Integer. Base seed (default: SIM_PARAMS$base_seed).
#' @return Integer. Unique seed for this (scenario, replication) pair.
get_rep_seed <- function(scenario_number, rep_id,
                         base_seed = SIM_PARAMS$base_seed) {
  as.integer(base_seed + scenario_number * 10000L + rep_id)
}


# =============================================================================
# Section 6: Grid Summary Display
# =============================================================================

#' Print a concise summary of the factorial design
#'
#' @param grid Data.frame from build_scenario_grid().
print_grid_summary <- function(grid = build_scenario_grid()) {
  cat("================================================================\n")
  cat("Simulation Study: Factorial Design Summary\n")
  cat("================================================================\n")
  cat(sprintf("Total scenarios:       %d\n", nrow(grid)))
  cat(sprintf("J values:              %s\n",
              paste(sort(unique(grid$J)), collapse = ", ")))
  cat(sprintf("CV_w values:           %s\n",
              paste(sort(unique(grid$cv_w)), collapse = ", ")))
  cat(sprintf("ICC values:            %s\n",
              paste(sort(unique(grid$icc)), collapse = ", ")))
  cat(sprintf("Informative weights:   %s\n",
              paste(sort(unique(grid$informative)), collapse = ", ")))
  cat(sprintf("Cluster size (n_j):    %d\n", SIM_PARAMS$n_j))
  cat(sprintf("Sample sizes (N):      %s\n",
              paste(sort(unique(grid$N)), collapse = ", ")))
  cat(sprintf("Replications:          %d (quick: %d)\n",
              SIM_PARAMS$R, SIM_PARAMS$R_quick))
  cat(sprintf("DER threshold (tau):   %.2f\n", SIM_PARAMS$tau))
  cat(sprintf("CI level:              %.0f%%\n", SIM_PARAMS$ci_level * 100))
  cat("----------------------------------------------------------------\n")
  cat("MCMC settings:\n")
  cat(sprintf("  chains = %d, warmup = %d, sampling = %d\n",
              SIM_PARAMS$mcmc$chains,
              SIM_PARAMS$mcmc$warmup,
              SIM_PARAMS$mcmc$sampling))
  cat(sprintf("  adapt_delta = %.2f, max_treedepth = %d\n",
              SIM_PARAMS$mcmc$adapt_delta,
              SIM_PARAMS$mcmc$max_treedepth))
  cat("----------------------------------------------------------------\n")
  cat("True coefficients (logit scale):\n")
  cat(sprintf("  beta = (%s)\n",
              paste(SIM_PARAMS$beta_true, collapse = ", ")))
  cat("----------------------------------------------------------------\n")
  cat("ICC -> random-effect SD mapping:\n")
  for (icc_val in sort(unique(grid$icc))) {
    cat(sprintf("  ICC = %.2f  =>  sigma^2_theta = %.4f,  sigma_theta = %.4f\n",
                icc_val, icc_to_sigma2(icc_val), icc_to_sigma(icc_val)))
  }
  cat("----------------------------------------------------------------\n")
  cat("Approximate DEFF from weight CV:\n")
  for (cv_val in sort(unique(grid$cv_w))) {
    cat(sprintf("  cv_w = %.1f  =>  DEFF approx %.2f\n",
                cv_val, 1 + cv_val^2))
  }
  cat("================================================================\n")
}


# =============================================================================
# Section 7: Self-Validation (runs when sourced interactively)
# =============================================================================

if (interactive()) {
  grid <- build_scenario_grid()
  print_grid_summary(grid)

  # Verify grid dimensions: 3 * 3 * 3 * 2 = 54
  stopifnot(nrow(grid) == 54)
  stopifnot(all(grid$n_j == 50))

  # Verify round-trip parsing of scenario IDs
  test_id <- "J050_CV100_ICC015_NI"
  parsed  <- parse_scenario_id(test_id)
  stopifnot(parsed$J == 50, parsed$cv_w == 1.0,
            parsed$icc == 0.15, parsed$informative == FALSE)

  # Verify ICC-to-variance mapping (analytic check)
  stopifnot(abs(icc_to_sigma2(0.05) - 0.05 * pi^2 / (3 * 0.95)) < 1e-10)

  # Verify seed uniqueness across all (scenario, rep) pairs
  all_seeds <- integer(0)
  for (sc in seq_len(nrow(grid))) {
    for (r in 1:5) {
      all_seeds <- c(all_seeds, get_rep_seed(sc, r))
    }
  }
  stopifnot(length(all_seeds) == length(unique(all_seeds)))

  cat("\nAll self-validation checks passed.\n")
}
