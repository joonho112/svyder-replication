# =============================================================================
# sim_00_config.R: Simulation Design Configuration (24-cell grid)
# =============================================================================
#
# Purpose : Define the factorial design for the Monte Carlo study, which uses
#           a genuine two-stage sampling mechanism (finite population,
#           inclusion probabilities depending on the response/random effects,
#           exact inverse-probability weights). This is the first stage of the
#           pipeline (config -> dgp -> fit -> postprocess -> run_single ->
#           run_batch -> results); it defines SIM_PARAMS and the scenario grid
#           consumed by every downstream stage.
#
# Design  : 2 (J) x 2 (nbar_j) x 3 (informativeness gamma) x 2 (structure)
#           = 24 cells, R = 200 replications each.
#
#   Factor 1 -- J in {20, 50}: number of model groups.
#   Factor 2 -- nbar_j in {10, 50}: expected sampled units per group.
#     The nbar_j = 10 cells stress the small-cluster regime raised in the
#     JSSAM review (higher-order inclusion literature).
#   Factor 3 -- informativeness in {none, moderate, strong}: strength of the
#     dependence of inclusion probabilities on (theta, y). Realized DEFF is
#     MEASURED per replication, never asserted from a nominal CV.
#   Factor 4 -- structure in {"sampled_groups", "psu_within_groups"}:
#     * sampled_groups  (c = j): stage 1 samples J of M population groups
#       with PPS proportional to exp(gamma_grp * theta_g + noise); stage 2
#       samples units within groups by Poisson sampling with probabilities
#       proportional to exp(gamma_unit * y + noise). Design PSU = model
#       group.
#     * psu_within_groups (c != j): ALL J groups are observed (the
#       "all groups in the sample" regime; cf. the NSECE application where
#       all 51 states appear). Within each group, k of K population PSUs
#       are sampled with PPS proportional to exp(gamma_grp * upsilon_c +
#       noise), where upsilon_c is a PSU-level outcome shock that the
#       ANALYSIS MODEL DOES NOT INCLUDE (deliberate, realistic
#       misspecification: the fitted HLR has group effects only, as in the
#       NSECE analysis). Stage 2 is Poisson within PSUs. This arm powers
#       the dual-target (model-group vs design-PSU) analysis with known
#       truth; the PSU stage is stratified BY GROUP, so the stratified
#           meat (strata = group, C_h = k per stratum) gets real use.
#
# ICC is fixed at 0.15 across all cells; the design instead varies the
# sampling mechanism (informativeness and structural arm).
#
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
# =============================================================================

SIM_PARAMS <- list(

  # -- True regression coefficients -------------------------------------------
  #   logit(p) = beta_0 + beta_1 * x_wc + beta_2 * z_bc [+ theta_g (+ upsilon_c)]
  beta_true = c(-0.5, 0.5, 0.3),

  # -- Latent ICC (fixed) -----------------------------------------------------
  icc = 0.15,

  # -- PSU-level shock (structure = "psu_within_groups" only) -----------------
  # Latent-scale variance share of the PSU level. Small but real: this is
  # the design layer the analysis model omits, and the reason the
  # design-PSU target exceeds the model-group target.
  icc_psu = 0.05,

  # -- Population sizes -------------------------------------------------------
  # sampled_groups arm: M population groups; J are sampled (f1 = J/M).
  #   N_g population units per group; expected f2 = nbar_j / N_g.
  M_groups     = 200L,
  N_g_factor   = 20L,     # N_g = N_g_factor * nbar_j  (f2 = 5%)

  # psu_within_groups arm: K population PSUs per group, k sampled;
  #   N_c population units per PSU; stage-2 expected per-PSU sample
  #   = nbar_j / k.
  K_psu        = 12L,
  k_psu        = 4L,
  N_c_factor   = 10L,     # N_c = N_c_factor * (nbar_j / k_psu)  (f2 = 10%)

  # -- Informativeness calibration -------------------------------------------
  # Stage-1 size measure: s proportional to exp(gamma_grp * shock + e1),
  #   e1 ~ N(0, tau1^2); stage-2: s proportional to exp(gamma_unit * y + e2),
  #   e2 ~ N(0, tau2^2). The noise SDs give the none-arm its baseline weight
  #   variation (so DEFF > 1 even without informativeness); the gammas add
  #   the (theta, y)-dependence. The values below are calibrated so that the
  #   realized DEFF and the post-normalization informativeness signal
  #   separate the gamma levels; dgp_acceptance_test() (sim_01_dgp.R) checks
  #   these properties on the analysis weights.
  tau1 = 0.40,
  tau2 = 0.60,
  gamma_levels = list(
    none     = list(gamma_grp = 0.0, gamma_unit = 0.0),
    moderate = list(gamma_grp = 0.6, gamma_unit = 0.5),
    strong   = list(gamma_grp = 1.2, gamma_unit = 1.0)
  ),

  # -- Weight convention (declared; part of the variance target) --------------
  # Exact inverse-probability weights, then global unit-mean normalization.
  normalize_weights = "global_unit_mean",

  # -- Monte Carlo replications -----------------------------------------------
  R       = 200L,
  R_quick = 20L,

  # -- MCMC settings ----------------------------------------------------------
  mcmc = list(
    chains        = 4L,
    warmup        = 1000L,
    sampling      = 1500L,
    adapt_delta   = 0.90,
    max_treedepth = 12L
  ),

  stan_model_path = "stan/hlr_weighted.stan",

  # -- DER threshold / CI level -------------------------------------------------
  tau      = 1.2,
  ci_level = 0.90,

  # -- Meat options (canonical: centered, DF-corrected) -----------------------
  meat = list(center = TRUE, df_correct = TRUE),

  # -- Convergence gate: strict (see sim_02_fit.R for definitions) ------------
  convergence_gate = list(
    version = "strict",
    strict  = list(rhat_max = 1.01, params = "all",
                   max_divergences = 0L, min_ess = 100L),
    v1_weak = list(rhat_max = 1.10, params = c("beta", "sigma_theta"),
                   max_div_pct = 0.05, min_ess = 100L)
  ),

  # -- Seeds --------------------------------------------------------------------
  # seed = base_seed + cell_number * 10000 + rep_id (non-overlapping streams)
  base_seed = 20260704L
)

# Derived: sigma_theta from fixed ICC (logistic latent-scale conversion)
SIM_PARAMS$sigma2_theta <- SIM_PARAMS$icc * pi^2 / (3 * (1 - SIM_PARAMS$icc))
SIM_PARAMS$sigma_theta  <- sqrt(SIM_PARAMS$sigma2_theta)
# PSU shock variance from icc_psu (share of latent variance at PSU level)
SIM_PARAMS$sigma2_psu   <- SIM_PARAMS$icc_psu * pi^2 / (3 * (1 - SIM_PARAMS$icc - SIM_PARAMS$icc_psu))
SIM_PARAMS$sigma_psu    <- sqrt(SIM_PARAMS$sigma2_psu)


#' Build the 24-cell scenario grid
#'
#' @return data.frame with one row per cell:
#'   cell_id, J, nbar_j, informativeness, structure, gamma_grp, gamma_unit,
#'   sigma_theta, sigma_psu, R
build_scenario_grid <- function() {

  factors <- expand.grid(
    J               = c(20L, 50L),
    nbar_j          = c(10L, 50L),
    informativeness = c("none", "moderate", "strong"),
    structure       = c("sampled_groups", "psu_within_groups"),
    stringsAsFactors = FALSE
  )

  factors$gamma_grp <- vapply(
    factors$informativeness,
    function(g) SIM_PARAMS$gamma_levels[[g]]$gamma_grp, numeric(1)
  )
  factors$gamma_unit <- vapply(
    factors$informativeness,
    function(g) SIM_PARAMS$gamma_levels[[g]]$gamma_unit, numeric(1)
  )

  factors$sigma_theta <- SIM_PARAMS$sigma_theta
  factors$sigma_psu   <- ifelse(factors$structure == "psu_within_groups",
                                SIM_PARAMS$sigma_psu, 0)
  factors$R <- SIM_PARAMS$R

  factors$cell_id <- make_cell_id(
    J = factors$J, nbar_j = factors$nbar_j,
    informativeness = factors$informativeness,
    structure = factors$structure
  )

  col_order <- c("cell_id", "J", "nbar_j", "informativeness", "structure",
                 "gamma_grp", "gamma_unit", "sigma_theta", "sigma_psu", "R")
  factors <- factors[, col_order]
  factors <- factors[order(factors$cell_id), ]
  rownames(factors) <- NULL
  factors
}


#' Cell ID: "J020_N10_GMOD_SG" / "J050_N50_GSTR_PW" etc.
make_cell_id <- function(J, nbar_j, informativeness, structure) {
  g_str <- toupper(substr(informativeness, 1, 3))
  g_str <- c(NON = "GNON", MOD = "GMOD", STR = "GSTR")[substr(g_str, 1, 3)]
  s_str <- ifelse(structure == "sampled_groups", "SG", "PW")
  sprintf("J%03d_N%02d_%s_%s", J, nbar_j, g_str, s_str)
}


#' Parse a cell ID back to factors
parse_cell_id <- function(cell_id) {
  parts <- strsplit(cell_id, "_")[[1]]
  stopifnot(length(parts) == 4)
  J      <- as.integer(sub("^J", "", parts[1]))
  nbar_j <- as.integer(sub("^N", "", parts[2]))
  informativeness <- switch(parts[3],
                            GNON = "none", GMOD = "moderate", GSTR = "strong",
                            stop("bad gamma code: ", parts[3]))
  structure <- switch(parts[4],
                      SG = "sampled_groups", PW = "psu_within_groups",
                      stop("bad structure code: ", parts[4]))
  list(J = J, nbar_j = nbar_j, informativeness = informativeness,
       structure = structure)
}


#' Deterministic per-(cell, rep) seed
get_rep_seed <- function(cell_number, rep_id,
                            base_seed = SIM_PARAMS$base_seed) {
  as.integer(base_seed + cell_number * 10000L + rep_id)
}


# -- Self-validation ----------------------------------------------------------
if (interactive()) {
  g <- build_scenario_grid()
  stopifnot(nrow(g) == 24L,
            length(unique(g$cell_id)) == 24L,
            all(g$sigma_psu[g$structure == "sampled_groups"] == 0))
  rt <- parse_cell_id("J020_N10_GMOD_PW")
  stopifnot(rt$J == 20L, rt$nbar_j == 10L,
            rt$informativeness == "moderate",
            rt$structure == "psu_within_groups")
  cat("sim_00_config: all self-validation checks passed.\n")
}
