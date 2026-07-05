# =============================================================================
# sim_04_run_single.R: Single-Replication Pipeline (dual-target DER)
# =============================================================================
#
# Purpose : Execute ONE (cell, rep) of the Monte Carlo study end to end. This
#           is the per-replication stage that run_batch drives in parallel:
#
#   1. Generate data with the two-stage DGP (generate_survey_data,
#      sim_01_dgp.R), seeded by get_rep_seed(cell_number, rep_id).
#   2. Fit hlr_weighted.stan via the sim_02 machinery (fit_hlr). The
#      PARALLEL_CHAINS environment variable is respected; if unset it
#      defaults to 1 (the worker-safe setting).
#   3. Run the STRICT convergence gate on ALL parameters and record the
#      full structured verdict. Gate-failing replications are NOT
#      discarded: they are persisted with gate$pass = FALSE and the
#      aggregation step decides how to treat them.
#   4. Compute the DUAL-TARGET sandwich/DER with the helpers
#      (canonical meat: center = TRUE, df_correct = TRUE):
#        - target_group: cluster = model group, single stratum;
#        - target_psu  : cluster = design PSU, strata = model group
#          (the PSU stage is stratified by group) -- only in the
#          psu_within_groups arm. In the sampled_groups arm the design
#          PSU IS the model group, so target_psu coincides with
#          target_group and is not recomputed (targets_coincide = TRUE).
#   5. Classify parameters at tau (classify_parameters) and apply, per
#      target, three interval strategies:
#        - naive     : uncorrected posterior draws;
#        - selective : block-Cholesky correction of the flagged set
#          (helpers selective_correct, method = "block_cholesky");
#        - blanket   : marginal per-parameter rescale of ALL p + J
#          parameters (see .blanket_marginal for why a full-matrix
#          Cholesky is not used at the group target).
#   6. Evaluate 90% equal-tail CI coverage and widths against the known
#      truth: beta_true (length 3) and truth$theta (length J).
#   7. Persist ONE compact .rds per rep at
#      out_dir/<cell_id>/rep_<0000>.rds (target <= ~100 KB; no draws, no
#      matrices larger than (p+J) x (p+J)), with skip-if-exists
#      resumability.
#
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models.
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Inputs  : code/helpers/sandwich_functions.R   (bread/meat/sandwich)
#           code/helpers/der_functions.R        (DER, classify, selective)
#           code/01_simulation/sim_00_config.R (SIM_PARAMS, grid, seeds)
#           code/01_simulation/sim_01_dgp.R  (two-stage DGP)
#           code/01_simulation/sim_02_fit.R     (fit_hlr, strict gate)
#           code/01_simulation/sim_03_postprocess.R (CI/coverage utilities)
# Called by: sim_05_run_batch.R
# =============================================================================


# =============================================================================
# Section 0: Package-root detection and dependency bootstrap
# =============================================================================

#' Locate the replication package root
#'
#' Checks, in order: the SVYDER_REPLICATION_ROOT environment variable, the
#' working directory and up to six of its parents, and here::here(). A
#' directory qualifies if it contains both the helpers and the config.
.find_root <- function() {
  probe <- function(d) {
    nzchar(d) &&
      file.exists(file.path(d, "code", "helpers", "sandwich_functions.R")) &&
      file.exists(file.path(d, "code", "01_simulation", "sim_00_config.R"))
  }
  env_root <- Sys.getenv("SVYDER_REPLICATION_ROOT", "")
  if (probe(env_root)) return(normalizePath(env_root))
  d <- getwd()
  for (i in 1:6) {
    if (probe(d)) return(normalizePath(d))
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  if (requireNamespace("here", quietly = TRUE)) {
    d <- tryCatch(here::here(), error = function(e) "")
    if (probe(d)) return(normalizePath(d))
  }
  stop("Cannot locate the svyder-replication package root. Run from the ",
       "package root or set SVYDER_REPLICATION_ROOT.", call. = FALSE)
}

# Source the pipeline dependencies once (idempotent: each file is skipped
# if its sentinel function already exists in the session). chdir = TRUE so
# sim_03_postprocess.R's own helper locator resolves ../helpers correctly.
local({
  root <- .find_root()
  src <- function(rel) source(file.path(root, rel), chdir = TRUE)
  if (!exists("build_J_cluster"))            src("code/helpers/sandwich_functions.R")
  if (!exists("selective_correct"))          src("code/helpers/der_functions.R")
  if (!exists("build_scenario_grid"))     src("code/01_simulation/sim_00_config.R")
  if (!exists("generate_survey_data"))    src("code/01_simulation/sim_01_dgp.R")
  if (!exists("fit_hlr"))                    src("code/01_simulation/sim_02_fit.R")
  if (!exists("compute_credible_intervals")) src("code/01_simulation/sim_03_postprocess.R")
})


# =============================================================================
# Section 1: Cached Stan model access
# =============================================================================

.v2_model_cache <- new.env(parent = emptyenv())

#' Get the compiled hlr_weighted CmdStanModel (memoized per session)
#'
#' cmdstanr itself caches the compiled binary next to the .stan file, so
#' repeated calls across sessions reuse the executable; this wrapper only
#' avoids re-constructing the R-side model object within a session.
#'
#' @param stan_path Path to the .stan file; NULL resolves
#'   SIM_PARAMS$stan_model_path against the detected package root.
#' @param quiet Suppress compilation messages.
#' @return A CmdStanModel object.
get_stan_model <- function(stan_path = NULL, quiet = TRUE) {
  if (is.null(stan_path)) {
    stan_path <- file.path(.find_root(), SIM_PARAMS$stan_model_path)
  }
  key <- normalizePath(stan_path, mustWork = TRUE)
  if (!is.null(.v2_model_cache$model) && identical(.v2_model_cache$key, key)) {
    return(.v2_model_cache$model)
  }
  model <- compile_stan_model(key, quiet = quiet)
  .v2_model_cache$model <- model
  .v2_model_cache$key   <- key
  model
}


# =============================================================================
# Section 2: Internal building blocks
# =============================================================================

#' Sandwich/DER for one aggregation target (internal)
#'
#' The bread H_obs and the posterior covariance sigma_mcmc do not depend on
#' the clustering, so they are computed once by the caller and shared; only
#' the meat differs between the model-group and design-PSU targets.
#'
#' @return list(V_sand, der, der_laplace, meat_rank_max = G_effective).
.target_der <- function(X, group, r, p, J, N, cluster, strata,
                           cluster_strata = NULL,
                           H_obs, sigma_mcmc, par_names,
                           meat = SIM_PARAMS$meat) {
  J_meat <- build_J_cluster(
    cluster_strata = cluster_strata,
    X = X, group = group, r = r, p = p, J = J, N = N,
    cluster = cluster, strata = strata,
    center = isTRUE(meat$center), df_correct = isTRUE(meat$df_correct)
  )
  sw <- compute_sandwich(H_obs, J_meat)
  dr <- compute_der(sw$V_sand, sigma_mcmc, H_obs_inv = sw$H_obs_inv,
                    param_names = par_names)
  list(
    V_sand      = sw$V_sand,
    der         = dr$der,
    der_laplace = dr$der_laplace,
    n_clusters  = if (!is.null(cluster_strata)) length(cluster_strata)
                  else max(cluster)
  )
}


#' Blanket correction: marginal per-parameter rescale (internal)
#'
#' Rescales EVERY parameter k by sqrt(diag(V_sand)[k] / diag(sigma_mcmc)[k])
#' around the posterior mean. A full-matrix (joint Cholesky) blanket is
#' generally INFEASIBLE at the model-group target: the centered meat is a
#' sum of G = J cluster outer products, so rank(J_meat) <= J - 1 < d = 3 + J,
#' V_sand is rank-deficient, and chol(V_sand) fails (a nearPD fallback would
#' silently fabricate the missing directions). The marginal operator matches
#' each parameter's sandwich variance without requiring a full-rank joint
#' target; the rep file records blanket_method = "marginal_per_parameter".
.blanket_marginal <- function(draws, V_sand, sigma_mcmc,
                                 point_est = colMeans(draws)) {
  diag_V <- pmax(diag(V_sand), 0)     # PSD numerics guard
  diag_m <- diag(sigma_mcmc)
  stopifnot(all(diag_m > 0), length(point_est) == ncol(draws))
  s <- sqrt(diag_V / diag_m)
  centered <- sweep(draws, 2L, point_est)
  scaled   <- sweep(centered, 2L, s, `*`)
  sweep(scaled, 2L, point_est, `+`)
}


#' Equal-tail CI coverage and width evaluation (internal)
#'
#' Thin adapter over compute_credible_intervals() / compute_coverage()
#' (sim_03_postprocess.R), returning named numeric vectors ready for
#' compact persistence.
.ci_eval <- function(draws, truth, level) {
  ci      <- compute_credible_intervals(draws, level = level)
  covered <- compute_coverage(truth, ci)
  list(
    cover = setNames(as.numeric(covered), ci$param),
    width = setNames(ci$width, ci$param)
  )
}


# =============================================================================
# Section 3: Single-replication pipeline
# =============================================================================

#' Run one replication end to end and persist a compact rep file
#'
#' @param cell One row of build_scenario_grid() (single-row data.frame),
#'   or a cell_id string (e.g., "J020_N10_GMOD_SG").
#' @param rep_id Integer replication index (1..R).
#' @param out_dir Base output directory; the rep file is written to
#'   out_dir/<cell_id>/rep_<0000>.rds.
#' @param mcmc MCMC settings list (default SIM_PARAMS$mcmc).
#' @param quiet Suppress Stan output (default TRUE).
#' @param stan_model Pre-built CmdStanModel or NULL (compile/load via
#'   get_stan_model(); parallel drivers pass a worker-local model
#'   re-created from the master's executable path).
#' @param grid Scenario grid or NULL (rebuilt deterministically). The
#'   cell_number used for seeding is the ROW INDEX of the cell in this grid.
#' @param tau DER flagging threshold (default SIM_PARAMS$tau).
#' @param ci_level Credible level (default SIM_PARAMS$ci_level = 0.90).
#' @param beta_prior_sd Prior SD on beta in hlr_weighted.stan (5).
#' @param overwrite If FALSE (default), an existing rep file is returned
#'   as-is with $skipped = TRUE (resumability).
#'
#' @return The persisted result list (invisibly). Hard failures (data
#'   generation or Stan sampling errors) raise conditions; the batch
#'   driver catches them and writes rep_<id>_ERROR.txt.
run_single_rep <- function(cell, rep_id, out_dir,
                              mcmc = SIM_PARAMS$mcmc, quiet = TRUE,
                              stan_model = NULL, grid = NULL,
                              tau = SIM_PARAMS$tau,
                              ci_level = SIM_PARAMS$ci_level,
                              beta_prior_sd = 5,
                              overwrite = FALSE) {

  t0     <- proc.time()["elapsed"]
  timing <- list()

  # -- Resolve cell, cell_number, and seed ------------------------------------
  if (is.null(grid)) grid <- build_scenario_grid()
  if (is.character(cell)) {
    stopifnot(length(cell) == 1L, cell %in% grid$cell_id)
    cell <- grid[grid$cell_id == cell, , drop = FALSE]
  }
  stopifnot(is.data.frame(cell), nrow(cell) == 1L)
  cell_id     <- cell$cell_id
  cell_number <- match(cell_id, grid$cell_id)
  stopifnot(!is.na(cell_number))
  rep_id <- as.integer(rep_id)
  seed   <- get_rep_seed(cell_number, rep_id)

  # -- Skip-if-exists (resumability) ------------------------------------------
  cell_dir <- file.path(out_dir, cell_id)
  out_file <- file.path(cell_dir, sprintf("rep_%04d.rds", rep_id))
  if (!overwrite && file.exists(out_file)) {
    res <- readRDS(out_file)
    res$skipped <- TRUE
    return(invisible(res))
  }

  # -- PARALLEL_CHAINS: respected if set, default 1 (worker-safe) -------------
  if (!nzchar(Sys.getenv("PARALLEL_CHAINS"))) {
    Sys.setenv(PARALLEL_CHAINS = "1")
    on.exit(Sys.unsetenv("PARALLEL_CHAINS"), add = TRUE)
  }
  pc_used <- as.integer(Sys.getenv("PARALLEL_CHAINS"))

  if (is.null(stan_model)) stan_model <- get_stan_model(quiet = quiet)

  # -- Step 1: Generate data (real two-stage sampling) -------------------------
  t1   <- proc.time()["elapsed"]
  data <- generate_survey_data(cell, seed)
  timing$data_gen <- as.numeric(proc.time()["elapsed"] - t1)

  p <- ncol(data$X)
  J <- as.integer(cell$J)
  N <- data$diagnostics$N
  stopifnot(
    p == 3L,
    max(data$group) == J,
    length(data$truth$theta) == J,
    length(data$truth$beta)  == p
  )

  # -- Step 2: Fit the Stan model (sim_02 machinery) ---------------------------
  fit_data <- list(y = data$y, X = data$X, group = data$group,
                   weights = data$weights, J = J, N = N)
  fit_res <- fit_hlr(fit_data, stan_model, mcmc_settings = mcmc,
                     seed = seed + 1000000L, quiet = quiet,
                     gate_config = SIM_PARAMS$convergence_gate)
  timing$fit <- fit_res$timing
  if (!fit_res$succeeded) {
    stop("Stan fitting failed for ", cell_id, " rep ", rep_id, ": ",
         fit_res$error, call. = FALSE)
  }

  # -- Step 3: Strict gate verdict (recorded, never used to discard) ----------
  gate_verdict <- fit_res$convergence
  gate_verdict$pass <- isTRUE(gate_verdict$passed)   # convenience alias

  # -- Step 4: Posterior draws, order-checked ----------------------------------
  par_names <- c(sprintf("beta[%d]", seq_len(p)),
                 sprintf("theta[%d]", seq_len(J)))
  draws_raw <- fit_res$fit$draws(variables = par_names, format = "matrix")
  stopifnot(identical(colnames(draws_raw), par_names))   # posterior pkg order
  draws <- matrix(as.numeric(draws_raw), nrow(draws_raw), ncol(draws_raw),
                  dimnames = list(NULL, par_names))

  sigma_theta_hat <- mean(as.numeric(
    fit_res$fit$draws(variables = "sigma_theta", format = "matrix")
  ))
  stopifnot(sigma_theta_hat > 0)
  beta_hat  <- colMeans(draws)[seq_len(p)]
  theta_hat <- colMeans(draws)[p + seq_len(J)]

  # -- Steps 5-7: sandwich, DER, corrections, coverage (warnings collected) ----
  t5 <- proc.time()["elapsed"]
  warn_log <- character(0)

  core <- withCallingHandlers({

    # Shared bread and posterior covariance (target-independent)
    eta  <- as.numeric(data$X %*% beta_hat) + theta_hat[data$group]
    q    <- plogis(eta)
    wtv  <- q * (1 - q)
    rres <- data$weights * (data$y - q)

    H_obs <- build_H_obs_logistic(data$X, data$group, data$weights, wtv,
                                  p, J, N, sigma_theta_hat, beta_prior_sd)
    sigma_mcmc <- cov(draws)

    # -- Dual targets ----------------------------------------------------------
    target_group <- .target_der(
      X = data$X, group = data$group, r = rres, p = p, J = J, N = N,
      cluster = data$group, strata = NULL,
      H_obs = H_obs, sigma_mcmc = sigma_mcmc, par_names = par_names
    )
    timing_sandwich_group <- as.numeric(proc.time()["elapsed"] - t5)

    targets_coincide <- identical(cell$structure, "sampled_groups")
    if (targets_coincide) {
      # Design PSU IS the model group: target_psu = target_group by
      # construction; do not recompute.
      stopifnot(identical(as.integer(data$psu), as.integer(data$group)))
      target_psu <- NULL
    } else {
      # PSU stage stratified by group: strata = data$strata (= group).
      stopifnot(!is.null(data$strata),
                identical(as.integer(data$strata), as.integer(data$group)))
      stopifnot(!is.null(data$psu_strata_map))
      target_psu <- .target_der(
        X = data$X, group = data$group, r = rres, p = p, J = J, N = N,
        cluster = data$psu, strata = data$strata,
        cluster_strata = data$psu_strata_map,
        H_obs = H_obs, sigma_mcmc = sigma_mcmc, par_names = par_names
      )
    }

    # -- Classification at tau (paper taxonomy) --------------------------------
    param_types <- c("fe_between", "fe_within", "fe_between", rep("re", J))
    cls_group <- classify_parameters(target_group$der, param_types, tau = tau)
    cls_psu   <- if (!is.null(target_psu)) {
      classify_parameters(target_psu$der, param_types, tau = tau)
    } else NULL

    # -- Truth vector for coverage ---------------------------------------------
    truth <- setNames(c(data$truth$beta, data$truth$theta), par_names)

    # -- Corrections and interval evaluation per target -------------------------
    naive_eval <- .ci_eval(draws, truth, ci_level)

    eval_target <- function(target, cls) {
      sel <- selective_correct(draws, target$der, target$V_sand, sigma_mcmc,
                               tau = tau, method = "block_cholesky")
      # The flagged set is DER > tau under both codepaths; keep them synced
      # (compare values only: selective_correct inherits names from der,
      # classify_parameters does not).
      stopifnot(identical(unname(as.integer(sel$corrected_idx)),
                          unname(as.integer(which(cls$flagged)))))
      blk_draws <- .blanket_marginal(draws, target$V_sand, sigma_mcmc)
      list(
        flagged   = setNames(sel$corrected_idx,
                             par_names[sel$corrected_idx]),
        n_flagged = sel$n_corrected,
        eval = list(
          naive     = naive_eval,
          selective = .ci_eval(sel$draws_corrected, truth, ci_level),
          blanket   = .ci_eval(blk_draws, truth, ci_level)
        )
      )
    }

    res_group <- eval_target(target_group, cls_group)
    res_psu   <- if (!is.null(target_psu)) eval_target(target_psu, cls_psu) else NULL

    # -- Small realized-design extras (used by the aggregation) ----------------
    deff_j <- tryCatch(compute_group_deff(data$group, data$weights, J),
                       error = function(e) NULL)
    B_j <- tryCatch(
      compute_shrinkage_factors(data$group, data$weights, wtv,
                                sigma_theta_hat, J),
      error = function(e) NULL
    )
    theory_group <- if (!is.null(deff_j) && !is.null(B_j)) {
      tryCatch(compute_theoretical_der(B_j, deff_j, J, sigma_theta_hat),
               error = function(e) NULL)
    } else NULL

    list(target_group = target_group, target_psu = target_psu,
         targets_coincide = targets_coincide,
         res_group = res_group, res_psu = res_psu,
         param_types = param_types, deff_j = deff_j, B_j = B_j,
         theory_group = theory_group,
         timing_sandwich_group = timing_sandwich_group)

  }, warning = function(w) {
    warn_log <<- c(warn_log, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  timing$postprocess <- as.numeric(proc.time()["elapsed"] - t5)
  timing$parallel_chains_used <- pc_used
  timing$total <- as.numeric(proc.time()["elapsed"] - t0)

  # -- Step 8: Assemble the compact per-rep record -----------------------------
  # Everything is O(p + J) except nothing: no draws, no d x d matrices.
  result <- list(
    pipeline    = "sim_04_run_single",
    cell_id     = cell_id,
    cell_number = cell_number,
    rep_id      = rep_id,
    seed        = seed,
    cell        = list(J = J, nbar_j = cell$nbar_j,
                       informativeness = cell$informativeness,
                       structure = cell$structure),
    gate        = gate_verdict,           # full structured strict verdict
    gate_pass   = gate_verdict$pass,
    diagnostics = data$diagnostics,       # realized N, n_j, DEFF, cv_w, ...

    # Dual-target DER
    targets_coincide = core$targets_coincide,
    der_group        = core$target_group$der,
    der_laplace_group = core$target_group$der_laplace,
    der_psu          = if (!core$targets_coincide) core$target_psu$der else NULL,
    der_laplace_psu  = if (!core$targets_coincide) core$target_psu$der_laplace else NULL,
    n_clusters       = c(group = core$target_group$n_clusters,
                         psu   = if (!core$targets_coincide)
                                   core$target_psu$n_clusters else NA_integer_),

    # Flagged sets and correction bookkeeping
    tau            = tau,
    param_types    = core$param_types,
    flagged        = list(group = core$res_group$flagged,
                          psu   = if (!core$targets_coincide)
                                    core$res_psu$flagged else NULL),
    n_flagged      = c(group = core$res_group$n_flagged,
                       psu   = if (!core$targets_coincide)
                                 core$res_psu$n_flagged else NA_integer_),
    blanket_method = "marginal_per_parameter",

    # Coverage indicators (0/1) and CI widths, per (target x correction)
    ci_level = ci_level,
    coverage = list(
      group = list(naive     = core$res_group$eval$naive$cover,
                   selective = core$res_group$eval$selective$cover,
                   blanket   = core$res_group$eval$blanket$cover),
      psu   = if (!core$targets_coincide) {
        list(naive     = core$res_psu$eval$naive$cover,
             selective = core$res_psu$eval$selective$cover,
             blanket   = core$res_psu$eval$blanket$cover)
      } else NULL
    ),
    ci_width = list(
      group = list(naive     = core$res_group$eval$naive$width,
                   selective = core$res_group$eval$selective$width,
                   blanket   = core$res_group$eval$blanket$width),
      psu   = if (!core$targets_coincide) {
        list(naive     = core$res_psu$eval$naive$width,
             selective = core$res_psu$eval$selective$width,
             blanket   = core$res_psu$eval$blanket$width)
      } else NULL
    ),

    # Point estimates
    sigma_theta_hat = sigma_theta_hat,
    beta_hat        = beta_hat,
    theta_hat       = theta_hat,

    # Realized per-group design quantities + closed-form predictions
    deff_j       = core$deff_j,
    B_j          = core$B_j,
    theory_group = core$theory_group,

    # Provenance
    meat     = SIM_PARAMS$meat,
    mcmc     = mcmc,
    warnings = if (length(warn_log) > 0) warn_log else NULL,
    timing   = timing,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  # -- Step 9: Persist atomically (write tmp, then rename) ---------------------
  dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
  tmp_file <- paste0(out_file, ".tmp")
  saveRDS(result, tmp_file)
  file.rename(tmp_file, out_file)

  invisible(result)
}
