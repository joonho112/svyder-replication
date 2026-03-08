# =============================================================================
# sim_04_run.R: Parallel Execution of the Simulation Study
# =============================================================================
#
# Purpose : Execute the full Monte Carlo simulation study across all 54
#           scenarios and 200 replications per scenario. Manages parallel
#           execution, incremental saving, resume capability, and progress
#           reporting.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Section 4 (Simulation Study)
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Track   : A (Full Replication)
# Inputs  : sim_00_config.R  (scenario grid, global parameters)
#           sim_01_dgp.R     (data generation)
#           sim_02_fit.R     (Stan model fitting)
#           sim_03_postprocess.R (DER, correction, coverage)
#           stan/hlr_weighted.stan (Stan model)
# Outputs : Per-replication .rds files in data/simulation/{scenario_id}/
#           Per-scenario summary .rds files in data/simulation/summary/
#           Run log in data/simulation/run_log.txt
#
# Usage (from replication package root):
#   Rscript code/01_simulation/sim_04_run.R [options]
#
# Options:
#   --J, --j_levels   Comma-separated J values (default: "20,50,100")
#   --reps, -R        Number of replications (default: 200)
#   --cores, -c       Number of parallel cores (default: detected)
#   --outdir          Output directory (default: "data/simulation")
#   --resume          Skip completed scenario/rep combos
#   --scenarios       Comma-separated scenario IDs (overrides --J)
#   --quick           Use reduced replication count (R_quick) for testing
#
# Estimated runtime:
#   J = 20  (18 scenarios): ~4-8 hours on 8 cores
#   J = 50  (18 scenarios): ~10-20 hours on 8 cores
#   J = 100 (18 scenarios): ~24-48 hours on 8 cores
#   Full grid (54 scenarios): ~40-80 hours on 8 cores
#
# The simulation can be interrupted and resumed with --resume.
# Each completed replication is saved immediately.
# =============================================================================

cat("================================================================\n")
cat("  SIMULATION STUDY: Parallel Execution\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")


# =============================================================================
# Section 0: Command-Line Argument Parsing
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Defaults
  parsed <- list(
    j_levels  = c(20L, 50L, 100L),
    reps      = 200L,
    cores     = max(1L, parallel::detectCores() - 2L),
    outdir    = "data/simulation",
    resume    = FALSE,
    scenarios = NULL,
    quick     = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c("--J", "--j_levels")) {
      i <- i + 1
      parsed$j_levels <- as.integer(strsplit(args[i], ",")[[1]])
    } else if (arg %in% c("--reps", "-R")) {
      i <- i + 1
      parsed$reps <- as.integer(args[i])
    } else if (arg %in% c("--cores", "-c")) {
      i <- i + 1
      parsed$cores <- as.integer(args[i])
    } else if (arg == "--outdir") {
      i <- i + 1
      parsed$outdir <- args[i]
    } else if (arg == "--resume") {
      parsed$resume <- TRUE
    } else if (arg == "--scenarios") {
      i <- i + 1
      parsed$scenarios <- strsplit(args[i], ",")[[1]]
    } else if (arg == "--quick") {
      parsed$quick <- TRUE
    } else {
      warning("Unknown argument: ", arg)
    }
    i <- i + 1
  }

  return(parsed)
}

cli_args <- parse_args()

cat("--- Configuration ---\n")
if (!is.null(cli_args$scenarios)) {
  cat(sprintf("  Scenarios:    %s\n",
              paste(cli_args$scenarios, collapse = ", ")))
} else {
  cat(sprintf("  J levels:     %s\n",
              paste(cli_args$j_levels, collapse = ", ")))
}
cat(sprintf("  Replications: %d\n", cli_args$reps))
cat(sprintf("  Cores:        %d\n", cli_args$cores))
cat(sprintf("  Output dir:   %s\n", cli_args$outdir))
cat(sprintf("  Resume:       %s\n", cli_args$resume))
cat(sprintf("  Quick mode:   %s\n", cli_args$quick))
cat("\n")


# =============================================================================
# Section 1: Setup -- Paths and Dependencies
# =============================================================================

# Detect project root: use here::here() if available, else fall back
# to the PROJ_ROOT environment variable or the script's grandparent dir.
if (requireNamespace("here", quietly = TRUE)) {
  proj_root <- here::here()
} else {
  proj_root <- Sys.getenv(
    "PROJ_ROOT",
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))
  )
}

# Make output directory absolute
if (!startsWith(cli_args$outdir, "/")) {
  cli_args$outdir <- file.path(proj_root, cli_args$outdir)
}

# Create output directories
dir.create(cli_args$outdir, recursive = TRUE, showWarnings = FALSE)
summary_dir <- file.path(cli_args$outdir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# Section 2: Load Packages and Source Pipeline Scripts
# =============================================================================

cat("--- Loading packages ---\n")
suppressPackageStartupMessages({
  library(cmdstanr)
  library(Matrix)
  library(parallel)
})
cat("  Packages loaded.\n\n")

cat("--- Sourcing pipeline scripts ---\n")
script_dir <- file.path(proj_root, "code/01_simulation")
pipeline_scripts <- c(
  "sim_00_config.R",
  "sim_01_dgp.R",
  "sim_02_fit.R",
  "sim_03_postprocess.R"
)

for (script in pipeline_scripts) {
  cat(sprintf("  Sourcing: %s\n", script))
  source(file.path(script_dir, script))
}
cat("  All pipeline scripts sourced.\n\n")


# =============================================================================
# Section 3: Compile Stan Model
# =============================================================================
#
# We compile once in the main process and pass the executable path to
# worker processes. This avoids redundant compilation and ensures all
# workers use the same model binary.
#
# IMPORTANT: cmdstanr model objects cannot be serialized across forked
# processes (mclapply). Workers must re-create the CmdStanModel from the
# pre-compiled executable path.
# =============================================================================

cat("--- Compiling Stan model ---\n")
stan_file  <- file.path(proj_root, SIM_PARAMS$stan_model_path)
stan_model <- compile_stan_model(stan_file, quiet = FALSE)
exe_path   <- stan_model$exe_file()
cat(sprintf("  Compiled. Executable: %s\n\n", exe_path))


# =============================================================================
# Section 4: Build Scenario List
# =============================================================================

cat("--- Building scenario list ---\n")
grid <- build_scenario_grid()

if (!is.null(cli_args$scenarios)) {
  # Use explicitly specified scenario IDs
  scenario_ids <- cli_args$scenarios
  invalid <- setdiff(scenario_ids, grid$scenario_id)
  if (length(invalid) > 0) {
    stop("Invalid scenario IDs: ", paste(invalid, collapse = ", "))
  }
  scenarios_df <- grid[grid$scenario_id %in% scenario_ids, ]
  scenarios_df <- scenarios_df[match(scenario_ids,
                                     scenarios_df$scenario_id), ]
} else {
  # Filter by J levels
  scenarios_df <- grid[grid$J %in% cli_args$j_levels, ]
}

n_scenarios <- nrow(scenarios_df)
cat(sprintf("  Total scenarios: %d\n", n_scenarios))
cat(sprintf("  Scenario IDs: %s\n",
            paste(scenarios_df$scenario_id, collapse = ", ")))
cat("\n")

# -- Quick mode: override replication count for debugging/testing -----------
if (cli_args$quick) {
  cli_args$reps <- SIM_PARAMS$R_quick
  cat(sprintf("  ** QUICK MODE: Replications reduced to %d **\n\n",
              cli_args$reps))
}


# =============================================================================
# Section 5: Initialize Logging
# =============================================================================

log_file <- file.path(cli_args$outdir, "run_log.txt")

write_log <- function(msg, ...) {
  formatted <- sprintf(msg, ...)
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  line <- paste(timestamp, formatted)
  cat(line, "\n", file = log_file, append = TRUE)
  cat(line, "\n")
}

cat("", file = log_file)   # create/truncate log file
write_log("Simulation Study: Production Run")
write_log("  Scenarios: %d", n_scenarios)
write_log("  Replications: %d", cli_args$reps)
write_log("  Cores: %d", cli_args$cores)
write_log("  Resume: %s", cli_args$resume)
write_log("  Output: %s", cli_args$outdir)
write_log("---")


# =============================================================================
# Section 6: Resume Support
# =============================================================================

#' Identify already-completed replications for a scenario
#'
#' Used when --resume is active to skip replications that have already
#' been saved to disk.
#'
#' @param scenario_id Character. Scenario identifier.
#' @param outdir Character. Base output directory.
#' @return Integer vector of completed replication numbers.
get_existing_reps <- function(scenario_id, outdir) {
  scenario_dir <- file.path(outdir, scenario_id)
  if (!dir.exists(scenario_dir)) return(integer(0))

  rds_files <- list.files(scenario_dir, pattern = "^rep_\\d+\\.rds$")
  if (length(rds_files) == 0) return(integer(0))

  rep_nums <- as.integer(gsub("rep_(\\d+)\\.rds", "\\1", rds_files))
  sort(rep_nums)
}


# =============================================================================
# Section 7: Main Simulation Loop
# =============================================================================
#
# The outer loop iterates over scenarios sequentially. Within each scenario,
# replications are executed in parallel using parallel::mclapply.
#
# Sequential scenarios + parallel replications is the recommended strategy
# because:
#   1. Each scenario uses the same Stan model but different data sizes,
#      so memory usage is predictable per-scenario.
#   2. cmdstanr model objects cannot be shared across forked processes;
#      we pass the exe_path and re-create models in each worker.
#   3. Incremental saving per-replication allows resume capability.
# =============================================================================

t_global_start       <- proc.time()["elapsed"]
cumulative_reps      <- 0L
cumulative_converged <- 0L
cumulative_failed    <- 0L
scenario_timings     <- numeric(n_scenarios)

for (sc_idx in seq_len(n_scenarios)) {

  scenario <- scenarios_df[sc_idx, ]
  sc_id    <- scenario$scenario_id
  t_sc_start <- proc.time()["elapsed"]

  cat("\n")
  cat("================================================================\n")
  cat(sprintf("Scenario %d/%d: %s\n", sc_idx, n_scenarios, sc_id))
  cat(sprintf("  J=%d, cv_w=%.1f, icc=%.2f, informative=%s\n",
              scenario$J, scenario$cv_w, scenario$icc, scenario$informative))
  cat(sprintf("  N=%d, approx_deff=%.2f, sigma_theta=%.4f\n",
              scenario$N, scenario$approx_deff, scenario$sigma_theta))
  cat("================================================================\n")

  write_log("Starting scenario %d/%d: %s", sc_idx, n_scenarios, sc_id)

  # -- Determine which replications to run ------------------------------------
  reps_to_run <- seq_len(cli_args$reps)

  if (cli_args$resume) {
    existing_reps <- get_existing_reps(sc_id, cli_args$outdir)
    if (length(existing_reps) > 0) {
      reps_to_run <- setdiff(reps_to_run, existing_reps)
      cat(sprintf("  Resume: %d existing, %d remaining\n",
                  length(existing_reps), length(reps_to_run)))
      write_log("  Resume: %d existing, %d remaining",
                length(existing_reps), length(reps_to_run))
    }
  }

  if (length(reps_to_run) == 0) {
    cat("  All reps already completed. Skipping.\n")
    write_log("  Skipped (all reps complete)")
    scenario_timings[sc_idx] <- 0
    next
  }

  cat(sprintf("  Running %d reps on %d cores...\n",
              length(reps_to_run), cli_args$cores))

  # -- Run replications in parallel -------------------------------------------
  # CRITICAL: Re-create the CmdStanModel in each worker from the
  # pre-compiled executable path. Do NOT pass the model object directly.

  rep_results <- parallel::mclapply(reps_to_run, function(rep_id) {

    # Re-create CmdStanModel in the forked worker process
    worker_model <- cmdstanr::cmdstan_model(
      exe_file  = exe_path,
      stan_file = stan_file
    )

    tryCatch({
      run_single_rep(
        scenario      = scenario,
        rep_id        = rep_id,
        stan_model    = worker_model,
        mcmc_settings = SIM_PARAMS$mcmc,
        output_dir    = cli_args$outdir,
        thresholds    = c(1.2, 1.5),
        ci_level      = 0.90,
        beta_true     = SIM_PARAMS$beta_true,
        quiet         = TRUE
      )
    }, error = function(e) {
      list(
        scenario_id        = sc_id,
        rep_id             = rep_id,
        convergence_failed = FALSE,
        fitting_failed     = TRUE,
        der = NULL, coverage = NULL, ci_widths = NULL,
        n_corrected = NULL, timing = list(total = NA_real_),
        data_diagnostics = NULL, param_types = NULL,
        seed = NA_integer_,
        error_message = conditionMessage(e)
      )
    })
  }, mc.cores = cli_args$cores)

  t_sc_elapsed <- as.numeric(proc.time()["elapsed"] - t_sc_start)
  scenario_timings[sc_idx] <- t_sc_elapsed

  # -- Tally results ----------------------------------------------------------
  n_run <- length(rep_results)
  n_converged <- sum(sapply(rep_results, function(r)
    !isTRUE(r$fitting_failed) && !isTRUE(r$convergence_failed)))
  n_fit_fail  <- sum(sapply(rep_results, function(r)
    isTRUE(r$fitting_failed)))
  n_conv_fail <- sum(sapply(rep_results, function(r)
    isTRUE(r$convergence_failed)))

  cumulative_reps      <- cumulative_reps + n_run
  cumulative_converged <- cumulative_converged + n_converged
  cumulative_failed    <- cumulative_failed + (n_fit_fail + n_conv_fail)

  cat(sprintf("\n  Results: %d converged, %d fit-failed, %d conv-failed (%.1f sec)\n",
              n_converged, n_fit_fail, n_conv_fail, t_sc_elapsed))

  # Print any error messages
  for (r in rep_results) {
    if (!is.null(r$error_message)) {
      cat(sprintf("    ERROR rep %d: %s\n", r$rep_id, r$error_message))
    }
  }

  write_log("  Completed: %d converged, %d fit-fail, %d conv-fail in %.1f sec",
            n_converged, n_fit_fail, n_conv_fail, t_sc_elapsed)

  # -- Scenario DER summary (converged reps only) ----------------------------
  converged_results <- Filter(
    function(r) !isTRUE(r$fitting_failed) && !isTRUE(r$convergence_failed),
    rep_results
  )

  if (length(converged_results) > 0) {
    der_betas  <- do.call(rbind, lapply(converged_results,
                                        function(r) r$der_beta))
    der_thetas <- sapply(converged_results, function(r) mean(r$der_theta))

    cat(sprintf("  Mean DER (beta):  [%.3f, %.3f, %.3f]\n",
                mean(der_betas[, 1]),
                mean(der_betas[, 2]),
                mean(der_betas[, 3])))
    cat(sprintf("  Mean DER (theta): %.3f\n", mean(der_thetas)))
  }

  # -- Aggregate and save scenario summary ------------------------------------
  all_scenario_rds <- tryCatch({
    scenario_dir <- file.path(cli_args$outdir, sc_id)
    rds_files <- sort(list.files(scenario_dir, pattern = "^rep_\\d+\\.rds$",
                                 full.names = TRUE))
    lapply(rds_files, readRDS)
  }, error = function(e) rep_results)

  if (length(all_scenario_rds) > 0) {
    agg <- tryCatch(
      aggregate_scenario_results(all_scenario_rds),
      error = function(e) {
        write_log("  WARNING: Aggregation failed: %s", e$message)
        NULL
      }
    )
    if (!is.null(agg)) {
      summary_file <- file.path(summary_dir,
                                paste0(sc_id, "_summary.rds"))
      saveRDS(agg, file = summary_file)
      cat(sprintf("  Summary saved: %s\n", basename(summary_file)))
    }
  }

  # -- Cumulative progress ----------------------------------------------------
  t_global_elapsed    <- as.numeric(proc.time()["elapsed"] - t_global_start)
  scenarios_remaining <- n_scenarios - sc_idx
  mean_scenario_time  <- t_global_elapsed / sc_idx
  eta_sec             <- scenarios_remaining * mean_scenario_time

  cat(sprintf("\n  --- Cumulative Progress ---\n"))
  cat(sprintf("    Scenarios: %d/%d completed\n", sc_idx, n_scenarios))
  cat(sprintf("    Reps: %d total (%d converged, %d failed)\n",
              cumulative_reps, cumulative_converged, cumulative_failed))
  cat(sprintf("    Wall time: %.1f min | ETA: %.1f min\n",
              t_global_elapsed / 60, eta_sec / 60))

  write_log("  Progress: %d/%d scenarios, %.1f min elapsed, ETA %.1f min",
            sc_idx, n_scenarios, t_global_elapsed / 60, eta_sec / 60)
}


# =============================================================================
# Section 8: Aggregation Helper (used in Section 7)
# =============================================================================
#
# This function aggregates per-replication results into a scenario-level
# summary for quick inspection and later analysis (sim_05_analyze.R).
# =============================================================================

#' Aggregate results across replications for one scenario
#'
#' @param results List of replication results.
#' @param nominal_level Numeric. Nominal CI level (default 0.90).
#' @return A list with scenario-level summaries.
aggregate_scenario_results <- function(results, nominal_level = 0.90) {
  stopifnot(is.list(results), length(results) > 0)

  scenario_id <- results[[1]]$scenario_id
  n_total     <- length(results)

  # Separate converged vs failed
  converged <- Filter(function(r)
    !r$convergence_failed && !r$fitting_failed, results)
  n_converged <- length(converged)
  n_failed    <- n_total - n_converged

  if (n_converged == 0) {
    warning("No converged replications for scenario: ", scenario_id)
    return(list(scenario_id = scenario_id, n_total = n_total,
                n_converged = 0L, n_failed = n_failed,
                coverage_summary = NULL, der_summary = NULL,
                timing_summary = NULL))
  }

  # -- Coverage aggregation ---------------------------------------------------
  strategies  <- names(converged[[1]]$coverage)
  d           <- length(converged[[1]]$der)
  param_types <- converged[[1]]$param_types

  # Build coverage matrices (R x d) per strategy
  coverage_matrices <- list()
  for (strat in strategies) {
    cov_mat <- matrix(NA, n_converged, d)
    for (r in seq_len(n_converged)) {
      cov_mat[r, ] <- as.logical(converged[[r]]$coverage[[strat]])
    }
    colnames(cov_mat) <- names(converged[[1]]$der)
    coverage_matrices[[strat]] <- cov_mat
  }

  # Compute coverage rates by type and strategy
  coverage_comparison <- do.call(rbind, lapply(names(coverage_matrices),
    function(strat) {
      cov_mat <- coverage_matrices[[strat]]
      cov_rate <- colMeans(cov_mat, na.rm = TRUE)

      types <- unique(param_types)
      do.call(rbind, lapply(types, function(tp) {
        idx <- which(param_types == tp)
        data.frame(
          strategy       = strat,
          type           = tp,
          n_params       = length(idx),
          mean_coverage  = mean(cov_rate[idx]),
          min_coverage   = min(cov_rate[idx]),
          max_coverage   = max(cov_rate[idx]),
          mean_deviation = mean(cov_rate[idx]) - nominal_level,
          stringsAsFactors = FALSE
        )
      }))
    }
  ))
  rownames(coverage_comparison) <- NULL

  # -- DER aggregation --------------------------------------------------------
  der_mat <- matrix(NA_real_, n_converged, d)
  for (r in seq_len(n_converged)) {
    der_mat[r, ] <- converged[[r]]$der
  }
  colnames(der_mat) <- names(converged[[1]]$der)

  der_summary <- data.frame(
    param    = colnames(der_mat),
    type     = param_types,
    der_mean = colMeans(der_mat, na.rm = TRUE),
    der_sd   = apply(der_mat, 2, sd, na.rm = TRUE),
    der_q025 = apply(der_mat, 2, quantile, 0.025, na.rm = TRUE),
    der_q500 = apply(der_mat, 2, quantile, 0.500, na.rm = TRUE),
    der_q975 = apply(der_mat, 2, quantile, 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # -- Timing aggregation -----------------------------------------------------
  total_times <- sapply(converged, function(r) r$timing$total)
  fit_times   <- sapply(converged, function(r) r$timing$fitting)

  timing_summary <- data.frame(
    stage      = c("fitting", "total"),
    mean_sec   = c(mean(fit_times, na.rm = TRUE),
                   mean(total_times, na.rm = TRUE)),
    median_sec = c(median(fit_times, na.rm = TRUE),
                   median(total_times, na.rm = TRUE)),
    max_sec    = c(max(fit_times, na.rm = TRUE),
                   max(total_times, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )

  list(
    scenario_id      = scenario_id,
    n_total          = n_total,
    n_converged      = n_converged,
    n_failed         = n_failed,
    coverage_summary = coverage_comparison,
    der_summary      = der_summary,
    timing_summary   = timing_summary
  )
}


# =============================================================================
# Section 9: Final Summary
# =============================================================================

t_global_total <- as.numeric(proc.time()["elapsed"] - t_global_start)

cat("\n\n")
cat("################################################################\n")
cat("#              SIMULATION RUN COMPLETE                         #\n")
cat("################################################################\n\n")

cat(sprintf("Total scenarios run:       %d\n", n_scenarios))
cat(sprintf("Total reps attempted:      %d\n", cumulative_reps))
cat(sprintf("Total converged:           %d (%.1f%%)\n",
            cumulative_converged,
            100 * cumulative_converged / max(cumulative_reps, 1)))
cat(sprintf("Total failed:              %d (%.1f%%)\n",
            cumulative_failed,
            100 * cumulative_failed / max(cumulative_reps, 1)))
cat(sprintf("Total wall-clock time:     %.1f min (%.2f hrs)\n",
            t_global_total / 60, t_global_total / 3600))
cat(sprintf("Mean time per scenario:    %.1f sec\n",
            mean(scenario_timings[scenario_timings > 0])))

write_log("---")
write_log("RUN COMPLETE: %d scenarios, %d reps (%d converged, %d failed), %.1f min",
          n_scenarios, cumulative_reps, cumulative_converged,
          cumulative_failed, t_global_total / 60)

cat(sprintf("\nFinished: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Output directory: %s\n", cli_args$outdir))
cat(sprintf("Log file: %s\n", log_file))
cat("================================================================\n")
