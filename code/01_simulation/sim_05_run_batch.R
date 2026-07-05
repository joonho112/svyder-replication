# =============================================================================
# sim_05_run_batch.R: CLI batch driver for the simulation grid
# =============================================================================
#
# Purpose : Run (cell x rep) tasks of the Monte Carlo study in parallel with
#           parallel::mclapply, with skip-if-exists resumability, per-task
#           error isolation, and periodic progress reporting. This is the
#           final compute stage: it drives run_single_rep (sim_04) across the
#           grid and writes one rep file per (cell, rep) for the results
#           scripts to aggregate. It runs replications in parallel across
#           local cores and is adaptable to an HPC job array.
#
# Usage (from the replication package root):
#
#   Rscript code/01_simulation/sim_05_run_batch.R \
#       --cells all --reps 1:200 --workers 8 --outdir output/simulation
#
# Options:
#   --cells    "all" (default) or a comma-separated list of cell ids
#              (e.g., "J020_N10_GMOD_SG,J050_N50_GSTR_PW") and/or 1-based
#              cell numbers (row indices of build_scenario_grid()).
#   --reps     Replication spec (default "1:200"): "a:b", "k", or a
#              comma-separated mix (e.g., "1:50,101,150:160").
#   --workers  Number of parallel workers (default 8). Each worker runs
#              ONE chain at a time: PARALLEL_CHAINS is set to 1 unless the
#              caller already exported it (effective threads =
#              workers x PARALLEL_CHAINS should not exceed the core count).
#   --outdir   Output directory (default "output/simulation"; relative
#              paths resolve against the package root).
#   --overwrite  Re-run tasks whose rep files already exist (default:
#              skip-if-exists).
#
# Behavior:
#   * Deterministic task ordering: cells in grid order, reps ascending.
#   * Skip-if-exists: existing rep_<0000>.rds files are counted as skipped.
#   * Hard errors never kill the batch: each task runs inside tryCatch and
#     writes out_dir/<cell_id>/rep_<0000>_ERROR.txt with the seed and
#     message, then the driver moves on.
#   * Gate-failing reps are persisted (gate$pass = FALSE), counted in the
#     tally, and left for the aggregation step to adjudicate.
#   * Progress line every ~25 completions; final tally of
#     done / failed-gate / skipped / errors.
#
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
# =============================================================================

suppressPackageStartupMessages({
  library(parallel)
})


# =============================================================================
# Section 0: Locate the package root and source the pipeline
# =============================================================================

.v2_driver_root <- function() {
  cand <- character(0)
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    script <- tryCatch(normalizePath(sub("^--file=", "", file_arg[1])),
                       error = function(e) NULL)
    if (!is.null(script)) cand <- c(cand, file.path(dirname(script), "..", ".."))
  }
  cand <- c(cand, Sys.getenv("SVYDER_REPLICATION_ROOT", ""), getwd())
  for (d in cand) {
    if (nzchar(d) &&
        file.exists(file.path(d, "code", "01_simulation",
                              "sim_04_run_single.R"))) {
      return(normalizePath(d))
    }
  }
  stop("Cannot locate the replication package root; run from the package ",
       "root or set SVYDER_REPLICATION_ROOT.", call. = FALSE)
}

proj_root <- .v2_driver_root()
setwd(proj_root)

# sim_04_run_single.R bootstraps helpers + config + dgp + sim_02 + sim_03.
source(file.path(proj_root, "code", "01_simulation", "sim_04_run_single.R"))


# =============================================================================
# Section 1: Command-line parsing
# =============================================================================

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- list(cells = "all", reps = "1:200", workers = 8L,
              outdir = "output/simulation", overwrite = FALSE)
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--cells") {
      opt$cells <- args[i + 1]; i <- i + 2
    } else if (a == "--reps") {
      opt$reps <- args[i + 1]; i <- i + 2
    } else if (a == "--workers") {
      opt$workers <- as.integer(args[i + 1]); i <- i + 2
    } else if (a == "--outdir") {
      opt$outdir <- args[i + 1]; i <- i + 2
    } else if (a == "--overwrite") {
      opt$overwrite <- TRUE; i <- i + 1
    } else {
      stop("Unknown argument: ", a, call. = FALSE)
    }
  }
  stopifnot(is.finite(opt$workers), opt$workers >= 1L)
  opt
}

#' Expand a reps spec like "1:200", "7", or "1:50,101,150:160"
parse_reps_spec <- function(spec) {
  parts <- strsplit(spec, ",", fixed = TRUE)[[1]]
  reps <- integer(0)
  for (tok in parts) {
    tok <- trimws(tok)
    if (grepl("^\\d+:\\d+$", tok)) {
      ab <- as.integer(strsplit(tok, ":", fixed = TRUE)[[1]])
      reps <- c(reps, seq.int(ab[1], ab[2]))
    } else if (grepl("^\\d+$", tok)) {
      reps <- c(reps, as.integer(tok))
    } else {
      stop("Cannot parse --reps token: '", tok, "'", call. = FALSE)
    }
  }
  reps <- sort(unique(reps))
  stopifnot(length(reps) >= 1L, all(reps >= 1L))
  reps
}

#' Resolve a cells spec to grid row indices ("all", cell ids, or numbers)
parse_cells_spec <- function(spec, grid) {
  if (identical(tolower(spec), "all")) return(seq_len(nrow(grid)))
  toks <- trimws(strsplit(spec, ",", fixed = TRUE)[[1]])
  idx <- integer(0)
  for (tok in toks) {
    if (grepl("^\\d+$", tok)) {
      k <- as.integer(tok)
      if (k < 1L || k > nrow(grid)) {
        stop("Cell number out of range (1..", nrow(grid), "): ", tok,
             call. = FALSE)
      }
      idx <- c(idx, k)
    } else {
      k <- match(tok, grid$cell_id)
      if (is.na(k)) stop("Unknown cell_id: ", tok, call. = FALSE)
      idx <- c(idx, k)
    }
  }
  sort(unique(idx))   # deterministic: grid order
}

opt <- parse_args()
if (!startsWith(opt$outdir, "/")) opt$outdir <- file.path(proj_root, opt$outdir)


# =============================================================================
# Section 2: Configuration banner and task list
# =============================================================================

# Worker-safe Stan threading unless the caller already chose otherwise
# (workers x PARALLEL_CHAINS must not exceed the available core count).
if (!nzchar(Sys.getenv("PARALLEL_CHAINS"))) Sys.setenv(PARALLEL_CHAINS = "1")

grid      <- build_scenario_grid()
cell_idx  <- parse_cells_spec(opt$cells, grid)
rep_ids   <- parse_reps_spec(opt$reps)

# Deterministic task ordering: cell-major (grid order), reps ascending.
tasks <- data.frame(
  cell_number = rep(cell_idx, each = length(rep_ids)),
  rep_id      = rep(rep_ids, times = length(cell_idx))
)
tasks$cell_id  <- grid$cell_id[tasks$cell_number]
tasks$rep_file <- file.path(opt$outdir, tasks$cell_id,
                            sprintf("rep_%04d.rds", tasks$rep_id))

cat("================================================================\n")
cat("  SIMULATION BATCH DRIVER\n")
cat("  Started :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")
cat(sprintf("  Cells           : %d (%s)\n", length(cell_idx),
            if (length(cell_idx) == nrow(grid)) "all" else
              paste(grid$cell_id[cell_idx], collapse = ", ")))
cat(sprintf("  Reps            : %d (%s)\n", length(rep_ids),
            paste(range(rep_ids), collapse = "..")))
cat(sprintf("  Tasks           : %d\n", nrow(tasks)))
cat(sprintf("  Workers         : %d\n", opt$workers))
cat(sprintf("  PARALLEL_CHAINS : %s\n", Sys.getenv("PARALLEL_CHAINS")))
cat(sprintf("  Output dir      : %s\n", opt$outdir))
cat(sprintf("  Overwrite       : %s\n", opt$overwrite))
cat("\n")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# Section 3: Compile the Stan model once in the master process
# =============================================================================
# cmdstanr model objects must not be shared across forked workers; workers
# re-create a CmdStanModel from the pre-compiled executable path (cheap).

cat("--- Compiling / loading Stan model ---\n")
stan_file  <- file.path(proj_root, SIM_PARAMS$stan_model_path)
stan_model <- compile_stan_model(stan_file, quiet = FALSE)
exe_path   <- stan_model$exe_file()
cat(sprintf("  Executable: %s\n\n", exe_path))


# =============================================================================
# Section 4: Task runner (error-isolated)
# =============================================================================

write_error_file <- function(cell_number, rep_id, msg, outdir, grid) {
  # Defensive by design: a failure to WRITE the error report (e.g., an
  # unwritable outdir) must never escalate into a batch-killing error.
  tryCatch({
    cell_id  <- grid$cell_id[cell_number]
    err_dir  <- file.path(outdir, cell_id)
    dir.create(err_dir, recursive = TRUE, showWarnings = FALSE)
    err_file <- file.path(err_dir, sprintf("rep_%04d_ERROR.txt", rep_id))
    writeLines(c(
      sprintf("timestamp : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      sprintf("cell_id   : %s", cell_id),
      sprintf("rep_id    : %d", rep_id),
      sprintf("seed      : %d", get_rep_seed(cell_number, rep_id)),
      sprintf("error     : %s", msg)
    ), err_file)
    err_file
  }, error = function(e) {
    message(sprintf("[batch] could not write ERROR file for %s rep %d: %s",
                    grid$cell_id[cell_number], rep_id, conditionMessage(e)))
    NA_character_
  })
}

run_task <- function(k) {
  cn  <- tasks$cell_number[k]
  rid <- tasks$rep_id[k]

  # Re-create the CmdStanModel in the forked worker from the executable.
  worker_model <- cmdstanr::cmdstan_model(exe_file = exe_path,
                                          stan_file = stan_file)

  res <- tryCatch(
    run_single_rep(cell = grid[cn, , drop = FALSE], rep_id = rid,
                      out_dir = opt$outdir, quiet = TRUE,
                      stan_model = worker_model, grid = grid,
                      overwrite = opt$overwrite),
    error = function(e) {
      write_error_file(cn, rid, conditionMessage(e), opt$outdir, grid)
      structure(list(msg = conditionMessage(e)), class = "v2_rep_error")
    }
  )

  if (inherits(res, "v2_rep_error")) return(list(status = "error"))
  if (isTRUE(res$skipped))           return(list(status = "skipped"))
  list(status = if (isTRUE(res$gate_pass)) "pass" else "fail_gate")
}


# =============================================================================
# Section 5: Skip-if-exists prescan + chunked parallel execution
# =============================================================================

n_tasks  <- nrow(tasks)
statuses <- character(n_tasks)

if (!opt$overwrite) {
  pre_exists <- file.exists(tasks$rep_file)
  statuses[pre_exists] <- "skipped"
  if (any(pre_exists)) {
    cat(sprintf("--- Resume: %d of %d rep files already exist (skipped) ---\n\n",
                sum(pre_exists), n_tasks))
  }
}

todo <- which(statuses == "")
t_start <- proc.time()["elapsed"]

print_progress <- function(n_done_total) {
  el <- as.numeric(proc.time()["elapsed"] - t_start)
  tab <- table(factor(statuses,
                      levels = c("pass", "fail_gate", "skipped", "error")))
  n_new <- sum(tab[c("pass", "fail_gate", "error")])
  rate  <- if (n_new > 0) el / n_new else NA_real_
  n_left <- sum(statuses == "")   # "" never matches as a names subscript
  eta_min <- if (!is.na(rate)) n_left * rate / 60 else NA_real_
  cat(sprintf(
    "[%s] %5d/%d done | pass %d | gate-fail %d | skipped %d | errors %d | %6.1f min elapsed | ETA %6.1f min\n",
    format(Sys.time(), "%H:%M:%S"), n_done_total, n_tasks,
    tab["pass"], tab["fail_gate"], tab["skipped"], tab["error"],
    el / 60, eta_min))
}

if (length(todo) > 0) {
  chunk_size <- max(25L, opt$workers)   # progress line every ~25 completions
  chunks <- split(todo, ceiling(seq_along(todo) / chunk_size))

  for (ch in chunks) {
    out <- parallel::mclapply(ch, run_task,
                              mc.cores = opt$workers,
                              mc.preschedule = FALSE)
    for (j in seq_along(ch)) {
      r <- out[[j]]
      if (is.null(r) || inherits(r, "try-error") || is.null(r$status)) {
        # Worker died (e.g., OOM/segfault): record from the master side.
        msg <- if (inherits(r, "try-error")) as.character(r) else
          "worker process died (mclapply returned no result)"
        write_error_file(tasks$cell_number[ch[j]], tasks$rep_id[ch[j]],
                            msg, opt$outdir, grid)
        statuses[ch[j]] <- "error"
      } else {
        statuses[ch[j]] <- r$status
      }
    }
    print_progress(sum(statuses != ""))
  }
} else {
  cat("Nothing to do: all requested rep files already exist.\n")
}


# =============================================================================
# Section 6: Final tally
# =============================================================================

t_total <- as.numeric(proc.time()["elapsed"] - t_start)
tab <- table(factor(statuses, levels = c("pass", "fail_gate", "skipped", "error")))

cat("\n================================================================\n")
cat("  BATCH COMPLETE\n")
cat("================================================================\n")
cat(sprintf("  Tasks             : %d\n", n_tasks))
cat(sprintf("  Done (gate pass)  : %d\n", tab["pass"]))
cat(sprintf("  Done (gate FAIL)  : %d  (persisted; aggregation decides)\n",
            tab["fail_gate"]))
cat(sprintf("  Skipped (existing): %d\n", tab["skipped"]))
cat(sprintf("  Hard errors       : %d  (see rep_<id>_ERROR.txt)\n",
            tab["error"]))
cat(sprintf("  Wall time         : %.1f min\n", t_total / 60))
cat(sprintf("  Finished          : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("================================================================\n")

quit(save = "no", status = if (tab["error"] > 0) 1L else 0L)
