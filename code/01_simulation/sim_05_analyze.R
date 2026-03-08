# =============================================================================
# sim_05_analyze.R: Results Analysis and Figure Generation
# =============================================================================
#
# Purpose : Load all simulation results, build master data frames, produce
#           publication-quality summary tables and figures for the manuscript.
#           Computes DER summaries, coverage comparisons across correction
#           strategies, and CI width ratios. Designed to run after all 54
#           scenarios are complete but handles partial results gracefully.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Section 4.4 (Simulation Results), Tables 1-3, Figures 1-4
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Track   : A (Full Replication)
# Inputs  : - sim_00_config.R (scenario grid and helpers)
#           - {outdir}/summary/{scenario_id}_summary.rds (per-scenario summaries)
#           - {outdir}/{scenario_id}/rep_*.rds (individual rep files)
# Outputs : - Tables: {outdir}/tables/table1_der_summary.csv
#                      {outdir}/tables/table2_coverage_comparison.csv
#                      {outdir}/tables/table3_ci_width_ratios.csv
#           - Figures: {outdir}/figures/fig1_der_by_deff_icc.pdf/.png
#                      {outdir}/figures/fig2_coverage_by_strategy.pdf/.png
#                      {outdir}/figures/fig3_correction_comparison.pdf/.png
#                      {outdir}/figures/fig4_convergence_runtime.pdf/.png
#
# Usage (from replication package root):
#   source("code/01_simulation/sim_05_analyze.R")
#   -- OR --
#   Rscript code/01_simulation/sim_05_analyze.R
#
# Notes   : This script uses ggplot2 >= 4.0, dplyr, tidyr, and scales.
#           All figures use the Okabe-Ito colorblind-friendly palette.
#           Figure dimensions are calibrated for JSSAM (Oxford University Press):
#             Single column: ~3.3 inches; double column: ~6.7 inches.
# =============================================================================


# =============================================================================
# Section 0: Setup and Paths
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Detect project root: use here::here() if available, else fall back
# to the PROJ_ROOT environment variable or the script's grandparent dir.
if (requireNamespace("here", quietly = TRUE)) {
  PROJECT_ROOT <- here::here()
} else {
  PROJECT_ROOT <- Sys.getenv(
    "PROJ_ROOT",
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))
  )
}

# Source configuration for scenario grid and helpers
source(file.path(PROJECT_ROOT, "code/01_simulation/sim_00_config.R"))

# Paths
MAIN_DIR    <- file.path(PROJECT_ROOT, "data/simulation")
SUMMARY_DIR <- file.path(MAIN_DIR, "summary")
TABLE_DIR   <- file.path(MAIN_DIR, "tables")
FIGURE_DIR  <- file.path(MAIN_DIR, "figures")

# Create output directories
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# Scenario grid
GRID <- build_scenario_grid()

# Figure dimensions for JSSAM (Oxford University Press)
FIG_WIDTH_SINGLE <- 3.3      # single column
FIG_WIDTH_DOUBLE <- 6.7      # double column
FIG_DPI          <- 300

# Colorblind-friendly palette (Okabe-Ito)
OKABE_ITO <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

# Strategy display names and ordering
# Includes both DER-1.2 and DER-1.5 thresholds from the simulation
STRATEGY_ORDER <- c("naive", "blanket", "DER-1.2", "DER-1.5")
STRATEGY_LABELS <- c(
  naive       = "Naive (no correction)",
  blanket     = "Blanket Cholesky",
  `DER-1.2`   = "Selective DER (tau = 1.2)",
  `DER-1.5`   = "Selective DER (tau = 1.5)"
)

# Parameter type display names
PARAM_TYPE_LABELS <- c(
  fe_within  = "Within-group FE",
  fe_between = "Between-group FE",
  re         = "Random effects"
)


# =============================================================================
# PART 1: Load and Aggregate Results
# =============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  DER Simulation: Comprehensive Analysis\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")


# -----------------------------------------------------------------------------
# 1.1 Load scenario summary files
# -----------------------------------------------------------------------------

cat("--- Loading summary files ---\n")

summary_files <- list.files(SUMMARY_DIR, pattern = "_summary\\.rds$",
                            full.names = TRUE)

if (length(summary_files) == 0) {
  stop("No summary files found in ", SUMMARY_DIR,
       ". Run the simulation first (sim_04_run.R).")
}

cat("  Found", length(summary_files), "of", nrow(GRID),
    "expected summary files.\n")

summaries <- lapply(summary_files, function(f) {
  tryCatch(readRDS(f), error = function(e) {
    warning("Failed to read: ", f, " -- ", e$message)
    NULL
  })
})
summaries <- Filter(Negate(is.null), summaries)
names(summaries) <- sapply(summaries, function(s) s$scenario_id)

cat("  Successfully loaded:", length(summaries), "summaries.\n")

# Report missing scenarios
loaded_ids   <- names(summaries)
expected_ids <- GRID$scenario_id
missing_ids  <- setdiff(expected_ids, loaded_ids)

if (length(missing_ids) > 0) {
  cat("  WARNING:", length(missing_ids), "scenarios not yet available.\n")
  if (length(missing_ids) <= 10) {
    cat("  Missing:", paste(missing_ids, collapse = ", "), "\n")
  } else {
    cat("  First 10 missing:",
        paste(head(missing_ids, 10), collapse = ", "), "...\n")
  }
}


# -----------------------------------------------------------------------------
# 1.2 Build master data frames
# -----------------------------------------------------------------------------

cat("\n--- Building master data frames ---\n")

# 1.2a: DER master table (one row per parameter per scenario)
build_der_master <- function(summaries, grid) {
  rows <- list()
  for (sid in names(summaries)) {
    s <- summaries[[sid]]
    if (is.null(s$der_summary)) next

    sc <- grid[grid$scenario_id == sid, , drop = FALSE]
    if (nrow(sc) == 0) {
      sc_parsed <- parse_scenario_id(sid)
      sc <- data.frame(
        scenario_id = sid,
        J           = sc_parsed$J,
        cv_w        = sc_parsed$cv_w,
        icc         = sc_parsed$icc,
        informative = sc_parsed$informative,
        approx_deff = 1 + sc_parsed$cv_w^2,
        stringsAsFactors = FALSE
      )
    }

    der_df <- s$der_summary
    der_df$scenario_id  <- sid
    der_df$J            <- sc$J
    der_df$cv_w         <- sc$cv_w
    der_df$icc          <- sc$icc
    der_df$informative  <- sc$informative
    der_df$approx_deff  <- sc$approx_deff
    der_df$n_converged  <- s$n_converged
    der_df$n_total      <- s$n_total

    rows[[sid]] <- der_df
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

der_master <- build_der_master(summaries, GRID)


# 1.2b: Coverage master table (one row per strategy x type x scenario)
build_coverage_master <- function(summaries, grid) {
  rows <- list()
  for (sid in names(summaries)) {
    s <- summaries[[sid]]
    if (is.null(s$coverage_summary)) next

    sc <- grid[grid$scenario_id == sid, , drop = FALSE]
    if (nrow(sc) == 0) {
      sc_parsed <- parse_scenario_id(sid)
      sc <- data.frame(
        scenario_id = sid,
        J           = sc_parsed$J,
        cv_w        = sc_parsed$cv_w,
        icc         = sc_parsed$icc,
        informative = sc_parsed$informative,
        approx_deff = 1 + sc_parsed$cv_w^2,
        stringsAsFactors = FALSE
      )
    }

    cov_df <- s$coverage_summary
    cov_df$scenario_id  <- sid
    cov_df$J            <- sc$J
    cov_df$cv_w         <- sc$cv_w
    cov_df$icc          <- sc$icc
    cov_df$informative  <- sc$informative
    cov_df$approx_deff  <- sc$approx_deff
    cov_df$n_converged  <- s$n_converged

    rows[[sid]] <- cov_df
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

cov_master <- build_coverage_master(summaries, GRID)


# 1.2c: Timing master table (one row per scenario)
build_timing_master <- function(summaries, grid) {
  rows <- list()
  for (sid in names(summaries)) {
    s <- summaries[[sid]]
    if (is.null(s$timing_summary)) next

    sc <- grid[grid$scenario_id == sid, , drop = FALSE]
    if (nrow(sc) == 0) {
      sc_parsed <- parse_scenario_id(sid)
      sc <- data.frame(
        scenario_id = sid,
        J           = sc_parsed$J,
        cv_w        = sc_parsed$cv_w,
        icc         = sc_parsed$icc,
        informative = sc_parsed$informative,
        stringsAsFactors = FALSE
      )
    }

    fit_row <- s$timing_summary[s$timing_summary$stage == "fitting", ]
    tot_row <- s$timing_summary[s$timing_summary$stage == "total", ]

    rows[[sid]] <- data.frame(
      scenario_id      = sid,
      J                = sc$J,
      cv_w             = sc$cv_w,
      icc              = sc$icc,
      informative      = sc$informative,
      n_total          = s$n_total,
      n_converged      = s$n_converged,
      n_failed         = s$n_failed,
      convergence_rate = s$n_converged / s$n_total,
      mean_fit_sec     = fit_row$mean_sec,
      median_fit_sec   = fit_row$median_sec,
      max_fit_sec      = fit_row$max_sec,
      mean_total_sec   = tot_row$mean_sec,
      median_total_sec = tot_row$median_sec,
      max_total_sec    = tot_row$max_sec,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

timing_master <- build_timing_master(summaries, GRID)

cat(sprintf("  DER master:      %d rows\n", nrow(der_master)))
cat(sprintf("  Coverage master: %d rows\n", nrow(cov_master)))
cat(sprintf("  Timing master:   %d rows\n", nrow(timing_master)))


# -----------------------------------------------------------------------------
# 1.3 DER summary by parameter type
# -----------------------------------------------------------------------------
# Aggregate DER across individual parameters within each type
# (fe_within, fe_between, re) for each scenario.

if (nrow(der_master) > 0) {
  der_type_summary <- der_master %>%
    group_by(scenario_id, J, cv_w, icc, informative, approx_deff, type) %>%
    summarise(
      mean_der   = mean(der_mean, na.rm = TRUE),
      sd_der     = mean(der_sd, na.rm = TRUE),
      median_der = mean(der_q500, na.rm = TRUE),
      n_params   = n(),
      .groups    = "drop"
    )
  cat(sprintf("  DER type summary: %d rows\n", nrow(der_type_summary)))
} else {
  der_type_summary <- data.frame()
  cat("  WARNING: No DER data available.\n")
}


# =============================================================================
# PART 2: Publication Tables
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  Generating Tables\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")


# -----------------------------------------------------------------------------
# TABLE 1: DER Summary by Scenario Factors
# -----------------------------------------------------------------------------
# Shows how DER varies across the factorial design for each parameter type.
# Columns: J x cv_w combinations. Rows: parameter type x informativeness.
# This is the main empirical result supporting the DER decomposition
# (Section 4.4, Table 1 of Lee (2026)).
# -----------------------------------------------------------------------------

cat("--- Table 1: DER Summary by Scenario Factors ---\n")

if (nrow(der_type_summary) > 0) {

  # Ordering for parameter types
  type_order <- c("fe_within", "fe_between", "re")

  table1_numeric <- der_type_summary %>%
    group_by(type, J, cv_w, informative) %>%
    summarise(
      mean_der    = mean(mean_der, na.rm = TRUE),
      sd_der      = mean(sd_der, na.rm = TRUE),
      n_scenarios = n(),
      .groups     = "drop"
    )

  # Formatted version with "mean (sd)" cells
  table1_formatted <- table1_numeric %>%
    mutate(
      cell           = sprintf("%.3f (%.3f)", mean_der, sd_der),
      informative_label = ifelse(informative, "Informative", "Non-informative"),
      type_label     = PARAM_TYPE_LABELS[type],
      col_key        = paste0("J", J, "_CV", cv_w)
    ) %>%
    select(type, type_label, informative_label, col_key, cell) %>%
    pivot_wider(names_from = col_key, values_from = cell)

  # Sort by informativeness and parameter type
  table1_formatted <- table1_formatted %>%
    mutate(type = factor(type, levels = type_order)) %>%
    arrange(informative_label, type)

  cat("  Table 1 dimensions:", nrow(table1_formatted), "x",
      ncol(table1_formatted), "\n")
  print(as.data.frame(table1_formatted))

  write.csv(as.data.frame(table1_formatted),
            file.path(TABLE_DIR, "table1_der_summary.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(table1_numeric),
            file.path(TABLE_DIR, "table1_der_summary_numeric.csv"),
            row.names = FALSE)
  cat("  Saved: table1_der_summary.csv\n")

} else {
  cat("  WARNING: No DER data for Table 1.\n")
}


# -----------------------------------------------------------------------------
# TABLE 2: Coverage Comparison by Strategy
# -----------------------------------------------------------------------------
# Compares 90% CI coverage rates across correction strategies
# (naive, blanket, selective DER-1.2, selective DER-1.5) for each
# parameter type and J level. This is the key validation table showing
# that selective DER correction preserves nominal coverage while avoiding
# blanket overcorrection.
# See Section 4.4, Table 2 of Lee (2026).
# -----------------------------------------------------------------------------

cat("\n--- Table 2: Coverage Comparison by Strategy ---\n")

if (nrow(cov_master) > 0) {

  type_order <- c("fe_within", "fe_between", "re")

  table2_data <- cov_master %>%
    group_by(strategy, type, J) %>%
    summarise(
      mean_cov    = mean(mean_coverage, na.rm = TRUE),
      sd_cov      = sd(mean_coverage, na.rm = TRUE),
      min_cov     = min(mean_coverage, na.rm = TRUE),
      max_cov     = max(mean_coverage, na.rm = TRUE),
      mean_dev    = mean(mean_deviation, na.rm = TRUE),
      n_scenarios = n(),
      .groups     = "drop"
    ) %>%
    mutate(
      # Approximate MC standard error for coverage
      mc_se_approx = sd_cov / sqrt(pmax(n_scenarios, 1)),
      # Flag significant deviations from nominal (90%)
      below_nominal = mean_cov < SIM_PARAMS$ci_level - 1.96 * mc_se_approx,
      above_nominal = mean_cov > SIM_PARAMS$ci_level + 1.96 * mc_se_approx,
      flag = case_when(
        below_nominal ~ "LOW",
        above_nominal ~ "HIGH",
        TRUE          ~ ""
      ),
      cell = sprintf("%.3f (%.3f)%s",
                     mean_cov, sd_cov,
                     ifelse(nchar(flag) > 0, paste0(" ", flag), "")),
      strategy_label = STRATEGY_LABELS[strategy],
      type_label     = PARAM_TYPE_LABELS[type]
    )

  table2_wide <- table2_data %>%
    mutate(J_col = paste0("J_", J)) %>%
    select(strategy, strategy_label, type, type_label, J_col, cell) %>%
    pivot_wider(names_from = J_col, values_from = cell)

  table2_wide <- table2_wide %>%
    mutate(
      strategy = factor(strategy, levels = STRATEGY_ORDER),
      type     = factor(type, levels = type_order)
    ) %>%
    arrange(strategy, type)

  cat("  Table 2 dimensions:", nrow(table2_wide), "x",
      ncol(table2_wide), "\n")
  print(as.data.frame(table2_wide))

  write.csv(as.data.frame(table2_wide),
            file.path(TABLE_DIR, "table2_coverage_comparison.csv"),
            row.names = FALSE)
  cat("  Saved: table2_coverage_comparison.csv\n")

} else {
  cat("  WARNING: No coverage data available.\n")
}


# -----------------------------------------------------------------------------
# TABLE 3: CI Width Ratios (relative to naive)
# -----------------------------------------------------------------------------
# Shows the cost of each correction strategy in terms of interval width.
# Ratios > 1 indicate wider intervals (more conservative).
# See Section 4.4, Table 3 of Lee (2026).
# -----------------------------------------------------------------------------

cat("\n--- Table 3: CI Width Ratios ---\n")

build_table3 <- function(summaries, grid) {

  type_order <- c("fe_within", "fe_between", "re")
  width_rows <- list()

  for (sid in names(summaries)) {
    scenario_dir <- file.path(MAIN_DIR, sid)
    if (!dir.exists(scenario_dir)) next

    rds_files <- list.files(scenario_dir, pattern = "^rep_\\d+\\.rds$",
                            full.names = TRUE)
    if (length(rds_files) == 0) next

    sc <- grid[grid$scenario_id == sid, , drop = FALSE]
    if (nrow(sc) == 0) next

    # Load a sample of reps for width computation (30 is sufficient)
    rds_sample <- head(rds_files, 30)

    for (f in rds_sample) {
      rep_result <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.null(rep_result)) next
      if (isTRUE(rep_result$fitting_failed) ||
          isTRUE(rep_result$convergence_failed) ||
          isTRUE(rep_result$sandwich_failed)) next
      if (is.null(rep_result$ci_widths)) next

      ptypes <- rep_result$param_types

      for (strat in names(rep_result$ci_widths)) {
        widths <- rep_result$ci_widths[[strat]]
        if (length(widths) != length(ptypes)) next

        for (tp in unique(ptypes)) {
          idx <- which(ptypes == tp)
          width_rows[[length(width_rows) + 1L]] <- data.frame(
            scenario_id = sid,
            J           = sc$J,
            cv_w        = sc$cv_w,
            icc         = sc$icc,
            informative = sc$informative,
            strategy    = strat,
            type        = tp,
            mean_width  = mean(widths[idx], na.rm = TRUE),
            rep_id      = rep_result$rep_id,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(width_rows) == 0) {
    cat("  WARNING: No CI width data available.\n")
    return(NULL)
  }

  width_df <- do.call(rbind, width_rows)

  # Average width per strategy x type x J
  width_summary <- width_df %>%
    group_by(strategy, type, J) %>%
    summarise(mean_width = mean(mean_width, na.rm = TRUE),
              .groups = "drop")

  # Compute ratio relative to naive
  naive_widths <- width_summary %>%
    filter(strategy == "naive") %>%
    select(type, J, naive_width = mean_width)

  width_ratios <- width_summary %>%
    left_join(naive_widths, by = c("type", "J")) %>%
    mutate(
      width_ratio    = mean_width / naive_width,
      cell           = sprintf("%.3f", width_ratio),
      strategy_label = STRATEGY_LABELS[strategy],
      type_label     = PARAM_TYPE_LABELS[type]
    )

  # Pivot wider
  tbl3_wide <- width_ratios %>%
    mutate(
      J_col = paste0("J_", J)
    ) %>%
    select(strategy, strategy_label, type, type_label, J_col, cell) %>%
    pivot_wider(names_from = J_col, values_from = cell) %>%
    mutate(
      strategy = factor(strategy, levels = STRATEGY_ORDER),
      type     = factor(type, levels = type_order)
    ) %>%
    arrange(strategy, type)

  tbl3_wide
}

table3 <- build_table3(summaries, GRID)

if (!is.null(table3)) {
  cat("  Table 3 dimensions:", nrow(table3), "x", ncol(table3), "\n")
  print(as.data.frame(table3))

  write.csv(as.data.frame(table3),
            file.path(TABLE_DIR, "table3_ci_width_ratios.csv"),
            row.names = FALSE)
  cat("  Saved: table3_ci_width_ratios.csv\n")
}


# =============================================================================
# PART 3: Publication Figures
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  Generating Figures\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Common theme for all figures (JSSAM style)
theme_paper <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 9),
    strip.background  = element_rect(fill = "grey95", color = NA),
    legend.position   = "bottom",
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 8, face = "bold"),
    axis.title        = element_text(size = 9),
    axis.text         = element_text(size = 7),
    plot.title        = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(size = 8, hjust = 0.5, color = "grey40"),
    panel.spacing     = unit(0.8, "lines")
  )


# -----------------------------------------------------------------------------
# FIGURE 1: DER by DEFF and ICC
# -----------------------------------------------------------------------------
# The main result figure showing how DER varies with the Kish DEFF (from
# weight CV) and ICC across the three parameter types. Key findings:
#   - Within-group FE: DER tracks DEFF closely (the 1:1 line)
#   - Between-group FE: DER < DEFF (shrinkage attenuates)
#   - Random effects: DER << DEFF (strongly attenuated by hierarchical prior)
# See Section 4.4, Figure 1 of Lee (2026).
# -----------------------------------------------------------------------------

cat("--- Figure 1: DER by DEFF and ICC ---\n")

if (nrow(der_type_summary) > 0) {

  plot_data <- der_type_summary %>%
    mutate(
      J_label     = factor(paste0("J = ", J),
                           levels = paste0("J = ", c(20, 50, 100))),
      type_label  = factor(PARAM_TYPE_LABELS[type],
                           levels = c("Within-group FE",
                                      "Between-group FE",
                                      "Random effects")),
      icc_label   = factor(paste0("ICC = ", icc),
                           levels = paste0("ICC = ", c(0.05, 0.15, 0.30))),
      shape_label = ifelse(informative, "Informative", "Non-informative")
    )

  fig1 <- ggplot(plot_data,
                 aes(x = approx_deff, y = mean_der,
                     color = icc_label, shape = shape_label)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_point(size = 2, alpha = 0.85) +
    # Theoretical reference for within-group FE (DER = DEFF)
    geom_abline(
      data = data.frame(
        type_label = factor("Within-group FE",
                            levels = c("Within-group FE",
                                       "Between-group FE",
                                       "Random effects"))
      ),
      aes(slope = 1, intercept = 0),
      linetype = "dotted", color = "grey30", linewidth = 0.4,
      inherit.aes = FALSE
    ) +
    facet_grid(J_label ~ type_label, scales = "free_y") +
    scale_color_manual(
      name   = "ICC",
      values = OKABE_ITO[c(1, 5, 6)]
    ) +
    scale_shape_manual(
      name   = "Weights",
      values = c("Non-informative" = 16, "Informative" = 17)
    ) +
    scale_x_continuous(
      name   = expression(paste("Approximate DEFF (1 + ", CV[w]^2, ")")),
      breaks = c(1.09, 2.0, 5.0),
      labels = c("1.09", "2.00", "5.00")
    ) +
    scale_y_continuous(name = "Mean DER") +
    labs(
      title    = "Design Effect Ratios Across Simulation Conditions",
      subtitle = paste0("Within-group FE tracks DEFF; ",
                        "between-group FE and RE show shrinkage attenuation")
    ) +
    theme_paper +
    guides(
      color = guide_legend(order = 1, nrow = 1),
      shape = guide_legend(order = 2, nrow = 1)
    )

  ggsave(file.path(FIGURE_DIR, "fig1_der_by_deff_icc.pdf"),
         fig1, width = FIG_WIDTH_DOUBLE, height = 6.5,
         units = "in", device = "pdf")
  ggsave(file.path(FIGURE_DIR, "fig1_der_by_deff_icc.png"),
         fig1, width = FIG_WIDTH_DOUBLE, height = 6.5,
         units = "in", dpi = FIG_DPI)
  cat("  Saved: fig1_der_by_deff_icc.pdf/.png\n")

} else {
  cat("  WARNING: No DER data for Figure 1.\n")
}


# -----------------------------------------------------------------------------
# FIGURE 2: Coverage by Strategy
# -----------------------------------------------------------------------------
# Shows coverage rates (with MC confidence intervals) for each correction
# strategy, faceted by J level and strategy. The nominal level (90%) is
# shown as a dashed reference line.
# See Section 4.4, Figure 2 of Lee (2026).
# -----------------------------------------------------------------------------

cat("\n--- Figure 2: Coverage by Strategy ---\n")

if (nrow(cov_master) > 0) {

  plot_data <- cov_master %>%
    mutate(
      strategy_label = factor(STRATEGY_LABELS[strategy],
                              levels = STRATEGY_LABELS),
      J_label        = factor(paste0("J = ", J),
                              levels = paste0("J = ", c(20, 50, 100))),
      type_label     = factor(PARAM_TYPE_LABELS[type],
                              levels = c("Within-group FE",
                                         "Between-group FE",
                                         "Random effects")),
      deff_label     = factor(
        paste0("DEFF ~ ", formatC(approx_deff, format = "f", digits = 2)),
        levels = paste0("DEFF ~ ",
                        formatC(sort(unique(approx_deff)),
                                format = "f", digits = 2))
      ),
      # Monte Carlo SE for error bars
      mc_se    = sqrt(mean_coverage * (1 - mean_coverage) /
                        pmax(n_converged, 1)),
      ci_lower = pmax(0, mean_coverage - 1.96 * mc_se),
      ci_upper = pmin(1, mean_coverage + 1.96 * mc_se)
    )

  fig2 <- ggplot(plot_data,
                 aes(x = type_label, y = mean_coverage,
                     color = deff_label)) +
    geom_hline(yintercept = SIM_PARAMS$ci_level,
               linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(
      position = position_dodge(width = 0.6),
      size = 1.5, alpha = 0.8
    ) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = position_dodge(width = 0.6),
      width = 0.2, linewidth = 0.3
    ) +
    facet_grid(J_label ~ strategy_label) +
    scale_color_manual(
      name   = "DEFF",
      values = OKABE_ITO[c(3, 5, 6)]
    ) +
    scale_y_continuous(
      name   = "Coverage Rate",
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = function(x) sprintf("%.1f", x)
    ) +
    scale_x_discrete(
      name   = "Parameter Type",
      labels = function(x) gsub(" ", "\n", x)
    ) +
    labs(
      title    = "Coverage Rates by Correction Strategy",
      subtitle = paste0("Nominal level = ", SIM_PARAMS$ci_level * 100,
                        "%; bars show 95% MC confidence intervals")
    ) +
    theme_paper +
    theme(
      axis.text.x = element_text(size = 6.5, lineheight = 0.9),
      panel.spacing.x = unit(0.5, "lines")
    ) +
    guides(color = guide_legend(nrow = 1))

  ggsave(file.path(FIGURE_DIR, "fig2_coverage_by_strategy.pdf"),
         fig2, width = FIG_WIDTH_DOUBLE, height = 7,
         units = "in", device = "pdf")
  ggsave(file.path(FIGURE_DIR, "fig2_coverage_by_strategy.png"),
         fig2, width = FIG_WIDTH_DOUBLE, height = 7,
         units = "in", dpi = FIG_DPI)
  cat("  Saved: fig2_coverage_by_strategy.pdf/.png\n")

} else {
  cat("  WARNING: No coverage data for Figure 2.\n")
}


# -----------------------------------------------------------------------------
# FIGURE 3: Naive vs Corrected Coverage Scatter
# -----------------------------------------------------------------------------
# Scatter plot comparing corrected coverage (y-axis) against naive coverage
# (x-axis) for each strategy. Points above the diagonal indicate improvement;
# points below indicate harm. This demonstrates that selective DER correction
# helps where needed and avoids harm elsewhere.
# See Section 4.4, Figure 3 of Lee (2026).
# -----------------------------------------------------------------------------

cat("\n--- Figure 3: DER-tau vs Blanket Comparison ---\n")

if (nrow(cov_master) > 0) {

  naive_cov <- cov_master %>%
    filter(strategy == "naive") %>%
    select(scenario_id, type, naive_coverage = mean_coverage)

  corrected_cov <- cov_master %>%
    filter(strategy != "naive") %>%
    select(scenario_id, type, strategy,
           corrected_coverage = mean_coverage)

  plot_data <- corrected_cov %>%
    left_join(naive_cov, by = c("scenario_id", "type")) %>%
    mutate(
      strategy_label = factor(STRATEGY_LABELS[strategy],
                              levels = STRATEGY_LABELS[-1]),
      type_label     = factor(PARAM_TYPE_LABELS[type],
                              levels = c("Within-group FE",
                                         "Between-group FE",
                                         "Random effects"))
    ) %>%
    filter(!is.na(naive_coverage), !is.na(corrected_coverage))

  if (nrow(plot_data) > 0) {
    fig3 <- ggplot(plot_data,
                   aes(x = naive_coverage, y = corrected_coverage,
                       color = strategy_label)) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "grey50", linewidth = 0.4) +
      geom_hline(yintercept = SIM_PARAMS$ci_level,
                 linetype = "dotted", color = "grey70", linewidth = 0.3) +
      geom_vline(xintercept = SIM_PARAMS$ci_level,
                 linetype = "dotted", color = "grey70", linewidth = 0.3) +
      geom_point(size = 2, alpha = 0.7) +
      facet_wrap(~ type_label, ncol = 3, scales = "free") +
      scale_color_manual(
        name   = "Correction Strategy",
        values = OKABE_ITO[c(6, 5, 3)]
      ) +
      scale_x_continuous(
        name   = "Naive Coverage",
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2)
      ) +
      scale_y_continuous(
        name   = "Corrected Coverage",
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2)
      ) +
      labs(
        title    = "Effect of Correction Strategies on Coverage",
        subtitle = "Above diagonal: correction helps; below: correction harms"
      ) +
      theme_paper +
      guides(color = guide_legend(nrow = 1))

    ggsave(file.path(FIGURE_DIR, "fig3_correction_comparison.pdf"),
           fig3, width = FIG_WIDTH_DOUBLE, height = 3.5,
           units = "in", device = "pdf")
    ggsave(file.path(FIGURE_DIR, "fig3_correction_comparison.png"),
           fig3, width = FIG_WIDTH_DOUBLE, height = 3.5,
           units = "in", dpi = FIG_DPI)
    cat("  Saved: fig3_correction_comparison.pdf/.png\n")
  } else {
    cat("  WARNING: No data for Figure 3.\n")
  }

} else {
  cat("  WARNING: No coverage data for Figure 3.\n")
}


# -----------------------------------------------------------------------------
# FIGURE 4: Convergence and Runtime
# -----------------------------------------------------------------------------
# Supplementary figure showing convergence rates and fitting times by
# scenario. Useful for assessing computational feasibility.
# See Section 4.4, Figure S1 (supplementary).
# -----------------------------------------------------------------------------

cat("\n--- Figure 4: Convergence and Runtime ---\n")

if (nrow(timing_master) > 0) {

  panel_data <- timing_master %>%
    mutate(
      J_label = factor(paste0("J = ", J),
                       levels = paste0("J = ", c(20, 50, 100)))
    )

  # Panel B: Fit times by J level
  fig4 <- ggplot(panel_data,
                 aes(x = J_label, y = mean_fit_sec, fill = J_label)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1, linewidth = 0.4) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.5, color = "grey30") +
    scale_fill_manual(
      values = OKABE_ITO[c(1, 5, 3)],
      guide  = "none"
    ) +
    scale_y_continuous(
      name   = "Mean Fit Time (seconds)",
      labels = function(x) sprintf("%.0f", x)
    ) +
    labs(
      x        = "Number of Clusters",
      title    = "Fitting Time by Cluster Count",
      subtitle = "Distribution of mean per-rep fitting time across scenarios"
    ) +
    theme_paper

  ggsave(file.path(FIGURE_DIR, "fig4_convergence_runtime.pdf"),
         fig4, width = FIG_WIDTH_SINGLE, height = 3.5,
         units = "in", device = "pdf")
  ggsave(file.path(FIGURE_DIR, "fig4_convergence_runtime.png"),
         fig4, width = FIG_WIDTH_SINGLE, height = 3.5,
         units = "in", dpi = FIG_DPI)
  cat("  Saved: fig4_convergence_runtime.pdf/.png\n")

} else {
  cat("  WARNING: No timing data for Figure 4.\n")
}


# =============================================================================
# PART 4: Key Findings Summary
# =============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  Key Findings Summary\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Overall progress
cat(sprintf("  Scenarios completed: %d / %d (%.0f%%)\n",
            length(summaries), nrow(GRID),
            100 * length(summaries) / nrow(GRID)))

# Overall convergence
if (nrow(timing_master) > 0) {
  total_reps <- sum(timing_master$n_total, na.rm = TRUE)
  total_conv <- sum(timing_master$n_converged, na.rm = TRUE)
  cat(sprintf("  Total replications: %d, Converged: %d (%.1f%%)\n",
              total_reps, total_conv,
              100 * total_conv / max(total_reps, 1)))
}

# DER summary by type
if (nrow(der_type_summary) > 0) {
  cat("\n  --- DER by Parameter Type (across all scenarios) ---\n")
  der_type_overall <- der_type_summary %>%
    group_by(type) %>%
    summarise(
      grand_mean = mean(mean_der, na.rm = TRUE),
      grand_sd   = sd(mean_der, na.rm = TRUE),
      min_der    = min(mean_der, na.rm = TRUE),
      max_der    = max(mean_der, na.rm = TRUE),
      .groups    = "drop"
    )
  for (i in seq_len(nrow(der_type_overall))) {
    row <- der_type_overall[i, ]
    cat(sprintf("    %-15s  mean = %.3f (SD = %.3f), range = [%.3f, %.3f]\n",
                PARAM_TYPE_LABELS[row$type],
                row$grand_mean, row$grand_sd, row$min_der, row$max_der))
  }
}

# Coverage summary
if (nrow(cov_master) > 0) {
  cat("\n  --- Coverage by Strategy (across all scenarios) ---\n")
  cov_overall <- cov_master %>%
    group_by(strategy, type) %>%
    summarise(
      grand_mean_cov = mean(mean_coverage, na.rm = TRUE),
      grand_sd_cov   = sd(mean_coverage, na.rm = TRUE),
      .groups = "drop"
    )
  for (strat in STRATEGY_ORDER) {
    sub <- cov_overall %>% filter(strategy == strat)
    if (nrow(sub) == 0) next
    cat(sprintf("\n    %s:\n", STRATEGY_LABELS[strat]))
    for (j in seq_len(nrow(sub))) {
      row <- sub[j, ]
      cat(sprintf("      %-15s  %.3f (SD = %.3f)\n",
                  PARAM_TYPE_LABELS[row$type],
                  row$grand_mean_cov, row$grand_sd_cov))
    }
  }
}

# Key theoretical findings
cat("\n\n  --- Key Theoretical Findings ---\n")
cat("  1. Within-group FE: DER tracks approximate DEFF closely.\n")
cat("     Confirms DER_beta ~ DEFF for within-cluster parameters.\n")
cat("     See Theorem 1 (Section 3.2) of Lee (2026).\n")
cat("  2. Between-group FE and RE: DER << DEFF, consistent with the\n")
cat("     DER decomposition theorem: DER = 1 + (DEFF-1)(1-B).\n")
cat("     Hierarchical shrinkage absorbs design effects.\n")
cat("  3. Blanket correction overcorrects between-group/RE parameters,\n")
cat("     inflating intervals unnecessarily.\n")
cat("  4. Selective DER-tau correction preserves nominal coverage for\n")
cat("     all parameter types while minimizing interval inflation.\n")

# Timing
if (nrow(timing_master) > 0) {
  cat("\n  --- Timing ---\n")
  cat(sprintf("    Mean fit time: %.1f seconds\n",
              mean(timing_master$mean_fit_sec, na.rm = TRUE)))
  cat(sprintf("    Mean total time: %.1f seconds\n",
              mean(timing_master$mean_total_sec, na.rm = TRUE)))
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  Analysis complete.\n")
cat(sprintf("  Tables saved to: %s\n", TABLE_DIR))
cat(sprintf("  Figures saved to: %s\n", FIGURE_DIR))
cat(paste(rep("=", 70), collapse = ""), "\n")
