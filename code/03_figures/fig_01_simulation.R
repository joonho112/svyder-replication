# =============================================================================
# fig_01_simulation.R: Simulation Study Figures (Figures 2-3)
# =============================================================================
#
# Purpose : Generate Figures 2 and 3 from the simulation study (Paper Section 4).
#           Figure 2: Formula accuracy (Laplace vs Exact, Simplified vs Exact)
#           Figure 3: Coverage comparison across correction strategies
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Track B : Uses pre-computed simulation summaries from precomputed/.
#           No restricted data or model refitting required.
#
# Inputs  :
#   data/precomputed/simulation/analysis/analysis_all_j.rds  -- Aggregated results
#
# Outputs :
#   output/figures/figure2_formula_accuracy.pdf / .png
#   output/figures/figure3_coverage.pdf / .png
#
# =============================================================================


###############################################################################
## Section 0 : Setup
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

## Source plotting helpers
source(file.path(here::here(), "code", "helpers", "plotting_helpers.R"))

cat("=== Generating Simulation Figures (Track B) ===\n\n")


###############################################################################
## Section 1 : Load Pre-Computed Data
###############################################################################

cat("Loading pre-computed simulation summaries ...\n")

PRECOMP_DIR <- file.path(here::here(), "data", "precomputed")

## Aggregated analysis containing der_table and cov_table
analysis_path <- file.path(PRECOMP_DIR, "simulation", "analysis",
                           "analysis_all_j.rds")

if (!file.exists(analysis_path)) {
  stop("Pre-computed simulation data not found:\n  ", analysis_path,
       "\nPlease ensure data/precomputed/simulation/analysis/analysis_all_j.rds exists.\n",
       "See data/DATA_ACCESS.md for instructions.")
}

analysis <- readRDS(analysis_path)

cat("  Loaded analysis_all_j.rds  [OK]\n\n")


###############################################################################
## Section 2 : Figure 2 -- Formula Accuracy
###############################################################################
##
## Two-panel scatter plot comparing DER computation methods:
##   Panel (a): Laplace-based DER vs Exact DER (Algorithm 1)
##              Expected: tight agreement along 45-degree line.
##   Panel (b): Simplified formulas (Theorems 1-2) vs Exact DER
##              Expected: agreement for within-FE, divergence for
##              between-FE and random effects.
##
## Paper Section 4.2, Figure 2.
##
###############################################################################

cat("--- Creating Figure 2: Formula Accuracy ---\n")

## Representative scenarios spanning the design space
fig2_scenarios <- c(
  "J020_CV100_ICC015_NI",
  "J050_CV200_ICC030_NI",
  "J100_CV030_ICC005_NI",
  "J100_CV100_ICC015_IN",
  "J020_CV200_ICC030_IN",
  "J100_CV200_ICC005_NI"
)

## Read individual replicate files (replicates 1-3 per scenario)
DATA_DIR <- file.path(here::here(), "data", "precomputed", "simulation")

fig2_data <- do.call(rbind, lapply(fig2_scenarios, function(sc) {
  rep_ids <- sprintf("rep_%04d", 1:3)
  do.call(rbind, lapply(rep_ids, function(rid) {
    fpath <- file.path(DATA_DIR, sc, paste0(rid, ".rds"))
    if (!file.exists(fpath)) return(NULL)
    rep_data <- readRDS(fpath)

    exact   <- rep_data$der
    laplace <- rep_data$der_laplace
    types   <- rep_data$param_types
    nms     <- names(exact)

    ## Build theory vector by mapping Theorem predictions to params
    theory <- numeric(length(exact))
    re_counter <- 0L
    for (i in seq_along(exact)) {
      if (types[i] == "fe_within") {
        theory[i] <- rep_data$der_theory$der_fe_within_theory
      } else if (types[i] == "fe_between") {
        theory[i] <- rep_data$der_theory$der_fe_between_theory
      } else {
        re_counter <- re_counter + 1L
        theory[i]  <- rep_data$der_theory$der_re_theory[re_counter]
      }
    }

    data.frame(
      scenario = sc,
      rep      = rid,
      param    = nms,
      type     = types,
      exact    = as.numeric(exact),
      laplace  = as.numeric(laplace),
      theory   = theory,
      stringsAsFactors = FALSE
    )
  }))
}))

if (is.null(fig2_data) || nrow(fig2_data) == 0L) {
  cat("  [SKIP] Individual replicate files not found in precomputed/simulation/.\n")
  cat("  Figure 2 requires replicate-level data.\n\n")
} else {

  ## Filter non-positive values (log scale requires > 0)
  fig2_data <- fig2_data %>%
    filter(exact > 0, laplace > 0, theory > 0)

  ## Summary statistics for annotations
  median_re_laplace <- fig2_data %>%
    mutate(rel_error = abs(laplace - exact) / exact * 100) %>%
    summarise(med_re = median(rel_error)) %>%
    pull(med_re)

  re_med_ratio <- fig2_data %>%
    filter(type == "re") %>%
    summarise(ratio = median(theory / exact)) %>%
    pull(ratio)

  fb_med_ratio <- fig2_data %>%
    filter(type == "fe_between") %>%
    summarise(ratio = median(theory / exact)) %>%
    pull(ratio)

  ## Shared axis range across both panels
  all_vals <- c(fig2_data$exact, fig2_data$laplace, fig2_data$theory)
  ax_lo <- 10^floor(log10(min(all_vals, na.rm = TRUE)))
  ax_hi <- 10^ceiling(log10(max(all_vals, na.rm = TRUE)))

  ## Panel (a): Laplace approximation vs Exact
  p2a <- ggplot(fig2_data, aes(x = exact, y = laplace,
                                color = type, shape = type)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey35", linewidth = 0.6) +
    geom_point(data = ~ filter(., type == "re"),
               size = 1.4, alpha = 0.30) +
    geom_point(data = ~ filter(., type != "re"),
               size = 2.2, alpha = 0.80) +
    geom_rug(data = ~ filter(., type == "re"),
             alpha = 0.12, linewidth = 0.3, sides = "bl",
             length = unit(0.02, "npc")) +
    geom_rug(data = ~ filter(., type != "re"),
             alpha = 0.35, linewidth = 0.4, sides = "bl",
             length = unit(0.025, "npc")) +
    annotate("text",
             x = ax_lo * 2, y = ax_hi / 2,
             label = paste0("Median |RE| = ",
                            formatC(median_re_laplace, format = "f",
                                    digits = 1), "%"),
             size = 3.0, fontface = "italic", color = "grey30",
             hjust = 0, vjust = 1) +
    scale_x_log10(limits = c(ax_lo, ax_hi),
                  labels = label_number(drop0trailing = TRUE)) +
    scale_y_log10(limits = c(ax_lo, ax_hi),
                  labels = label_number(drop0trailing = TRUE)) +
    scale_color_manual(values = TYPE_COLORS, labels = TYPE_LABELS,
                       name = "Parameter type") +
    scale_shape_manual(values = TYPE_SHAPES, labels = TYPE_LABELS,
                       name = "Parameter type") +
    labs(
      title = "(a) General formula",
      x = "Exact DER (Algorithm 1)",
      y = "Laplace-based DER"
    ) +
    theme_pub() +
    theme(legend.position = "none")

  ## Panel (b): Simplified formula vs Exact
  p2b <- ggplot(fig2_data, aes(x = exact, y = theory,
                                color = type, shape = type)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey35", linewidth = 0.6) +
    geom_point(data = ~ filter(., type == "re"),
               size = 1.4, alpha = 0.25) +
    geom_point(data = ~ filter(., type != "re"),
               size = 2.2, alpha = 0.80) +
    geom_rug(data = ~ filter(., type == "re"),
             alpha = 0.12, linewidth = 0.3, sides = "bl",
             length = unit(0.02, "npc")) +
    geom_rug(data = ~ filter(., type != "re"),
             alpha = 0.35, linewidth = 0.4, sides = "bl",
             length = unit(0.025, "npc")) +
    annotate("text",
             x = ax_hi / 2, y = ax_lo * 3,
             label = paste0("RE: ",
                            formatC(re_med_ratio, format = "f", digits = 0),
                            "x overestimate\n",
                            "Between FE: ",
                            formatC(fb_med_ratio, format = "f", digits = 0),
                            "x overestimate"),
             size = 2.8, fontface = "italic", color = "grey30",
             hjust = 1, vjust = 0) +
    scale_x_log10(limits = c(ax_lo, ax_hi),
                  labels = label_number(drop0trailing = TRUE)) +
    scale_y_log10(limits = c(ax_lo, ax_hi),
                  labels = label_number(drop0trailing = TRUE)) +
    scale_color_manual(values = TYPE_COLORS, labels = TYPE_LABELS,
                       name = "Parameter type") +
    scale_shape_manual(values = TYPE_SHAPES, labels = TYPE_LABELS,
                       name = "Parameter type") +
    labs(
      title = "(b) Simplified formulas",
      x = "Exact DER (Algorithm 1)",
      y = "Simplified DER"
    ) +
    theme_pub() +
    theme(legend.position = "none")

  ## Combined Figure 2 with shared legend
  p2a_leg <- p2a + theme(legend.position = "bottom")
  p2_final <- p2a_leg + p2b +
    plot_layout(ncol = 2, widths = c(1, 1))

  save_figure(p2_final, "figure2_formula_accuracy", width = 7.0, height = 4.0)
}


###############################################################################
## Section 3 : Figure 3 -- Coverage Comparison
###############################################################################
##
## Two-panel coverage plot comparing correction strategies:
##   Panel (a): Target parameter (within-cluster FE)
##              Shows that DER-based selective correction maintains nominal
##              coverage while naive deteriorates with increasing DEFF.
##   Panel (b): Non-target parameters (between-cluster FE + RE)
##              Demonstrates the blanket correction catastrophe:
##              inflating all parameters collapses RE coverage to ~20%.
##
## Paper Section 4.3, Figure 3.
##
###############################################################################

cat("--- Creating Figure 3: Coverage Comparison ---\n")

if (!"cov_table" %in% names(analysis)) {
  cat("  [SKIP] cov_table not found in analysis object.\n\n")
} else {

  cov <- analysis$cov_table

  ## Aggregate across all J and scenario conditions
  cov_agg <- cov %>%
    group_by(deff, strategy, type) %>%
    summarise(
      mean_coverage = mean(mean_coverage, na.rm = TRUE),
      se_coverage   = sd(mean_coverage, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  ## Strategy factor levels
  strategy_levels <- c("naive", "blanket", "DER-1.2", "DER-1.5")
  strategy_labels <- c("Naive", "Blanket",
                        "DER (tau = 1.2)", "DER (tau = 1.5)")

  cov_agg$strategy <- factor(cov_agg$strategy,
                              levels = strategy_levels,
                              labels = strategy_labels)
  cov_agg$deff_f <- factor(cov_agg$deff)

  ## Monte Carlo confidence band around nominal 90% coverage
  n_reps  <- 200
  mc_band <- 1.96 * sqrt(0.90 * 0.10 / n_reps)

  ## Panel (a): Within-cluster FE (target parameter)
  cov_fw <- cov_agg %>% filter(type == "fe_within")

  p3a <- ggplot(cov_fw, aes(x = deff_f, y = mean_coverage,
                             color = strategy, linetype = strategy,
                             shape = strategy, group = strategy)) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = 0.90 - mc_band, ymax = 0.90 + mc_band,
             fill = "grey50", alpha = 0.10) +
    geom_hline(yintercept = 0.90, linetype = "solid", color = "grey50",
               linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean_coverage - 1.96 * se_coverage,
                      ymax = mean_coverage + 1.96 * se_coverage),
                  width = 0.08, linewidth = 0.4) +
    scale_color_manual(values = STRATEGY_COLORS, name = "Strategy") +
    scale_linetype_manual(values = STRATEGY_LINETYPES, name = "Strategy") +
    scale_shape_manual(values = STRATEGY_SHAPES, name = "Strategy") +
    scale_y_continuous(limits = c(0.55, 1.0), breaks = seq(0.6, 1.0, 0.1),
                       labels = label_percent(accuracy = 1)) +
    labs(
      title = "(a) Target (within-cluster FE)",
      x = "Design effect (DEFF)",
      y = "Mean 90% CI coverage"
    ) +
    theme_pub() +
    theme(legend.position = "none")

  ## Panel (b): Non-target parameters (between-cluster FE + RE)
  cov_nontarget <- cov_agg %>%
    filter(type %in% c("fe_between", "re")) %>%
    group_by(deff_f, strategy) %>%
    summarise(
      mean_coverage = mean(mean_coverage, na.rm = TRUE),
      se_coverage   = sqrt(sum(se_coverage^2, na.rm = TRUE)) / n(),
      .groups = "drop"
    )

  p3b <- ggplot(cov_nontarget, aes(x = deff_f, y = mean_coverage,
                                    color = strategy, linetype = strategy,
                                    shape = strategy, group = strategy)) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = 0.90 - mc_band, ymax = 0.90 + mc_band,
             fill = "grey50", alpha = 0.10) +
    geom_hline(yintercept = 0.90, linetype = "solid", color = "grey50",
               linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = pmax(mean_coverage - 1.96 * se_coverage, 0),
                      ymax = pmin(mean_coverage + 1.96 * se_coverage, 1)),
                  width = 0.08, linewidth = 0.4) +
    ## Catastrophe annotation for blanket correction
    annotate("text",
             x = 2.8, y = 0.10,
             label = "~20% coverage",
             size = 2.8, fontface = "bold.italic",
             color = COL_RED) +
    annotate("segment",
             x = 2.55, xend = 2.25,
             y = 0.13, yend = 0.19,
             color = COL_RED, linewidth = 0.4,
             arrow = arrow(length = unit(0.06, "inches"), type = "closed")) +
    scale_color_manual(values = STRATEGY_COLORS, name = "Strategy") +
    scale_linetype_manual(values = STRATEGY_LINETYPES, name = "Strategy") +
    scale_shape_manual(values = STRATEGY_SHAPES, name = "Strategy") +
    scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2),
                       labels = label_percent(accuracy = 1)) +
    labs(
      title = "(b) Non-target (between FE + RE)",
      x = "Design effect (DEFF)",
      y = "Mean 90% CI coverage"
    ) +
    theme_pub() +
    theme(legend.position = "none")

  ## Combined Figure 3 with shared legend
  p3_final <- (p3a | p3b) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box.margin = ggplot2::margin(t = -6),
          legend.key.width = unit(1.2, "cm"))

  save_figure(p3_final, "figure3_coverage", width = 7.0, height = 4.5)
}


###############################################################################
## Section 4 : Summary
###############################################################################

cat("\n=== Simulation figures complete ===\n")
cat("Output directory:", file.path(here::here(), "output", "figures"), "\n")
