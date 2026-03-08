# =============================================================================
# fig_02_application.R: Application Figure (Figure 4)
# =============================================================================
#
# Purpose : Generate Figure 4 from the NSECE empirical application
#           (Paper Section 5). Shows the DER profile by parameter type,
#           illustrating the "separation chasm" between survey-sensitive
#           and survey-robust parameters.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Track B : Uses pre-computed NSECE DER results from precomputed/.
#           No restricted data access required.
#
# Inputs  :
#   data/precomputed/application/phase3_classification.rds  -- DER classification
#   data/precomputed/application/phase3_sandwich.rds         -- Sandwich components
#   data/precomputed/application/phase3_summary.rds          -- Posterior summary
#
# Outputs :
#   output/figures/figure4_der_profile.pdf / .png
#
# =============================================================================


###############################################################################
## Section 0 : Setup
###############################################################################

library(ggplot2)
library(dplyr)
library(scales)

## Source plotting helpers
source(file.path(here::here(), "code", "helpers", "plotting_helpers.R"))

cat("=== Generating Application Figure (Track B) ===\n\n")


###############################################################################
## Section 1 : Load Pre-Computed Results
###############################################################################

cat("Loading pre-computed NSECE DER results ...\n")

PRECOMP_DIR <- file.path(here::here(), "data", "precomputed", "application")

## Load classification results (contains DER, tiers, param types)
class_path <- file.path(PRECOMP_DIR, "phase3_classification.rds")
sand_path  <- file.path(PRECOMP_DIR, "phase3_sandwich.rds")

if (!file.exists(class_path)) {
  stop("Pre-computed classification not found:\n  ", class_path,
       "\nPlease ensure data/precomputed/application/ files exist.\n",
       "See data/DATA_ACCESS.md for instructions.")
}

classification <- readRDS(class_path)
sandwich_data  <- if (file.exists(sand_path)) readRDS(sand_path) else NULL

cat("  Loaded phase3_classification.rds  [OK]\n")
if (!is.null(sandwich_data)) cat("  Loaded phase3_sandwich.rds  [OK]\n")
cat("\n")

## Extract needed components from classification object
der_values  <- classification$der          # named numeric [54]
param_types <- classification$param_types  # character [54]
param_names <- classification$param_names  # character [54]
deff_j      <- if (!is.null(sandwich_data)) sandwich_data$deff_j else
               classification$deff_j       # numeric [51]
B_j         <- if (!is.null(sandwich_data)) sandwich_data$B_j else
               classification$B_j          # numeric [51]

p <- classification$n_fe        # 3
J <- classification$n_states    # 51


###############################################################################
## Section 2 : Prepare Data for Plotting
###############################################################################

## Build the figure data frame
fig4_data <- data.frame(
  parameter = param_names,
  type      = param_types,
  der       = as.numeric(der_values),
  stringsAsFactors = FALSE
)

## Assign human-readable type labels for y-axis
## (Paper Figure 4 uses multi-line labels)
fig4_data$type_label <- factor(
  fig4_data$type,
  levels = c("re", "fe_between", "fe_within"),
  labels = c(
    paste0("Random effects\n(theta[1], ..., theta[", J, "])"),
    "Between-cluster FE\n(intercept, tiered_reim)",
    "Within-cluster FE\n(poverty_cwc)"
  )
)

## Summary statistics by parameter type
fig4_summary <- fig4_data %>%
  group_by(type, type_label) %>%
  summarise(
    med = median(der),
    lo  = min(der),
    hi  = max(der),
    n   = n(),
    .groups = "drop"
  )

cat("  DER summary by parameter type:\n")
for (i in seq_len(nrow(fig4_summary))) {
  cat(sprintf("    %-15s: n=%3d, median=%.4f, range=[%.4f, %.4f]\n",
              fig4_summary$type[i], fig4_summary$n[i],
              fig4_summary$med[i], fig4_summary$lo[i], fig4_summary$hi[i]))
}
cat("\n")


###############################################################################
## Section 3 : Figure 4 -- DER Profile / Separation Chasm
###############################################################################
##
## Horizontal strip plot showing the DER for all 54 parameters, grouped
## by type. The "separation chasm" refers to the orders-of-magnitude gap
## between the within-cluster FE DER (~DEFF) and the random effects DER
## (~0.01-0.10). This gap is the empirical manifestation of Theorem 2:
## hierarchical shrinkage dramatically reduces the sensitivity of random
## effects to survey design.
##
## Paper Section 5.2, Figure 4.
##
###############################################################################

cat("--- Creating Figure 4: DER Profile / Separation Chasm ---\n")

## DER threshold (Paper Section 3.2)
tau <- 1.20

## Mean Kish DEFF (for subtitle annotation)
mean_deff <- round(mean(deff_j), 2)

p4 <- ggplot(fig4_data, aes(x = der, y = type_label, color = type)) +

  ## Background zone shading
  annotate("rect",
           xmin = tau, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = COL_FW, alpha = 0.04) +
  annotate("rect",
           xmin = -Inf, xmax = tau, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.06) +

  ## Threshold line (tau = 1.20)
  geom_vline(xintercept = tau, linetype = "dashed", color = "grey25",
             linewidth = 0.65) +

  ## Random effects: jittered individual points (many overlapping)
  geom_point(data = ~ filter(., type == "re"),
             position = position_jitter(height = 0.18, width = 0, seed = 42),
             size = 1.4, alpha = 0.45, shape = 16) +

  ## Random effects: median indicator
  stat_summary(data = ~ filter(., type == "re"),
               fun = median, geom = "point",
               shape = "|", size = 5, color = "grey20") +

  ## Between-cluster FE: individual points (prominent)
  geom_point(data = ~ filter(., type == "fe_between"),
             size = 3.0, alpha = 0.85, shape = 15) +

  ## Within-cluster FE: single point (most prominent)
  geom_point(data = ~ filter(., type == "fe_within"),
             size = 4.0, alpha = 0.90, shape = 18) +

  ## Threshold label
  annotate("text", x = tau, y = 3.45,
           label = expression(paste(tau, " = 1.2")),
           size = 3.3, fontface = "bold", color = "grey25",
           hjust = -0.15) +

  ## Zone labels
  annotate("text", x = 0.006, y = 0.55,
           label = "No correction needed",
           size = 2.9, fontface = "italic", color = "grey45",
           hjust = 0) +
  annotate("text", x = 3.0, y = 0.55,
           label = "Correction needed",
           size = 2.9, fontface = "italic", color = COL_FW,
           hjust = 1) +

  ## Separation chasm annotation (double-headed arrow)
  annotate("segment",
           x = 0.10, xend = 1.10,
           y = 2.0, yend = 2.0,
           arrow = arrow(ends = "both", length = unit(0.08, "inches"),
                         type = "closed"),
           color = "grey50", linewidth = 0.4) +
  annotate("label", x = 0.33, y = 2.0,
           label = "Separation chasm",
           size = 2.7, fontface = "italic",
           fill = "white", alpha = 0.9,
           label.padding = unit(0.15, "lines"),
           color = "grey40") +

  ## Log scale x-axis with informative breaks
  scale_x_log10(
    breaks = c(0.01, 0.03, 0.1, 0.3, 1.2, 3.0),
    labels = c("0.01", "0.03", "0.1", "0.3", "1.2", "3.0"),
    limits = c(0.005, 4)
  ) +
  scale_color_manual(values = TYPE_COLORS, guide = "none") +

  labs(
    x = "Design effect ratio (DER)",
    y = NULL,
    title = "DER profile by parameter type",
    subtitle = paste0("NSECE 2019: J = ", J, " states, mean DEFF = ", mean_deff)
  ) +
  theme_pub() +
  theme(
    axis.text.y        = element_text(size = 10, face = "plain",
                                      lineheight = 1.1),
    panel.grid.major.y = element_blank(),
    plot.subtitle      = element_text(size = 9, color = "grey45",
                                      margin = margin(b = 8))
  )

save_figure(p4, "figure4_der_profile", width = 6.5, height = 4.0)


###############################################################################
## Section 4 : Summary
###############################################################################

cat("\n=== Application figure complete ===\n")
cat("Output directory:", file.path(here::here(), "output", "figures"), "\n")
cat("Files:\n")
cat("  figure4_der_profile.pdf / .png\n")
