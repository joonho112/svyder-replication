# =============================================================================
# plotting_helpers.R: Plotting Utilities and Theme for Manuscript Figures
# =============================================================================
#
# Purpose : Defines the custom ggplot2 theme, color palettes for DER tier
#           classification, and shared figure-saving utilities used across
#           all manuscript figures.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Requires : ggplot2 (>= 3.4.0), scales, grid
#
# Contents:
#   1. theme_pub()         -- Publication-quality ggplot2 theme
#   2. Color palettes      -- Tier classification and parameter type colors
#   3. Shape encodings     -- Redundant encoding for colorblind accessibility
#   4. save_figure()       -- Save figures in PDF + PNG format
#
# =============================================================================


###############################################################################
## Section 1 : Project Paths
###############################################################################

## Set project root (use here::here() if available, otherwise detect)
if (requireNamespace("here", quietly = TRUE)) {
  PROJECT_ROOT <- here::here()
} else {
  PROJECT_ROOT <- getwd()
}


###############################################################################
## Section 2 : Publication-Quality ggplot2 Theme
###############################################################################

#' Custom ggplot2 theme for manuscript figures
#'
#' Based on theme_classic() with:
#'   - Subtle major grid lines for readability
#'   - Thin panel border
#'   - Bottom-positioned legend
#'   - Balanced margins for multi-panel layouts (patchwork)
#'
#' Designed for JSSAM print journal formatting (single-column width ~3.3 in,
#' double-column width ~6.9 in).
#'
#' @param base_size Numeric. Base font size (default: 11).
#'
#' @return A ggplot2 theme object.
theme_pub <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      ## Grid: subtle major lines for value reading
      panel.grid.major   = ggplot2::element_line(color = "grey92",
                                                  linewidth = 0.3),
      ## Panel border: thin frame
      panel.border       = ggplot2::element_rect(color = "grey40", fill = NA,
                                                  linewidth = 0.4),
      ## Axes
      axis.title         = ggplot2::element_text(size = base_size),
      axis.text          = ggplot2::element_text(size = base_size - 1,
                                                  color = "grey25"),
      axis.ticks         = ggplot2::element_line(color = "grey40",
                                                  linewidth = 0.3),
      ## Legend
      legend.background  = ggplot2::element_blank(),
      legend.key         = ggplot2::element_blank(),
      legend.title       = ggplot2::element_text(size = base_size - 1,
                                                  face = "bold"),
      legend.text        = ggplot2::element_text(size = base_size - 1.5),
      legend.position    = "bottom",
      legend.key.size    = grid::unit(0.8, "lines"),
      ## Facet strips
      strip.text         = ggplot2::element_text(size = base_size - 1,
                                                  face = "bold"),
      ## Titles
      plot.title         = ggplot2::element_text(size = base_size,
                                                  face = "bold",
                                                  margin = ggplot2::margin(b = 4)),
      plot.subtitle      = ggplot2::element_text(size = base_size - 2,
                                                  color = "grey40",
                                                  margin = ggplot2::margin(b = 4)),
      ## Margins: tight but balanced for patchwork composition
      plot.margin        = ggplot2::margin(6, 6, 4, 4)
    )
}


###############################################################################
## Section 3 : Color Palettes
###############################################################################

## --- Parameter type colors (colorblind-friendly) ---
## Used in Figures 2-4 for distinguishing parameter types.
## Selected from the ColorBrewer qualitative palette.

COL_FW  <- "#2166AC"   # blue   -- within-cluster fixed effects
COL_FB  <- "#4DAF4A"   # green  -- between-cluster fixed effects
COL_RE  <- "#FF7F00"   # orange -- random effects
COL_RED <- "#E41A1C"   # red    -- error/catastrophe indicator

TYPE_COLORS <- c(
  "fe_within"  = COL_FW,
  "fe_between" = COL_FB,
  "re"         = COL_RE
)

TYPE_LABELS <- c(
  "fe_within"  = "Within FE",
  "fe_between" = "Between FE",
  "re"         = "Random effects"
)


## --- DER tier classification colors ---
## Used in DER profile plots and classification summaries.

TIER_COLORS <- c(
  "Tier I"   = "#4DAF4A",   # green  -- Survey-Robust (no action needed)
  "Tier II"  = "#FF7F00",   # orange -- Survey-Sensitive (monitor)
  "Tier III" = "#E41A1C"    # red    -- Survey-Inflated (correct)
)

TIER_LABELS <- c(
  "Tier I"   = "Tier I: Survey-Robust",
  "Tier II"  = "Tier II: Survey-Sensitive",
  "Tier III" = "Tier III: Survey-Inflated"
)


## --- Correction strategy colors ---
## Used in Figure 3 (coverage comparison).

STRATEGY_COLORS <- c(
  "Naive"            = "black",
  "Blanket"          = COL_RED,
  "DER (tau = 1.2)"  = COL_FW,
  "DER (tau = 1.5)"  = COL_FB
)

STRATEGY_LINETYPES <- c(
  "Naive"            = "dashed",
  "Blanket"          = "dotted",
  "DER (tau = 1.2)"  = "solid",
  "DER (tau = 1.5)"  = "longdash"
)

STRATEGY_SHAPES <- c(
  "Naive"            = 1,
  "Blanket"          = 2,
  "DER (tau = 1.2)"  = 16,
  "DER (tau = 1.5)"  = 17
)


###############################################################################
## Section 4 : Shape Encodings
###############################################################################

## Redundant shape encoding ensures figures remain interpretable in
## grayscale printing (colorblind-friendly design principle).

TYPE_SHAPES <- c(
  "fe_within"  = 17,   # filled triangle
  "fe_between" = 15,   # filled square
  "re"         = 16    # filled circle
)


###############################################################################
## Section 5 : Figure Save Helpers
###############################################################################

#' Save a ggplot figure in PDF and PNG formats
#'
#' Saves the figure at publication resolution (300 DPI for PNG).
#' PDF is preferred for journal submission; PNG for supplementary materials
#' and web display.
#'
#' @param plot       A ggplot object.
#' @param name       Character. File name stem (without extension).
#' @param width      Numeric. Width in inches (default: 7.0).
#' @param height     Numeric. Height in inches (default: 4.5).
#' @param output_dir Character. Directory for output files. Defaults to
#'   PROJECT_ROOT/output/figures.
save_figure <- function(plot, name, width = 7.0, height = 4.5,
                        output_dir = file.path(PROJECT_ROOT, "output", "figures")) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  pdf_path <- file.path(output_dir, paste0(name, ".pdf"))
  png_path <- file.path(output_dir, paste0(name, ".png"))

  ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                  device = "pdf")
  ggplot2::ggsave(png_path, plot, width = width, height = height,
                  dpi = 300)

  cat(sprintf("  [SAVED] %s.pdf / .png (%g x %g in)\n", name, width, height))
}


#' Save a data frame as CSV and LaTeX table
#'
#' @param df         A data frame.
#' @param name       Character. File name stem (without extension).
#' @param caption    Character. Table caption for LaTeX.
#' @param output_dir Character. Directory for output files. Defaults to
#'   PROJECT_ROOT/output/tables.
#' @param ...        Additional arguments passed to xtable::xtable().
save_table <- function(df, name, caption = "",
                       output_dir = file.path(PROJECT_ROOT, "output", "tables"),
                       ...) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  csv_path <- file.path(output_dir, paste0(name, ".csv"))
  tex_path <- file.path(output_dir, paste0(name, ".tex"))

  write.csv(df, csv_path, row.names = FALSE)

  if (requireNamespace("xtable", quietly = TRUE)) {
    xt <- xtable::xtable(df, caption = caption, ...)
    print(xt, file = tex_path, include.rownames = FALSE,
          booktabs = TRUE, floating = FALSE)
  }

  cat(sprintf("  [SAVED] %s.csv / .tex\n", name))
}


###############################################################################
## Section 6 : DER Threshold Reference Line
###############################################################################

#' Add DER threshold reference lines to a ggplot
#'
#' Adds vertical or horizontal dashed lines at the DER thresholds (lo, hi)
#' with optional zone shading. Used in DER profile plots (Figure 4).
#'
#' @param lo  Numeric. Lower threshold (default 0.80).
#' @param hi  Numeric. Upper threshold (default 1.20).
#' @param direction Character. "vertical" or "horizontal" (default "vertical").
#'
#' @return A list of ggplot2 annotation layers.
add_der_thresholds <- function(lo = 0.80, hi = 1.20, direction = "vertical") {

  if (direction == "vertical") {
    list(
      ggplot2::geom_vline(xintercept = hi, linetype = "dashed",
                          color = "grey25", linewidth = 0.65),
      ggplot2::geom_vline(xintercept = lo, linetype = "dotted",
                          color = "grey50", linewidth = 0.4)
    )
  } else {
    list(
      ggplot2::geom_hline(yintercept = hi, linetype = "dashed",
                          color = "grey25", linewidth = 0.65),
      ggplot2::geom_hline(yintercept = lo, linetype = "dotted",
                          color = "grey50", linewidth = 0.4)
    )
  }
}
