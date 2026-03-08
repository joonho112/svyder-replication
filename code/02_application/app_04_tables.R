# =============================================================================
# app_04_tables.R: Tables for NSECE Application Section
# =============================================================================
#
# Purpose : Generate publication-quality tables for the NSECE application:
#           Table 5 (DER diagnostics for all parameters) and supplementary
#           tables (full 54-row DER profile, tau sensitivity, MCMC diagnostics,
#           summary statistics).
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Section : Section 5 (Application: NSECE 2019) and Online Supplementary
#           Materials D (OSM-D)
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Track   : A (Full Replication) / B (Pre-computed)
# Inputs  : Track A: output/application/phase3_*.rds  (from app_01--app_03)
#           Track B: data/precomputed/application/phase3_*.rds (pre-computed)
#           Required: phase3_sandwich.rds, phase3_classification.rds
#           Optional: phase3_summary.rds, phase3_data.rds, phase3_theorem_verification.rds
# Outputs : output/tables/tab5_nsece_der.tex             (Table 5, main text)
#           output/tables/tabSM_full_der.tex             (full 54-row DER)
#           output/tables/tabSM_tau_sensitivity.tex      (tau sensitivity)
#           output/tables/tabSM_mcmc_diagnostics.tex     (MCMC diagnostics)
#           output/tables/tab_summary_stats.tex          (summary statistics)
#
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------

# Track B uses the same input data as Track A for table generation.
# The key difference is in the upstream scripts (app_01 to app_03), which
# either compute from scratch (Track A) or load pre-computed results (Track B).
# By this point, the output/application/ directory contains all needed inputs
# regardless of which track was used upstream.

REPLICATION_ROOT <- here::here()

# Input directory: check output/application first (Track A), then fall back
# to data/precomputed/application (Track B)
OUTPUT_APP_DIR  <- file.path(REPLICATION_ROOT, "output", "application")
PRECOMP_APP_DIR <- file.path(REPLICATION_ROOT, "data", "precomputed", "application")

if (dir.exists(OUTPUT_APP_DIR) &&
    file.exists(file.path(OUTPUT_APP_DIR, "phase3_sandwich.rds"))) {
  INPUT_DIR <- OUTPUT_APP_DIR
  cat("  Using Track A output: ", INPUT_DIR, "\n")
} else if (dir.exists(PRECOMP_APP_DIR) &&
           file.exists(file.path(PRECOMP_APP_DIR, "phase3_sandwich.rds"))) {
  INPUT_DIR <- PRECOMP_APP_DIR
  cat("  Using Track B pre-computed: ", INPUT_DIR, "\n")
} else {
  stop("No application results found in:\n  ", OUTPUT_APP_DIR,
       "\n  ", PRECOMP_APP_DIR,
       "\nRun app_03_der_analysis.R first (Track A) or ensure ",
       "precomputed data exists (Track B).")
}

# Output directory for tables
TABLE_DIR <- file.path(REPLICATION_ROOT, "output", "tables")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# Constants
P_FIXED  <- 3L
J_STATES <- 51L
D_TOTAL  <- P_FIXED + J_STATES


# =============================================================================
# 1. Load Input Data
# =============================================================================

cat("=== Table Generation for NSECE Application ===\n\n")
cat("Step 1: Loading input data ...\n")

sandwich <- readRDS(file.path(INPUT_DIR, "phase3_sandwich.rds"))
classify <- readRDS(file.path(INPUT_DIR, "phase3_classification.rds"))

# Load posterior summary if available (for MCMC diagnostics table)
summary_path <- file.path(INPUT_DIR, "phase3_summary.rds")
has_summary <- file.exists(summary_path)
if (has_summary) {
  phase3_summary <- readRDS(summary_path)
}

# Load Stan data if available (for summary statistics)
data_path <- file.path(INPUT_DIR, "phase3_data.rds")
has_data <- file.exists(data_path)
if (has_data) {
  stan_data <- readRDS(data_path)
}

# Extract key quantities
der_all     <- sandwich$der
param_names <- sandwich$param_names
param_types <- sandwich$param_types
beta_hat    <- sandwich$beta_hat
deff_j      <- sandwich$deff_j
B_j         <- sandwich$B_j

# Credible intervals from classification
ci_naive     <- classify$ci_naive
ci_selective <- classify$ci_selective

# Covariate names for display
covariate_names <- c("intercept", "poverty_cwc", "tiered_reim")

cat(sprintf("  Loaded DER profile: %d parameters\n", length(der_all)))
cat(sprintf("  Classification: %d flagged of %d total\n",
            classify$n_flagged, classify$n_params))


# =============================================================================
# 2. Table 5: DER Diagnostics for Fixed-Effect Parameters (Main Text)
# =============================================================================
#
# This is the primary results table in Section 5.2 (Table 5 in the paper).
# Reports DER, R_k, tier, action, and credible interval comparison for the
# 3 fixed effects plus a summary row for the 51 random effects.

cat("\nStep 2: Generating Table 5 (main text DER diagnostics) ...\n")

# Extract R_k values if available from theorem verification
theorem_path <- file.path(INPUT_DIR, "phase3_theorem_verification.rds")
if (file.exists(theorem_path)) {
  theorem_v <- readRDS(theorem_path)
  R_k <- theorem_v$theorem1_state$R_k
} else {
  # Compute R_k from H_obs blocks if theorem verification not available
  H <- sandwich$H_obs
  H_bb <- H[1:P_FIXED, 1:P_FIXED]
  H_bt <- H[1:P_FIXED, (P_FIXED + 1):D_TOTAL]
  H_tt <- H[(P_FIXED + 1):D_TOTAL, (P_FIXED + 1):D_TOTAL]
  H_tt_inv <- solve(H_tt)

  R_k <- vapply(seq_len(P_FIXED), function(k) {
    H_bt_k <- H_bt[k, ]
    num <- as.numeric(t(H_bt_k) %*% H_tt_inv %*% H_bt_k)
    num / H_bb[k, k]
  }, numeric(1))
  names(R_k) <- covariate_names
}

# Fixed-effect DER values
fe_der <- der_all[1:P_FIXED]

# Width change percentages
width_change_pct <- (ci_selective$width[1:P_FIXED] -
                       ci_naive$width[1:P_FIXED]) /
  ci_naive$width[1:P_FIXED] * 100

# Random-effect DER summary
re_der <- der_all[(P_FIXED + 1):D_TOTAL]
re_range_str <- paste0("$[", formatC(min(re_der), format = "f", digits = 3),
                        ",\\; ", formatC(max(re_der), format = "f", digits = 3),
                        "]$")

# Format helper for credible intervals
fmt_ci <- function(lower, upper) {
  paste0("$[", formatC(lower, format = "f", digits = 3), ",\\; ",
         formatC(upper, format = "f", digits = 3), "]$")
}

# Pre-format numeric values
der_str <- formatC(fe_der, format = "f", digits = 3)
Rk_str  <- formatC(R_k, format = "f", digits = 3)
ci_n_str <- vapply(1:P_FIXED, function(i)
  fmt_ci(ci_naive$lower[i], ci_naive$upper[i]), "")
ci_s_str <- vapply(1:P_FIXED, function(i)
  fmt_ci(ci_selective$lower[i], ci_selective$upper[i]), "")

# Build LaTeX table rows
# Row 1: intercept (beta_0, between-state, Tier I-b, retain)
row1 <- paste0(
  "$\\beta_0$ (intercept) & Between & ", der_str[1], " & ",
  Rk_str[1], " & I-b & Retain & ", ci_n_str[1], " & ",
  ci_s_str[1], " & $0\\%$ \\\\[2pt]"
)

# Row 2: poverty (beta_1, within-state, Tier I-a, CORRECT) -- bolded
wchg_str <- sprintf("%+.1f\\%%", width_change_pct[2])
row2 <- paste0(
  "$\\boldsymbol{\\beta_1}$ \\textbf{(poverty)} & \\textbf{Within} & \\textbf{",
  der_str[2], "} & \\textbf{", Rk_str[2],
  "} & \\textbf{I-a} & \\textbf{Correct} & $\\mathbf{",
  gsub("\\$", "", ci_n_str[2]),
  "}$ & $\\mathbf{",
  gsub("\\$", "", ci_s_str[2]),
  "}$ & $\\mathbf{", wchg_str, "}$ \\\\[2pt]"
)

# Row 3: tiered reimbursement (beta_2, between-state, Tier I-b, retain)
row3 <- paste0(
  "$\\beta_2$ (tiered reim.) & Between & ", der_str[3], " & ",
  Rk_str[3], " & I-b & Retain & ", ci_n_str[3], " & ",
  ci_s_str[3], " & $0\\%$ \\\\"
)

# Row 4: Random effects summary
row_re <- paste0(
  "$\\theta_j$ (51 state REs) & RE & ", re_range_str,
  " & --- & II & Retain & \\multicolumn{2}{c}{\\textit{all unchanged}} & $0\\%$ \\\\"
)

# Assemble complete LaTeX table
tab5_lines <- c(
  "\\begin{table}[t]",
  "\\centering",
  paste0("\\caption{Design effect ratio diagnostics for fixed-effect parameters ",
         "in the NSECE application. The selective correction applies sandwich ",
         "variance adjustment only to parameters with $\\mathrm{DER}_p > \\tau ",
         "= 1.2$ (Tier~I-a), leaving protected parameters unchanged.}"),
  "\\label{tab:nsece_der}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{@{}llcccl cc r@{}}",
  "\\toprule",
  " &  &  &  &  &  & \\multicolumn{2}{c}{90\\% Credible Interval} & \\\\",
  "\\cmidrule(lr){7-8}",
  paste0("Parameter & Type & DER & $R_k$ & Tier & Action & ",
         "Naive & Selective & $\\Delta w$ \\\\"),
  "\\midrule",
  row1,
  row2,
  row3,
  "\\midrule",
  row_re,
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "",
  "\\vspace{4pt}",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  paste0("\\textit{Notes.} $\\mathrm{DER}_p = \\mathrm{diag}(",
         "\\mathbf{V}_{\\mathrm{sand}})_p / \\mathrm{diag}(",
         "\\boldsymbol{\\Sigma}_{\\mathrm{MCMC}})_p$. ",
         "$R_k$ quantifies the share of between-cluster information for ",
         "covariate~$k$ (high $R_k \\Rightarrow$ protected). ",
         "Tier~I-a: within-cluster covariates with $\\mathrm{DER} > \\tau$ ",
         "(survey-dominated); Tier~I-b: between-cluster covariates (protected ",
         "by confounding with cluster means); Tier~II: random effects ",
         "(protected by hierarchical shrinkage). ",
         "$\\Delta w = (w_{\\mathrm{selective}} - w_{\\mathrm{naive}}) / ",
         "w_{\\mathrm{naive}} \\times 100\\%$. ",
         "Only 1 of 54 parameters (1.9\\%) requires correction under ",
         "$\\tau = 1.2$."),
  "\\end{minipage}",
  "",
  "\\end{table}"
)

tab5_path <- file.path(TABLE_DIR, "tab5_nsece_der.tex")
writeLines(tab5_lines, tab5_path)
cat(sprintf("  Saved: %s\n", tab5_path))


# =============================================================================
# 3. Supplementary Table: Full 54-Row DER Profile (OSM-D)
# =============================================================================
#
# Reports DER, DER (Laplace), DEFF_j, B_j, tier, and action for all 54
# parameters. This table appears in Online Supplementary Material D.

cat("\nStep 3: Generating full 54-row DER table (OSM-D) ...\n")

der_laplace <- sandwich$der_laplace
classification_df <- classify$classification

# Build LaTeX table rows
full_rows <- character(D_TOTAL)

for (k in seq_len(D_TOTAL)) {
  nm  <- param_names[k]
  typ <- param_types[k]
  dr  <- der_all[k]
  dl  <- der_laplace[k]
  tier <- classification_df$tier[k]
  act  <- classification_df$action[k]

  # Format parameter name for LaTeX
  if (k <= P_FIXED) {
    if (k == 1) latex_nm <- "$\\beta_0$ (intercept)"
    else if (k == 2) latex_nm <- "$\\beta_1$ (poverty)"
    else latex_nm <- "$\\beta_2$ (tiered reim.)"

    row_str <- sprintf(
      "%s & %s & %.4f & %.4f & --- & --- & %s & %s",
      latex_nm, typ, dr, dl, tier, act
    )
  } else {
    j_idx <- k - P_FIXED
    latex_nm <- sprintf("$\\theta_{%d}$", j_idx)

    row_str <- sprintf(
      "%s & re & %.4f & %.4f & %.4f & %.4f & %s & %s",
      latex_nm, dr, dl, deff_j[j_idx], B_j[j_idx], tier, act
    )
  }

  # Bold if flagged
  if (classification_df$flagged[k]) {
    row_str <- paste0("\\textbf{", row_str, "}")
  }

  full_rows[k] <- paste0(row_str, " \\\\")
}

# Add midrule between FE and RE sections
full_table_lines <- c(
  "\\begin{longtable}{@{}llrrrr ll@{}}",
  "\\caption{Full DER profile for all 54 parameters in the NSECE application.}",
  "\\label{tab:full_der} \\\\",
  "\\toprule",
  "Parameter & Type & DER & DER$_L$ & DEFF$_j$ & $B_j$ & Tier & Action \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\multicolumn{8}{l}{\\textit{Continued from previous page}} \\\\",
  "\\toprule",
  "Parameter & Type & DER & DER$_L$ & DEFF$_j$ & $B_j$ & Tier & Action \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\multicolumn{8}{r}{\\textit{Continued on next page}} \\\\",
  "\\endfoot",
  "\\bottomrule",
  "\\endlastfoot",
  "",
  "\\multicolumn{8}{l}{\\textit{Fixed Effects}} \\\\[2pt]",
  full_rows[1:P_FIXED],
  "\\midrule",
  "\\multicolumn{8}{l}{\\textit{Random Effects (State)}} \\\\[2pt]",
  full_rows[(P_FIXED + 1):D_TOTAL],
  "",
  "\\end{longtable}"
)

full_der_path <- file.path(TABLE_DIR, "tabSM_full_der.tex")
writeLines(full_table_lines, full_der_path)
cat(sprintf("  Saved: %s (%d rows)\n", full_der_path, D_TOTAL))


# =============================================================================
# 4. Supplementary Table: Tau Sensitivity Analysis (OSM-D)
# =============================================================================
#
# Section 5.4: "The DER classification is robust to the choice of tau:
# across six thresholds from 0.80 to 2.00, only at tau = 0.80 does one
# additional parameter (a single random effect) become flagged."

cat("\nStep 4: Generating tau sensitivity table ...\n")

tau_values <- c(0.80, 1.00, 1.10, 1.20, 1.50, 2.00)

tau_rows <- character(length(tau_values))

for (t_idx in seq_along(tau_values)) {
  tau_val <- tau_values[t_idx]

  # Count flagged at this threshold
  n_flag_fe <- sum(der_all[1:P_FIXED] > tau_val)
  n_flag_re <- sum(der_all[(P_FIXED + 1):D_TOTAL] > tau_val)
  n_flag_total <- n_flag_fe + n_flag_re

  # Identify flagged parameter names
  flagged_names <- param_names[der_all > tau_val]
  if (length(flagged_names) == 0) {
    flagged_str <- "---"
  } else {
    # Format for LaTeX
    flagged_str <- paste(gsub("\\[", "$_{", gsub("\\]", "}$", flagged_names)),
                         collapse = ", ")
  }

  # Percentage
  pct <- sprintf("%.1f\\%%", 100 * n_flag_total / D_TOTAL)

  # Mark recommended threshold
  marker <- ifelse(tau_val == 1.2, "\\textbf{", "")
  close_marker <- ifelse(tau_val == 1.2, "}", "")

  tau_rows[t_idx] <- sprintf(
    "%s%.2f%s & %s%d%s & %s%d%s & %s%d%s & %s%s%s & %s%s%s \\\\",
    marker, tau_val, close_marker,
    marker, n_flag_fe, close_marker,
    marker, n_flag_re, close_marker,
    marker, n_flag_total, close_marker,
    marker, pct, close_marker,
    marker, flagged_str, close_marker
  )
}

tau_table_lines <- c(
  "\\begin{table}[t]",
  "\\centering",
  paste0("\\caption{Sensitivity of DER classification to threshold ",
         "$\\tau$. The recommended threshold $\\tau = 1.2$ is shown ",
         "in bold. The poverty coefficient ($\\beta_1$) is flagged at ",
         "every threshold; one random effect is additionally flagged ",
         "only at $\\tau = 0.80$.}"),
  "\\label{tab:tau_sensitivity}",
  "\\begin{tabular}{@{}cccccl@{}}",
  "\\toprule",
  "$\\tau$ & FE Flagged & RE Flagged & Total & \\% Flagged & Flagged Parameters \\\\",
  "\\midrule",
  tau_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

tau_path <- file.path(TABLE_DIR, "tabSM_tau_sensitivity.tex")
writeLines(tau_table_lines, tau_path)
cat(sprintf("  Saved: %s\n", tau_path))


# =============================================================================
# 5. Supplementary Table: MCMC Diagnostics (OSM-D)
# =============================================================================
#
# Reports convergence diagnostics for key parameters: posterior mean, SD,
# Rhat, ESS_bulk, ESS_tail. Corresponds to Section 5.1 diagnostics.

cat("\nStep 5: Generating MCMC diagnostics table ...\n")

if (has_summary) {

  summ_df <- phase3_summary$summary
  diag    <- phase3_summary$diagnostics

  # Select key parameters for display
  key_params <- c("beta[1]", "beta[2]", "beta[3]", "sigma_theta")

  diag_rows <- character(length(key_params))

  for (k in seq_along(key_params)) {
    pname <- key_params[k]
    row <- summ_df[summ_df$variable == pname, ]

    if (nrow(row) == 0) next

    # Format parameter name for LaTeX
    if (pname == "beta[1]") latex_nm <- "$\\beta_0$ (intercept)"
    else if (pname == "beta[2]") latex_nm <- "$\\beta_1$ (poverty)"
    else if (pname == "beta[3]") latex_nm <- "$\\beta_2$ (tiered reim.)"
    else latex_nm <- "$\\sigma_\\theta$"

    diag_rows[k] <- sprintf(
      "%s & %.3f & %.3f & %.4f & %.0f & %.0f \\\\",
      latex_nm, row$mean, row$sd, row$rhat,
      row$ess_bulk, row$ess_tail
    )
  }

  # Add summary row for random effects
  theta_rows <- summ_df[grepl("^theta\\[", summ_df$variable), ]
  if (nrow(theta_rows) > 0) {
    re_summary_row <- sprintf(
      "$\\theta_j$ (51 states) & [%.3f, %.3f] & [%.3f, %.3f] & [%.4f, %.4f] & [%.0f, %.0f] & [%.0f, %.0f] \\\\",
      min(theta_rows$mean), max(theta_rows$mean),
      min(theta_rows$sd), max(theta_rows$sd),
      min(theta_rows$rhat), max(theta_rows$rhat),
      min(theta_rows$ess_bulk), max(theta_rows$ess_bulk),
      min(theta_rows$ess_tail), max(theta_rows$ess_tail)
    )
  } else {
    re_summary_row <- "$\\theta_j$ (51 states) & --- & --- & --- & --- & --- \\\\"
  }

  diag_table_lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    paste0("\\caption{MCMC convergence diagnostics for the NSECE ",
           "hierarchical logistic regression. Four chains of 2{,}000 ",
           "post-warmup iterations (2{,}000 warmup, ",
           "\\texttt{adapt\\_delta} $= 0.95$). Zero divergent ",
           "transitions across all chains.}"),
    "\\label{tab:mcmc_diagnostics}",
    "\\begin{tabular}{@{}lrrrrr@{}}",
    "\\toprule",
    "Parameter & Mean & SD & $\\hat{R}$ & ESS$_{\\mathrm{bulk}}$ & ESS$_{\\mathrm{tail}}$ \\\\",
    "\\midrule",
    diag_rows,
    "\\midrule",
    re_summary_row,
    "\\bottomrule",
    "\\end{tabular}",
    "",
    "\\vspace{4pt}",
    "\\begin{minipage}{\\textwidth}",
    "\\footnotesize",
    paste0("\\textit{Notes.} Ranges shown for random effects. ",
           "All parameters satisfy $\\hat{R} < 1.01$ and ",
           "ESS$_{\\mathrm{bulk}} > 400$."),
    "\\end{minipage}",
    "",
    "\\end{table}"
  )

  diag_path <- file.path(TABLE_DIR, "tabSM_mcmc_diagnostics.tex")
  writeLines(diag_table_lines, diag_path)
  cat(sprintf("  Saved: %s\n", diag_path))

} else {
  cat("  Skipped (phase3_summary.rds not available).\n")
}


# =============================================================================
# 6. Summary Statistics Table
# =============================================================================
#
# Descriptive statistics for the NSECE dataset: sample sizes, outcome
# prevalence, covariate distributions, survey design characteristics.

cat("\nStep 6: Generating summary statistics table ...\n")

if (has_data) {

  N <- stan_data$N
  J <- stan_data$J

  state_sizes <- as.integer(table(stan_data$group))
  prevalence  <- mean(stan_data$y)

  # Covariate summaries
  poverty_cwc_sd <- sd(stan_data$X[, 2])
  tiered_prop    <- mean(stan_data$X[, 3])

  # Weight summaries
  w_mean <- mean(stan_data$w)
  w_cv   <- sd(stan_data$w) / mean(stan_data$w)

  summary_table_lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    paste0("\\caption{Summary statistics for the NSECE 2019 application. ",
           "$N = 6{,}785$ center-based childcare providers across ",
           "$J = 51$ states (including DC).}"),
    "\\label{tab:nsece_summary}",
    "\\begin{tabular}{@{}lr@{}}",
    "\\toprule",
    "Characteristic & Value \\\\",
    "\\midrule",
    sprintf("Observations ($N$) & %s \\\\", formatC(N, big.mark = ",")),
    sprintf("States ($J$) & %d \\\\", J),
    sprintf("State sizes: min / median / max & %d / %d / %d \\\\",
            min(state_sizes), as.integer(median(state_sizes)), max(state_sizes)),
    "\\midrule",
    sprintf("IT enrollment prevalence & %.1f\\%% \\\\", prevalence * 100),
    sprintf("Poverty (CWC): SD & %.4f \\\\", poverty_cwc_sd),
    sprintf("Tiered reimbursement: proportion & %.1f\\%% \\\\", tiered_prop * 100),
    "\\midrule",
    sprintf("Normalized weight mean & %.4f \\\\", w_mean),
    sprintf("Normalized weight CV & %.4f \\\\", w_cv),
    sprintf("Global Kish DEFF & %.2f \\\\", sandwich$kish_deff_global %||%
              (N * sum(stan_data$w^2) / sum(stan_data$w)^2)),
    sprintf("Per-state DEFF: mean [min, max] & %.2f [%.2f, %.2f] \\\\",
            mean(deff_j), min(deff_j), max(deff_j)),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  summary_path <- file.path(TABLE_DIR, "tab_summary_stats.tex")
  writeLines(summary_table_lines, summary_path)
  cat(sprintf("  Saved: %s\n", summary_path))

} else {
  cat("  Skipped (phase3_data.rds not available).\n")
}


# =============================================================================
# 7. Verification
# =============================================================================

cat("\n=== Verification ===\n\n")

# Check all output files exist and have non-zero size
output_files <- c(
  tab5_path,
  full_der_path,
  tau_path
)

# Add conditional files
if (has_summary) output_files <- c(output_files, diag_path)
if (has_data) output_files <- c(output_files, summary_path)

all_ok <- TRUE
for (f in output_files) {
  info <- file.info(f)
  exists_ok <- !is.na(info$size) && info$size > 0
  cat(sprintf("  [%s] %s  (%.1f KB)\n",
              ifelse(exists_ok, "OK", "FAIL"),
              basename(f),
              ifelse(!is.na(info$size), info$size / 1024, 0)))
  if (!exists_ok) all_ok <- FALSE
}

# Data consistency checks
cat("\n  Data consistency:\n")
cat(sprintf("    FE DER: [%.4f, %.4f, %.4f]\n",
            fe_der[1], fe_der[2], fe_der[3]))
cat(sprintf("    R_k:    [%.3f, %.3f, %.3f]\n", R_k[1], R_k[2], R_k[3]))
cat(sprintf("    RE DER: range [%.4f, %.4f], mean = %.4f\n",
            min(re_der), max(re_der), mean(re_der)))
cat(sprintf("    Flagged: %d of %d at tau = %.1f\n",
            classify$n_flagged, classify$n_params, classify$tau_threshold))

# Key paper result: beta_1 CI width change
width_change_beta1 <- (ci_selective$width[2] - ci_naive$width[2]) /
  ci_naive$width[2] * 100
cat(sprintf("    beta_1 CI width change: +%.1f%% (paper reports +62.6%%)\n",
            width_change_beta1))


# =============================================================================
# Summary
# =============================================================================

cat(sprintf("\n=== Table Generation Complete ===\n"))
cat(sprintf("  Output directory: %s\n", TABLE_DIR))
cat(sprintf("  Tables generated: %d files\n", length(output_files)))
cat(sprintf("  All files OK: %s\n", ifelse(all_ok, "YES", "NO")))
cat("=== Proceed to 03_figures/ for figure generation ===\n")
