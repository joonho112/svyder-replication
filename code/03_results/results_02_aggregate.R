# =============================================================================
# results_02_aggregate.R
# Author : JoonHo Lee (jlee296@ua.edu)
# -----------------------------------------------------------------------------
# Reproduce the simulation and application figures and tables of the paper from
# the shipped results. This is the reporting stage: it fits no models and needs
# no restricted data.
#
# Reads
#   output/simulation/<cell>/rep_XXXX.rds              the 4,800 simulation reps
#                                                      (24 design cells x 200;
#                                                      built by code/01_simulation/)
#   data/precomputed/application/nsece_v2_refit.rds    the NSECE v2 model fit
#   data/precomputed/application/nsece_supp.rds        NSECE supplement (weight-
#                                                      convention sensitivity and
#                                                      structural factors R_k)
# Writes
#   output/figures/fig2_der_separation.pdf             Figure 2
#   output/figures/fig3_coverage.pdf                   Figure 3
#   output/figures/fig4_nsece_dualtarget.pdf           Figure 4
#     (Figure 1 is the theory panel, drawn by results_01_theory.R)
#   output/tables/tab2_sim_coverage.tex                Table 2
#   output/tables/tab3_nsece_der.tex                   Table 3
#   output/tables/tab4_nsece_2x2.tex                   Table 4
#   output/tables/tab5_c0_sweep.tex                    Table 5
#   output/tables/simulation_cell_summary.{rds,csv}    per-cell summary
#
# Gate-failing replications are excluded from every summary but counted and
# reported. Figures use a colorblind-safe Okabe-Ito palette; tables use booktabs.
#
# Usage : Rscript code/03_results/results_02_aggregate.R
#         (run from the package root, or open svyder-replication.Rproj first).
# =============================================================================

t_run0 <- proc.time()["elapsed"]

# =============================================================================
# Section 0: Root detection, paths, palettes, formatting helpers
# =============================================================================

ROOT    <- tryCatch(here::here(), error = function(e) getwd())
SIM_DIR <- file.path(ROOT, "output", "simulation")
FIG_DIR <- file.path(ROOT, "output", "figures")
TAB_DIR <- file.path(ROOT, "output", "tables")
PRECOMP <- file.path(ROOT, "data", "precomputed")
for (d in c(FIG_DIR, TAB_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
message("Package root: ", ROOT)

# ggplot2 is allowed for figures if it loads; this run requires it.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is not available; install it or port the figure code to base ",
       "graphics.", call. = FALSE)
}
library(ggplot2)

# Okabe-Ito palette (colorblind-safe)
OI <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00",
        purple = "#CC79A7", black = "#000000")

C0_DEFAULT <- 1.2
C0_GRID    <- c(0.8, 1.0, 1.2, 1.5, 2.0)   # tab5 sweep
CI_LEVEL   <- 0.90

fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d)
pctf <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
lint <- function(n) format(as.integer(round(n)), big.mark = "{,}",
                           scientific = FALSE, trim = TRUE)

# =============================================================================
# Section 1: Simulation aggregation core (Part A)
# =============================================================================

message("\n[1/6] Loading ", length(list.files(SIM_DIR)), " simulation cells ...")

cells <- sort(list.dirs(SIM_DIR, recursive = FALSE, full.names = FALSE))
stopifnot(length(cells) == 24L)

parse_cell <- function(cell_id) {
  parts <- strsplit(cell_id, "_")[[1]]
  list(J = as.integer(sub("^J", "", parts[1])),
       nbar_j = as.integer(sub("^N", "", parts[2])),
       informativeness = c(GNON = "none", GMOD = "moderate",
                           GSTR = "strong")[parts[3]],
       structure = c(SG = "sampled_groups", PW = "psu_within_groups")[parts[4]])
}

# Per-rep coverage summary for one strategy vector (named 0/1, length 3 + J)
.class_cov <- function(cv, J) {
  c(w  = unname(cv[["beta[2]"]]),
    b  = mean(c(cv[["beta[1]"]], cv[["beta[3]"]])),
    re = mean(cv[3 + seq_len(J)]))
}

rep_rows   <- vector("list", length(cells))   # per-rep long records
cell_extra <- vector("list", length(cells))   # per-cell param-level DER means
file_count <- integer(length(cells))
names(cell_extra) <- cells

first_checked <- FALSE

for (ci in seq_along(cells)) {
  cell_id <- cells[ci]
  fac <- parse_cell(cell_id)
  files <- sort(list.files(file.path(SIM_DIR, cell_id),
                           pattern = "^rep_\\d{4}\\.rds$", full.names = TRUE))
  file_count[ci] <- length(files)

  J <- fac$J
  der_sum   <- NULL   # running sum of der_group over gate-pass reps
  der_n     <- 0L
  re_pool_g <- NULL   # pooled RE DERs (group target), gate-pass
  re_pool_p <- NULL   # pooled RE DERs (psu target), gate-pass, PW arm

  rows <- vector("list", length(files))
  for (fi in seq_along(files)) {
    r <- readRDS(files[fi])

    if (!first_checked) {
      # One-time structural assertions on names/order
      pn <- c(sprintf("beta[%d]", 1:3), sprintf("theta[%d]", seq_len(J)))
      stopifnot(identical(names(r$der_group), pn),
                identical(names(r$coverage$group$naive), pn),
                identical(names(r$ci_width$group$naive), pn),
                identical(r$param_types[1:3],
                          c("fe_between", "fe_within", "fe_between")),
                all(r$param_types[-(1:3)] == "re"),
                r$tau == C0_DEFAULT, r$ci_level == CI_LEVEL)
      first_checked <- TRUE
    }

    dg  <- as.numeric(r$der_group)          # length 3 + J
    idx_re <- 3 + seq_len(J)
    idx_nt <- setdiff(seq_len(3 + J), 2L)   # non-target: fe_between + re
    gate_pass <- isTRUE(r$gate_pass)

    pw <- !isTRUE(r$targets_coincide)
    dp <- if (pw) as.numeric(r$der_psu) else rep(NA_real_, 3 + J)

    # Coverage per class x strategy, both targets
    cg <- lapply(r$coverage$group, .class_cov, J = J)
    cp <- if (pw) lapply(r$coverage$psu, .class_cov, J = J) else
      list(naive = c(w = NA_real_, b = NA_real_, re = NA_real_),
           selective = c(w = NA_real_, b = NA_real_, re = NA_real_),
           blanket = c(w = NA_real_, b = NA_real_, re = NA_real_))

    # Flagged-RE coverage under the design-PSU target (PW arm only):
    # the theta components of the flagged set, naive vs selective.
    n_fre <- 0L; fre_nv <- 0; fre_sel <- 0
    if (pw && length(r$flagged$psu) > 0) {
      fre <- r$flagged$psu[r$flagged$psu > 3]
      n_fre <- length(fre)
      if (n_fre > 0) {
        fre_nv  <- sum(r$coverage$psu$naive[fre])
        fre_sel <- sum(r$coverage$psu$selective[fre])
      }
    }

    max_nt_g <- max(dg[idx_nt])
    rows[[fi]] <- data.frame(
      cell_id = cell_id, rep_id = r$rep_id, gate_pass = gate_pass,
      J = J, nbar_j = fac$nbar_j, informativeness = fac$informativeness,
      structure = fac$structure,
      deff = r$diagnostics$deff_kish_global, cv_w = r$diagnostics$cv_w,
      N = r$diagnostics$N,
      sigma_theta_hat = r$sigma_theta_hat,
      timing_total = r$timing$total,
      # DER by class, group target
      der_w_g = dg[2], der_b_g = max(dg[c(1, 3)]),
      der_int_g = dg[1], der_bc_g = dg[3],
      der_re_med_g = stats::median(dg[idx_re]),
      max_nt_g = max_nt_g,
      # non-target (protected) flag decisions unchanged across c0 in [1.1, 1.6]
      stab_nt_g = sum(dg[idx_nt] > 1.1) == sum(dg[idx_nt] > 1.6),
      # DER by class, psu target (PW arm)
      der_w_p = dp[2],
      der_b_p = if (pw) max(dp[c(1, 3)]) else NA_real_,
      der_re_med_p = if (pw) stats::median(dp[idx_re]) else NA_real_,
      # flagged counts (stored) and threshold-derived counts, group target
      nflag_g = unname(r$n_flagged["group"]),
      nflag_p = unname(r$n_flagged["psu"]),
      nflag_g_c11 = sum(dg > 1.1), nflag_g_c16 = sum(dg > 1.6),
      # coverage: group target
      cov_g_nv_w = cg$naive["w"],  cov_g_nv_b = cg$naive["b"],  cov_g_nv_re = cg$naive["re"],
      cov_g_sel_w = cg$selective["w"], cov_g_sel_b = cg$selective["b"], cov_g_sel_re = cg$selective["re"],
      cov_g_blk_w = cg$blanket["w"], cov_g_blk_b = cg$blanket["b"], cov_g_blk_re = cg$blanket["re"],
      # coverage: psu target (PW arm)
      cov_p_nv_w = cp$naive["w"],  cov_p_nv_b = cp$naive["b"],  cov_p_nv_re = cp$naive["re"],
      cov_p_sel_w = cp$selective["w"], cov_p_sel_b = cp$selective["b"], cov_p_sel_re = cp$selective["re"],
      cov_p_blk_w = cp$blanket["w"], cov_p_blk_b = cp$blanket["b"], cov_p_blk_re = cp$blanket["re"],
      # CI width ratios for beta[2] (selective / naive)
      wr_bw_g = unname(r$ci_width$group$selective["beta[2]"] /
                       r$ci_width$group$naive["beta[2]"]),
      wr_bw_p = if (pw) unname(r$ci_width$psu$selective["beta[2]"] /
                               r$ci_width$psu$naive["beta[2]"]) else NA_real_,
      # flagged-RE coverage bookkeeping (psu target)
      n_flagged_re_p = n_fre, cov_fre_nv = fre_nv, cov_fre_sel = fre_sel,
      row.names = NULL
    )

    if (gate_pass) {
      der_sum <- if (is.null(der_sum)) dg else der_sum + dg
      der_n <- der_n + 1L
      re_pool_g <- c(re_pool_g, dg[idx_re])
      if (pw) re_pool_p <- c(re_pool_p, dp[idx_re])
    }
  }

  rep_rows[[ci]] <- do.call(rbind, rows)
  cell_extra[[ci]] <- list(
    der_param_mean = der_sum / der_n,       # per-parameter mean over pass reps
    n_pass = der_n,
    re_pooled_med_g = stats::median(re_pool_g),
    re_pooled_med_p = if (!is.null(re_pool_p)) stats::median(re_pool_p) else NA_real_
  )
  message(sprintf("  %-18s %d files, %d gate-pass", cell_id,
                  length(files), der_n))
}

reps <- do.call(rbind, rep_rows)
stopifnot(nrow(reps) == 24L * 200L)

# ---- Cell-level summary data frame -----------------------------------------
message("\n[2/6] Building the cell-level summary ...")

q95 <- function(x) unname(stats::quantile(x, 0.95, type = 7))

cell_summary <- do.call(rbind, lapply(seq_along(cells), function(ci) {
  cell_id <- cells[ci]
  a  <- reps[reps$cell_id == cell_id, ]
  p  <- a[a$gate_pass, ]                       # gate-pass reps only
  ex <- cell_extra[[cell_id]]
  pm <- ex$der_param_mean
  is_pw <- p$structure[1] == "psu_within_groups"

  # Separation: mean within-DER vs the largest per-parameter mean DER among
  # non-target parameters (fe_between + re), group target.
  max_nt_param_mean <- max(pm[-2])
  sep_factor <- pm[2] / max_nt_param_mean

  data.frame(
    cell_id = cell_id, J = p$J[1], nbar_j = p$nbar_j[1],
    informativeness = p$informativeness[1], structure = p$structure[1],
    n_reps = nrow(a), n_gate_fail = sum(!a$gate_pass), n_used = nrow(p),
    deff_mean_all = mean(a$deff), deff_mean_pass = mean(p$deff),
    cv_w_mean = mean(a$cv_w), N_mean = mean(a$N),
    sigma_theta_hat_mean = mean(p$sigma_theta_hat),
    timing_total_mean = mean(a$timing_total),
    # DER by class, group target
    der_w_g_mean = mean(p$der_w_g), der_w_g_med = stats::median(p$der_w_g),
    der_w_g_p05 = unname(stats::quantile(p$der_w_g, 0.05)),
    der_w_g_p95 = unname(stats::quantile(p$der_w_g, 0.95)),
    der_b_g_mean = mean(p$der_b_g), der_b_g_med = stats::median(p$der_b_g),
    der_b_g_p05 = unname(stats::quantile(p$der_b_g, 0.05)),
    der_b_g_p95 = unname(stats::quantile(p$der_b_g, 0.95)),
    der_int_g_mean = mean(p$der_int_g), der_bc_g_mean = mean(p$der_bc_g),
    der_re_med_g_mean = mean(p$der_re_med_g),
    der_re_med_g_p05 = unname(stats::quantile(p$der_re_med_g, 0.05)),
    der_re_med_g_p95 = unname(stats::quantile(p$der_re_med_g, 0.95)),
    der_re_pooled_med_g = ex$re_pooled_med_g,
    # DER by class, psu target (PW arm; NA in SG)
    der_w_p_mean = if (is_pw) mean(p$der_w_p) else NA_real_,
    der_b_p_mean = if (is_pw) mean(p$der_b_p) else NA_real_,
    der_re_med_p_mean = if (is_pw) mean(p$der_re_med_p) else NA_real_,
    der_re_pooled_med_p = ex$re_pooled_med_p,
    # flagged counts per rep (mean), both targets
    nflag_g_mean = mean(p$nflag_g),
    nflag_p_mean = if (is_pw) mean(p$nflag_p) else NA_real_,
    # flagged-RE coverage bookkeeping, PSU target (per-cell SUMS over
    # gate-pass reps), exported so the pooled naive/selective coverage of the
    # flagged random effects is recomputable from this CSV alone:
    # pooled naive = sum(fre_cov_nv_sum)/sum(fre_n_sum), etc.
    fre_n_sum       = if (is_pw) sum(p$n_flagged_re_p) else NA_real_,
    fre_cov_nv_sum  = if (is_pw) sum(p$cov_fre_nv) else NA_real_,
    fre_cov_sel_sum = if (is_pw) sum(p$cov_fre_sel) else NA_real_,
    # coverage, group target
    cov_g_nv_w = mean(p$cov_g_nv_w), cov_g_sel_w = mean(p$cov_g_sel_w), cov_g_blk_w = mean(p$cov_g_blk_w),
    cov_g_nv_b = mean(p$cov_g_nv_b), cov_g_sel_b = mean(p$cov_g_sel_b), cov_g_blk_b = mean(p$cov_g_blk_b),
    cov_g_nv_re = mean(p$cov_g_nv_re), cov_g_sel_re = mean(p$cov_g_sel_re), cov_g_blk_re = mean(p$cov_g_blk_re),
    # coverage, psu target (PW arm)
    cov_p_nv_w = if (is_pw) mean(p$cov_p_nv_w) else NA_real_,
    cov_p_sel_w = if (is_pw) mean(p$cov_p_sel_w) else NA_real_,
    cov_p_blk_w = if (is_pw) mean(p$cov_p_blk_w) else NA_real_,
    cov_p_nv_b = if (is_pw) mean(p$cov_p_nv_b) else NA_real_,
    cov_p_sel_b = if (is_pw) mean(p$cov_p_sel_b) else NA_real_,
    cov_p_blk_b = if (is_pw) mean(p$cov_p_blk_b) else NA_real_,
    cov_p_nv_re = if (is_pw) mean(p$cov_p_nv_re) else NA_real_,
    cov_p_sel_re = if (is_pw) mean(p$cov_p_sel_re) else NA_real_,
    cov_p_blk_re = if (is_pw) mean(p$cov_p_blk_re) else NA_real_,
    # width ratios (selective/naive) for beta[2]
    wr_bw_g_mean = mean(p$wr_bw_g),
    wr_bw_p_mean = if (is_pw) mean(p$wr_bw_p) else NA_real_,
    # max non-target DER distribution (group target)
    maxnt_min = min(p$max_nt_g), maxnt_mean = mean(p$max_nt_g),
    maxnt_p95 = q95(p$max_nt_g), maxnt_max = max(p$max_nt_g),
    maxnt_var = stats::var(p$max_nt_g),
    # false-positive rates at the c0 grid (group target)
    fp_c08 = mean(p$max_nt_g > 0.8), fp_c10 = mean(p$max_nt_g > 1.0),
    fp_c12 = mean(p$max_nt_g > 1.2), fp_c15 = mean(p$max_nt_g > 1.5),
    fp_c20 = mean(p$max_nt_g > 2.0),
    # flagged-set stability on c0 in [1.1, 1.6] (sets are nested in c0)
    stab_share_c11_c16 = mean(p$nflag_g_c11 == p$nflag_g_c16),
    stab_nt_share = mean(p$stab_nt_g),
    # separation factor (group target)
    sep_factor = sep_factor,
    row.names = NULL
  )
}))

saveRDS(cell_summary, file.path(TAB_DIR, "simulation_cell_summary.rds"))
utils::write.csv(cell_summary, file.path(TAB_DIR, "simulation_cell_summary.csv"),
                 row.names = FALSE)

# ---- Pooled quantities feeding the slots ------------------------------------
pass <- reps[reps$gate_pass, ]
pw_cells <- cell_summary[cell_summary$structure == "psu_within_groups", ]

n_gate_fail_total <- sum(!reps$gate_pass)
gate_pass_pct <- 100 * mean(reps$gate_pass)

# Realized DEFF by informativeness (all reps: a realized-design quantity)
deff_by_g <- tapply(reps$deff, reps$informativeness, mean)
deff_cell_rng <- lapply(split(cell_summary, cell_summary$informativeness),
                        function(d) range(d$deff_mean_all))

# Mean within-DER by informativeness (gate-pass, group target, pooled)
derw_by_g <- tapply(pass$der_w_g, pass$informativeness, mean)

# Between-class and RE summaries (group target).
# Per-coefficient cell means for the two between-group fixed effects (the
# manuscript sentence refers to "the between-group coefficient", i.e. beta_2
# in manuscript numbering = z_bc); the per-rep max-of-two class convention is
# kept as its own summary column.
der_b_percoef_max <- max(cell_summary$der_int_g_mean, cell_summary$der_bc_g_mean)
der_b_cellmean_max <- max(cell_summary$der_b_g_mean)   # per-rep max-of-two
re_med_rng <- range(cell_summary$der_re_med_g_mean)    # mean of per-rep medians
re_med_pooled_rng <- range(cell_summary$der_re_pooled_med_g)

# Dual-target flagged counts, PW arm (cell means)
flag_g_rng <- range(pw_cells$nflag_g_mean)
flag_p_rng <- range(pw_cells$nflag_p_mean)

# Flagged-RE coverage under the design-PSU target (pooled over PW pass reps)
pwp <- pass[pass$structure == "psu_within_groups", ]
n_fre_total  <- sum(pwp$n_flagged_re_p)
cov_fre_nv   <- sum(pwp$cov_fre_nv)  / n_fre_total
cov_fre_sel  <- sum(pwp$cov_fre_sel) / n_fre_total

# beta[2] coverage across informative cells (moderate + strong), group target
inf_cells <- cell_summary[cell_summary$informativeness != "none", ]
cov_w_nv_rng  <- range(inf_cells$cov_g_nv_w)
cov_w_sel_rng <- range(inf_cells$cov_g_sel_w)

# Blanket RE coverage range (group target, all 24 cells)
cov_blk_re_rng <- range(cell_summary$cov_g_blk_re)

# Separation factor and non-target tail (group target)
sep_min <- min(cell_summary$sep_factor)
sep_max <- max(cell_summary$sep_factor)
max_nt_pass <- max(pass$max_nt_g)
max_nt_all  <- max(reps$max_nt_g)

# Pooled false-positive rates at the c0 grid (group target, all cells)
fp_pooled <- vapply(C0_GRID, function(t) mean(pass$max_nt_g > t), numeric(1))
names(fp_pooled) <- paste0("c", C0_GRID)

# Flagged-set stability across c0 in [1.1, 1.6]: full set, and the
# protected (non-target) subset, whose flag decisions are the false positives.
stab_share    <- mean(pass$nflag_g_c11 == pass$nflag_g_c16)
stab_nt_share <- mean(pass$stab_nt_g)

# =============================================================================
# Section 2: NSECE -- refit readout and weight-convention comparison
# =============================================================================

message("\n[3/6] NSECE: reading the refit and supplement ...")

refit <- readRDS(file.path(PRECOMP, "application", "nsece_v2_refit.rds"))
supp  <- readRDS(file.path(PRECOMP, "application", "nsece_supp.rds"))
stopifnot(isTRUE(refit$gate$pass))

der_grp <- as.numeric(refit$fit_grp$der)   # beta[1:3], theta[1:51]
der_psu <- as.numeric(refit$fit_psu$der)
J_ns <- refit$stan_data_meta$J
stopifnot(length(der_grp) == 3 + J_ns)
idx_re_ns <- 3 + seq_len(J_ns)

flag_grp <- sum(der_grp > C0_DEFAULT)
flag_psu <- sum(der_psu > C0_DEFAULT)
flag_re_psu <- sum(der_psu[idx_re_ns] > C0_DEFAULT)
flag_re_grp <- sum(der_grp[idx_re_ns] > C0_DEFAULT)
re_med_psu <- stats::median(der_psu[idx_re_ns])
re_med_grp <- stats::median(der_grp[idx_re_ns])
shielded_max <- max(der_grp[c(1, 3)], der_psu[c(1, 3)])

# beta_1 (poverty) 90% CI width factor under the design-PSU correction
ciq <- function(x) unname(stats::quantile(x, c(0.05, 0.95)))
b2_nv  <- ciq(refit$draws[, 2])
b2_psu <- ciq(refit$fit_psu$corrected_draws[, 2])
width_factor_psu <- diff(b2_psu) / diff(b2_nv)

# Blanket damage under the model-group target: marginal rescale by
# sqrt(DER_k), so the width ratio is exactly sqrt(DER_k).
wr_blanket_grp <- sqrt(der_grp)
narrow_idx   <- which(wr_blanket_grp < 1)
n_narrowed   <- length(narrow_idx)
narrow_mean  <- mean(wr_blanket_grp[narrow_idx])
narrow_worst <- min(wr_blanket_grp[narrow_idx])

# c0 sweep, NSECE (both targets; all 54 params and the RE-only count)
sweep_grp    <- vapply(C0_GRID, function(t) sum(der_grp > t), integer(1))
sweep_psu    <- vapply(C0_GRID, function(t) sum(der_psu > t), integer(1))
sweep_psu_re <- vapply(C0_GRID, function(t) sum(der_psu[idx_re_ns] > t), integer(1))

# ---- Weight-convention sensitivity and structural factors -------------------
# Table 4 contrasts the declared global unit-mean convention with the
# alternative within-state scaling. The within-state quantities, the effective
# per-state sample sizes, and the structural protection factors R_k are formed
# from the restricted microdata in
# code/02_application/app_03_der_supplement.R and read here from the shipped
# (non-identifiable) supplement, so this reporting stage needs no restricted
# data and no model refitting.
conv <- list(
  grp_der_b1     = supp$convention$within_der_b1_grp,
  psu_der_b1     = supp$convention$within_der_b1_psu,
  grp_flagged    = supp$convention$within_flag_grp,
  psu_flagged    = supp$convention$within_flag_psu,
  psu_flagged_re = supp$convention$within_flag_re_psu,
  grp_flagged_re = supp$convention$within_flag_re_grp
)
eff_glob   <- supp$eff_global   # named numeric c(min = ., max = .)
eff_within <- supp$eff_within
rk_within  <- supp$R_k_within
rk_between <- supp$R_k_between

# =============================================================================
# Section 3: Figures (Part B)
# =============================================================================

message("\n[4/6] Figures ...")

theme_s4 <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom",
        plot.title = element_blank())

lab_structure <- c(sampled_groups   = "Sampled~groups~(c == j)",
                   psu_within_groups = "PSUs~within~groups~(c != j)")
lab_inf <- c(none = "None", moderate = "Moderate", strong = "Strong")

# ---- fig2: DER separation profile (group target) ----------------------------
fig2_dat <- do.call(rbind, lapply(c("fe_within", "fe_between", "re"), function(cl) {
  v <- switch(cl,
    fe_within  = cell_summary[, c("der_w_g_mean", "der_w_g_p05", "der_w_g_p95")],
    fe_between = cell_summary[, c("der_b_g_mean", "der_b_g_p05", "der_b_g_p95")],
    re         = cell_summary[, c("der_re_med_g_mean", "der_re_med_g_p05",
                                  "der_re_med_g_p95")])
  names(v) <- c("y", "lo", "hi")
  cbind(cell_summary[, c("cell_id", "J", "nbar_j", "informativeness",
                         "structure")],
        v, class = cl)
}))
fig2_dat$class <- factor(fig2_dat$class, levels = c("fe_within", "fe_between", "re"),
  labels = c("Within FE (mean)", "Between FE (mean of max)", "RE (mean of medians)"))
fig2_dat$size_lab <- factor(sprintf("J%d, %d", fig2_dat$J, fig2_dat$nbar_j),
  levels = c("J20, 10", "J20, 50", "J50, 10", "J50, 50"))
fig2_dat$structure_f <- factor(lab_structure[fig2_dat$structure],
                               levels = lab_structure)
fig2_dat$inf_f <- factor(lab_inf[fig2_dat$informativeness], levels = lab_inf)

fig2 <- ggplot(fig2_dat, aes(x = size_lab, y = y, colour = class, shape = class)) +
  geom_hline(yintercept = C0_DEFAULT, linetype = "dashed",
             colour = "grey35", linewidth = 0.45) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0,
                position = position_dodge(width = 0.62), linewidth = 0.45,
                alpha = 0.85) +
  geom_point(position = position_dodge(width = 0.62), size = 2.1) +
  facet_grid(structure_f ~ inf_f, labeller = labeller(
    structure_f = label_parsed, inf_f = label_value)) +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3),
                labels = c("0.01", "0.03", "0.1", "0.3", "1", "3")) +
  scale_colour_manual(values = unname(OI[c("vermillion", "blue", "green")])) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = expression(paste("Design cell:  ", J, ",  ", bar(n)[j])),
       y = expression(paste(DER, "  (model-group target, log scale)")),
       colour = NULL, shape = NULL) +
  theme_s4 +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# c0 label once, in the top-left panel (just below the line, clear of points)
fig2 <- fig2 + geom_text(
  data = data.frame(x = 0.42, y = 0.78,
                    structure_f = factor(lab_structure[1], levels = lab_structure),
                    inf_f = factor("None", levels = lab_inf)),
  aes(x = x, y = y), label = "c[0] == 1.2", parse = TRUE,
  colour = "grey25", size = 3.2, hjust = 0, inherit.aes = FALSE)

ggsave(file.path(FIG_DIR, "fig2_der_separation.pdf"), fig2,
       width = 9.0, height = 5.6, device = "pdf", useDingbats = FALSE)
ggsave(file.path(FIG_DIR, "fig2_der_separation.png"), fig2,
       width = 9.0, height = 5.6, dpi = 150)

# ---- fig3: coverage triptych (group target) ---------------------------------
mk_cov_long <- function() {
  strategies <- c(nv = "Naive", sel = "Selective", blk = "Blanket")
  classes <- c(w = "beta1_within", b = "fe_between", re = "re_avg")
  out <- list()
  for (s in names(strategies)) for (cl in names(classes)) {
    col <- paste0("cov_g_", s, "_", cl)
    agg <- aggregate(pass[[col]],
                     by = list(informativeness = pass$informativeness,
                               structure = pass$structure), FUN = mean)
    n_agg <- aggregate(rep(1, nrow(pass)),
                       by = list(informativeness = pass$informativeness,
                                 structure = pass$structure), FUN = sum)
    agg$n <- n_agg$x
    agg$strategy <- strategies[[s]]
    agg$class <- classes[[cl]]
    out[[paste(s, cl)]] <- agg
  }
  do.call(rbind, out)
}
cov_long <- mk_cov_long()
cov_long$se <- sqrt(cov_long$x * (1 - cov_long$x) / cov_long$n)
cov_long$class <- factor(cov_long$class,
  levels = c("beta1_within", "fe_between", "re_avg"),
  labels = c("beta[1]~(within)", "Between~FE", "RE~average"))
cov_long$strategy <- factor(cov_long$strategy,
                            levels = c("Naive", "Selective", "Blanket"))
cov_long$inf_f <- factor(lab_inf[cov_long$informativeness], levels = lab_inf)
cov_long$structure_f <- factor(lab_structure[cov_long$structure],
                               levels = lab_structure)

fig3 <- ggplot(cov_long, aes(x = inf_f, y = x, colour = strategy,
                             shape = strategy, group = strategy)) +
  geom_hline(yintercept = CI_LEVEL, linetype = "dashed", colour = "grey35",
             linewidth = 0.45) +
  geom_line(position = position_dodge(width = 0.35), linewidth = 0.4,
            alpha = 0.6) +
  geom_errorbar(aes(ymin = x - 1.96 * se, ymax = x + 1.96 * se), width = 0.12,
                position = position_dodge(width = 0.35), linewidth = 0.45) +
  geom_point(position = position_dodge(width = 0.35), size = 2.2) +
  facet_grid(structure_f ~ class, labeller = label_parsed) +
  scale_colour_manual(values = unname(OI[c("vermillion", "green", "skyblue")])) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_y_continuous(limits = c(0.12, 1.0),
                     breaks = seq(0.2, 1.0, by = 0.2)) +
  labs(x = "Informativeness of the sampling mechanism",
       y = "90% interval coverage (model-group target)",
       colour = "Strategy", shape = "Strategy") +
  theme_s4

ggsave(file.path(FIG_DIR, "fig3_coverage.pdf"), fig3,
       width = 9.0, height = 5.4, device = "pdf", useDingbats = FALSE)
ggsave(file.path(FIG_DIR, "fig3_coverage.png"), fig3,
       width = 9.0, height = 5.4, dpi = 150)

# ---- fig4: NSECE dual-target profile ----------------------------------------
ord_re <- order(der_grp[idx_re_ns])
fig4_dat <- rbind(
  data.frame(pos = 1:3, param = c("beta0", "beta1", "beta2"),
             der = der_grp[1:3], target = "Model-group (51 states)",
             type = "FE"),
  data.frame(pos = 1:3, param = c("beta0", "beta1", "beta2"),
             der = der_psu[1:3], target = "Design-PSU (415 PSUs, 30 strata)",
             type = "FE"),
  data.frame(pos = 3 + seq_len(J_ns), param = paste0("theta", ord_re),
             der = der_grp[idx_re_ns][ord_re],
             target = "Model-group (51 states)", type = "RE"),
  data.frame(pos = 3 + seq_len(J_ns), param = paste0("theta", ord_re),
             der = der_psu[idx_re_ns][ord_re],
             target = "Design-PSU (415 PSUs, 30 strata)", type = "RE")
)
fig4_dat$target <- factor(fig4_dat$target,
  levels = c("Design-PSU (415 PSUs, 30 strata)", "Model-group (51 states)"))

fe_lab <- data.frame(pos = 1:3, der = pmax(der_grp[1:3], der_psu[1:3]),
                     lab = c("beta[0]", "beta[1]", "beta[2]"))

fig4 <- ggplot(fig4_dat, aes(x = pos, y = der, colour = target, shape = target)) +
  geom_hline(yintercept = C0_DEFAULT, linetype = "dashed", colour = "grey35",
             linewidth = 0.45) +
  geom_vline(xintercept = 3.5, colour = "grey80", linewidth = 0.4) +
  geom_point(size = 1.9, alpha = 0.9) +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1.2, 3),
                labels = c("0.01", "0.03", "0.1", "0.3", "1.2", "3"),
                limits = c(min(c(der_grp, der_psu)) * 0.7, 40)) +
  scale_x_continuous(breaks = c(2, 29),
                     labels = c("Fixed effects", "State random effects (sorted by model-group DER)"),
                     expand = expansion(add = c(1.2, 1.2))) +
  scale_colour_manual(values = unname(OI[c("vermillion", "blue")])) +
  scale_shape_manual(values = c(17, 16)) +
  annotate("text", x = 54, y = 1.2, label = "c[0] == 1.2", parse = TRUE,
           colour = "grey25", size = 3.2, hjust = 1, vjust = -0.5) +
  annotate("text", x = 30, y = 22,
           label = sprintf("%d of %d state effects flagged under the design-PSU target",
                           flag_re_psu, J_ns),
           colour = OI[["vermillion"]], size = 3.4) +
  annotate("text", x = 33, y = 0.0055,
           label = sprintf("%d of %d flagged under the model-group target (median DER = %s)",
                           flag_re_grp, J_ns, fmt(re_med_grp, 2)),
           colour = OI[["blue"]], size = 3.4) +
  geom_text(data = fe_lab, aes(x = pos, y = der * 1.7, label = lab),
            parse = TRUE, inherit.aes = FALSE, size = 3.4, colour = "grey15") +
  labs(x = NULL, y = "DER (log scale)", colour = "Variance target",
       shape = "Variance target") +
  theme_s4 +
  theme(axis.ticks.x = element_blank())

ggsave(file.path(FIG_DIR, "fig4_nsece_dualtarget.pdf"), fig4,
       width = 9.0, height = 4.8, device = "pdf", useDingbats = FALSE)
ggsave(file.path(FIG_DIR, "fig4_nsece_dualtarget.png"), fig4,
       width = 9.0, height = 4.8, dpi = 150)

# Figure 1 (the closed-form theory curves) is produced by results_01_theory.R.

# =============================================================================
# Section 4: Tables (Part C)
# =============================================================================

message("\n[5/6] Tables ...")

tex_header <- function(script = "s4_aggregate.R") {
  sprintf("%%%% Generated by %s on %s -- do not hand-edit.\n",
          script, format(Sys.time(), "%Y-%m-%d"))
}

# ---- tab2: simulation coverage summary --------------------------------------
tab2_block <- function(struct) {
  rows <- character(0)
  for (g in c("none", "moderate", "strong")) {
    sub <- pass[pass$structure == struct & pass$informativeness == g, ]
    n_used <- nrow(sub)
    n_fail <- 4 * 200 - n_used
    for (s in c("nv", "sel", "blk")) {
      sname <- c(nv = "Naive", sel = "Selective", blk = "Blanket")[[s]]
      cw  <- mean(sub[[paste0("cov_g_", s, "_w")]])
      cb  <- mean(sub[[paste0("cov_g_", s, "_b")]])
      cre <- mean(sub[[paste0("cov_g_", s, "_re")]])
      lead <- if (s == "nv") {
        sprintf("%s & %d", c(none = "None", moderate = "Moderate",
                             strong = "Strong")[[g]], n_used)
      } else " & "
      rows <- c(rows, sprintf("%s & %s & %s & %s & %s \\\\",
                              lead, sname, fmt(cw, 3), fmt(cb, 3), fmt(cre, 3)))
    }
    if (g != "strong") rows <- c(rows, "\\addlinespace")
  }
  rows
}

tab2 <- c(
  tex_header(),
  "\\begin{table}[t]",
  "\\centering",
  paste0("\\caption{Simulation study: 90\\% interval coverage under the ",
         "model-group target by informativeness, structural arm, and interval ",
         "strategy, pooled over $J \\in \\{20, 50\\}$ and $\\bar n_j \\in ",
         "\\{10, 50\\}$ (four cells $\\times$ 200 replications per row block). ",
         "``Reps'' counts gate-passing replications out of 800; gate failures ",
         sprintf("(%d of 4{,}800 in total) are excluded from all summaries. ",
                 n_gate_fail_total),
         "The selective strategy corrects the flagged block only (threshold ",
         "$c_0 = 1.2$); the blanket strategy rescales every parameter ",
         "marginally to the target.}"),
  "\\label{tab:sim-coverage}",
  "\\begin{tabular}{@{}llcccc@{}}",
  "\\toprule",
  " & & & \\multicolumn{3}{c}{Coverage (nominal 0.90)} \\\\",
  "\\cmidrule(lr){4-6}",
  paste0("Informativeness & Reps & Strategy & $\\beta_1$ (within) & ",
         "FE between & RE average \\\\"),
  "\\midrule",
  "\\multicolumn{6}{@{}l}{\\itshape Sampled groups ($c = j$)} \\\\",
  tab2_block("sampled_groups"),
  "\\midrule",
  "\\multicolumn{6}{@{}l}{\\itshape PSUs within groups ($c \\neq j$)} \\\\",
  tab2_block("psu_within_groups"),
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{4pt}",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  paste0("\\textit{Notes.} ``FE between'' averages the intercept and the ",
         "between-group coefficient $\\beta_2$; ``RE average'' averages the ",
         "$J$ random intercepts. Monte Carlo standard errors are at most ",
         "about 2 percentage points per cell (200 replications) and about 1 ",
         "percentage point for the pooled entries shown here ($\\approx 800$ ",
         "replications)."),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(unlist(tab2), file.path(TAB_DIR, "tab2_sim_coverage.tex"))

# ---- tab3: NSECE headline dual-target table ---------------------------------
flagmark <- function(flag) if (flag) "$\\checkmark$" else "---"
tab3 <- c(
  tex_header(),
  "\\begin{table}[t]",
  "\\centering",
  paste0("\\caption{NSECE 2019: $\\DER$ diagnostics for all 54 parameters ",
         "under both declared variance targets (weights under the declared ",
         "global unit-mean convention on the raw design weights, raw-weight ",
         "Kish $\\DEFF = 3.76$; strict convergence gate passed: ",
         "$\\widehat{R} \\le 1.004$ on all parameters, zero divergent ",
         "transitions). Flags mark $\\DER > c_0 = 1.2$.}"),
  "\\label{tab:nsece_der}",
  "\\begin{tabular}{@{}llcccc@{}}",
  "\\toprule",
  paste0(" & & \\multicolumn{2}{c}{Model-group target} & ",
         "\\multicolumn{2}{c}{Design-PSU target} \\\\"),
  "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}",
  "Parameter & Type & $\\DER$ & Flag & $\\DER$ & Flag \\\\",
  "\\midrule",
  sprintf("$\\beta_0$ (intercept) & Between & %s & %s & %s & %s \\\\",
          fmt(der_grp[1], 3), flagmark(der_grp[1] > 1.2),
          fmt(der_psu[1], 3), flagmark(der_psu[1] > 1.2)),
  sprintf("$\\beta_1$ (poverty, CWC) & Within & %s & %s & %s & %s \\\\",
          fmt(der_grp[2], 3), flagmark(der_grp[2] > 1.2),
          fmt(der_psu[2], 3), flagmark(der_psu[2] > 1.2)),
  sprintf("$\\beta_2$ (tiered reimb.) & Between & %s & %s & %s & %s \\\\",
          fmt(der_grp[3], 3), flagmark(der_grp[3] > 1.2),
          fmt(der_psu[3], 3), flagmark(der_psu[3] > 1.2)),
  "\\addlinespace",
  sprintf(paste0("$\\theta_j$ (51 state effects) & RE & %s & %d of 51 & ",
                 "%s & %d of 51 \\\\"),
          paste0("median ", fmt(re_med_grp, 3)), flag_re_grp,
          paste0("median ", fmt(re_med_psu, 2)), flag_re_psu),
  "\\midrule",
  sprintf(paste0("Flagged total ($\\DER > 1.2$) & & ",
                 "\\multicolumn{2}{c}{%d of 54} & ",
                 "\\multicolumn{2}{c}{%d of 54} \\\\"),
          flag_grp, flag_psu),
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{4pt}",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  paste0("\\textit{Notes.} Model-group target: clusters $=$ 51 states, single ",
         "stratum. Design-PSU target: clusters $=$ 415 design PSUs within 30 ",
         "strata. Both use the centered, DF-corrected meat and the declared ",
         "global unit-mean weight convention; the within-state alternative is ",
         "examined in Table~\\ref{tab:nsece_convention}. ",
         sprintf("For $\\beta_1$, the selective correction under the design-PSU target widens the 90\\%% credible interval from $[%s, %s]$ to $[%s, %s]$ (factor %s).",
                 fmt(b2_nv[1], 3), fmt(b2_nv[2], 3),
                 fmt(b2_psu[1], 3), fmt(b2_psu[2], 3),
                 fmt(width_factor_psu, 2))),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tab3, file.path(TAB_DIR, "tab3_nsece_der.tex"))

# ---- tab4: weight-convention 2x2 --------------------------------------------
tab4 <- c(
  tex_header(),
  "\\begin{table}[t]",
  "\\centering",
  paste0("\\caption{NSECE 2019: sensitivity of the diagnostic to the declared ",
         "weight convention. The two rows are two \\emph{different} ",
         "pseudo-posteriors---the weights enter the pseudo-likelihood itself, ",
         "so changing the convention changes the fit, not merely the ",
         "target---each evaluated under both aggregation units with identical ",
         "target settings (centered, DF-corrected meat; unit-mean ",
         "normalization).}"),
  "\\label{tab:nsece_convention}",
  "\\begin{tabular}{@{}lcccc@{}}",
  "\\toprule",
  paste0(" & \\multicolumn{2}{c}{$\\DER(\\beta_1)$ (poverty)} & ",
         "\\multicolumn{2}{c}{Flagged (of 54)} \\\\"),
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Weight convention & Model-group & Design-PSU & Model-group & Design-PSU \\\\",
  "\\midrule",
  sprintf("Global unit-mean (declared) & %s & %s & %d & %d \\\\",
          fmt(der_grp[2], 2), fmt(der_psu[2], 2), flag_grp, flag_psu),
  sprintf("Within-state scaling & %s & %s & %d & %d \\\\",
          fmt(conv$grp_der_b1, 2), fmt(conv$psu_der_b1, 2),
          conv$grp_flagged, conv$psu_flagged),
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{4pt}",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  paste0("\\textit{Notes.} Within-state scaling rescales weights to sum to ",
         "each state's sample size (the multilevel scaling convention of ",
         "\\citealp{PfeffermannEtAl1998b}). Under it the state pseudo-likelihood ",
         "sample sizes are the realized $n_j$ ",
         sprintf("($[%d, %s]$), whereas the global convention allocates them by weighted share ($[%d, %d]$): ",
                 round(eff_within["min"]), lint(eff_within["max"]),
                 round(eff_glob["min"]), round(eff_glob["max"])),
         "the convention reallocates effective information across states, ",
         "which reverses the ordering of the two targets for $\\beta_1$. ",
         sprintf("The poverty coefficient is flagged under all four combinations; the random-effect contrast (flagged only under the design-PSU target: %d and %d of 51 states) persists.",
                 flag_re_psu, conv$psu_flagged_re)),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tab4, file.path(TAB_DIR, "tab4_nsece_2x2.tex"))

# ---- tab5: c0 sweep ----------------------------------------------------------
tab5_rows <- vapply(seq_along(C0_GRID), function(i) {
  c0 <- C0_GRID[i]
  bold <- function(x) if (c0 == C0_DEFAULT) paste0("\\textbf{", x, "}") else x
  sprintf("%s & %s & %s & %s \\\\",
          bold(fmt(c0, 1)), bold(sprintf("%d", sweep_grp[i])),
          bold(sprintf("%d", sweep_psu[i])),
          bold(pctf(fp_pooled[i], 1)))
}, character(1))

tab5 <- c(
  tex_header(),
  "\\begin{table}[t]",
  "\\centering",
  paste0("\\caption{Threshold sensitivity. NSECE: number of parameters (of ",
         "54) flagged at each threshold $c_0$ under both declared targets. ",
         "Simulation: per-replication false-positive rate at the same ",
         "thresholds under the model-group target (share of gate-passing ",
         "replications in which any non-target parameter---between-group ",
         "fixed effect or random effect---exceeds $c_0$; pooled over all 24 ",
         "cells). The recommended default $c_0 = 1.2$ is shown in bold.}"),
  "\\label{tab:c0_sensitivity}",
  "\\begin{tabular}{@{}cccc@{}}",
  "\\toprule",
  paste0(" & \\multicolumn{2}{c}{NSECE flagged (of 54)} & Simulation \\\\"),
  "\\cmidrule(lr){2-3}",
  paste0("$c_0$ & Model-group & Design-PSU & False-positive rate (\\%) \\\\"),
  "\\midrule",
  tab5_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{4pt}",
  "\\begin{minipage}{\\textwidth}",
  "\\footnotesize",
  paste0("\\textit{Notes.} NSECE counts are recomputed from the full $\\DER$ ",
         "vectors of the v2 refit. In the NSECE the poverty coefficient ",
         "$\\beta_1$ is flagged at every threshold shown, under both targets. ",
         sprintf("Simulation false positives are based on %s gate-passing replications; ", lint(nrow(pass))),
         "a false positive costs a wider interval for a protected ",
         "(unflagged) parameter, not invalid inference."),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tab5, file.path(TAB_DIR, "tab5_c0_sweep.tex"))

# =============================================================================
# Section 5: Integrity summary
# =============================================================================

message("\n[6/6] Integrity summary")

ok <- function(label, pass, detail = "") {
  message(sprintf("  [%s] %-50s %s", if (isTRUE(pass)) "OK" else "!!",
                  label, detail))
  invisible(pass)
}

ok("24 design cells x 200 replications loaded", nrow(reps) == 24L * 200L,
   sprintf("(%d reps)", nrow(reps)))
ok("gate-pass rate reported", TRUE,
   sprintf("%.1f%% (%d of %d excluded)", gate_pass_pct, n_gate_fail_total,
           nrow(reps)))
ok("NSECE flagged, model-group target = 1", flag_grp == 1L,
   sprintf("(= %d)", flag_grp))
ok("NSECE flagged, design-PSU target = 28", flag_psu == 28L,
   sprintf("(= %d)", flag_psu))
ok("NSECE refit re-derivation gap ~ 0", max(supp$reproduce_gap) < 1e-6,
   sprintf("(max |diff| = %.1e)", max(supp$reproduce_gap)))

message("\nReproduced headline quantities:")
message(sprintf("  Realized Kish DEFF: none %.2f | moderate %.2f | strong %.2f",
                deff_by_g[["none"]], deff_by_g[["moderate"]], deff_by_g[["strong"]]))
message(sprintf("  Within-coefficient DER (group target): none %.1f -> strong %.1f",
                derw_by_g[["none"]], derw_by_g[["strong"]]))
message(sprintf("  Separation factor (cell range): %.1f to %.1f", sep_min, sep_max))
message(sprintf("  False-positive rate at c0 = 1.2 (pooled): %.1f%%",
                100 * fp_pooled[["c1.2"]]))
message(sprintf("  beta_1 selective CI width factor (design-PSU): %.2f",
                width_factor_psu))

message(sprintf("\nFigures -> %s", FIG_DIR))
message(sprintf("Tables  -> %s", TAB_DIR))
message(sprintf("Done in %.1f s.", proc.time()["elapsed"] - t_run0))
