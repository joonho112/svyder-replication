# =============================================================================
# app_03_der_supplement.R -- Supplementary NSECE quantities (Track A)
# Author : JoonHo Lee (jlee296@ua.edu)
# =============================================================================
#
# The headline NSECE table (Table 3) and the dual-target figure (Figure 4) are
# built directly from the shipped model output "nsece_v2_refit.rds". Two
# further results in the paper need quantities that can only be formed from the
# restricted-use microdata, so we compute them here, once, and store the
# results (a handful of scalars and short vectors -- nothing respondent-level)
# for the reporting stage to read:
#
#   * Table 4 (weight-convention sensitivity): the Design Effect Ratio of the
#     poverty coefficient, and the flagged-parameter counts, under the
#     alternative "within-state" weight-scaling convention. Because the weights
#     enter the pseudo-likelihood itself, the alternative convention is a
#     different pseudo-posterior; here we evaluate the two declared targets on
#     the stored draws of that fit.
#   * The structural protection factors R_k of the fixed effects (Theorem 1),
#     and the effective per-state sample sizes under each convention.
#
# This is a Track A step: it requires the restricted NSECE derivatives. Readers
# working from the shipped outputs (Track B) can skip it entirely -- the result
# file "data/precomputed/application/nsece_supp.rds" is included in the package.
#
# Inputs  : nsece_v2_refit.rds (shipped model output)
#           <restricted>/phase3_data.rds       y, X, group, within-state w
#           <restricted>/phase3_analysis.rds   raw weight, psu / stratum indices
#           <restricted>/phase3_corrected_draws.rds  draws of the within-state fit
#           <restricted>/phase3_sandwich.rds   plug-in sigma_theta
# Output  : data/precomputed/application/nsece_supp.rds
#
# The location of the restricted derivatives is taken from the environment
# variable NSECE_RESTRICTED_DIR (see data/DATA_ACCESS.md); it defaults to
# "data/nsece", which is git-ignored.
# =============================================================================

suppressMessages(library(svyder))

here_root <- tryCatch(here::here(), error = function(e) getwd())
restricted_dir <- Sys.getenv("NSECE_RESTRICTED_DIR",
                             file.path(here_root, "data", "nsece"))
refit_path <- file.path(here_root, "data", "precomputed", "application",
                        "nsece_v2_refit.rds")
out_path   <- file.path(here_root, "data", "precomputed", "application",
                        "nsece_supp.rds")

need <- file.path(restricted_dir,
                  c("phase3_data.rds", "phase3_analysis.rds",
                    "phase3_corrected_draws.rds", "phase3_sandwich.rds"))
if (!all(file.exists(need))) {
  stop("Restricted NSECE derivatives not found under '", restricted_dir,
       "'. This is a Track A step; see data/DATA_ACCESS.md. Track B users can ",
       "use the shipped data/precomputed/application/nsece_supp.rds instead.",
       call. = FALSE)
}

C0 <- 1.2
pt_fe <- c("fe_between", "fe_within", "fe_between")   # intercept, poverty, tiered

refit <- readRDS(refit_path)
d  <- readRDS(file.path(restricted_dir, "phase3_data.rds"))            # y, X, group, w_within
a  <- readRDS(file.path(restricted_dir, "phase3_analysis.rds"))        # weight, psu_idx, stratum_idx
cd <- readRDS(file.path(restricted_dir, "phase3_corrected_draws.rds")) # draws_naive
sw <- readRDS(file.path(restricted_dir, "phase3_sandwich.rds"))        # sigma_theta_hat
stopifnot(all(d$group == a$state_idx))

J  <- refit$stan_data_meta$J
idx_re <- 3 + seq_len(J)

# ---- Table 4: the two declared targets under the within-state convention ----
# "d$w" carries phase3_data$w: the within-state (per-state-scaled) weights of
# the alternative fit. The headline v2 convention is rebuilt below from
# phase3_analysis$weight as w_global.
conv_psu <- der_compute(
  cd$draws_naive, y = d$y, X = d$X, group = d$group,
  weights = d$w, cluster = a$psu_idx, strata = a$stratum_idx,
  family = "binomial", sigma_theta = sw$sigma_theta_hat,
  normalize = "unit_mean", beta_prior_sd = Inf, param_types = pt_fe)
conv_grp <- der_compute(
  cd$draws_naive, y = d$y, X = d$X, group = d$group,
  weights = d$w, cluster = d$group,
  family = "binomial", sigma_theta = sw$sigma_theta_hat,
  normalize = "unit_mean", beta_prior_sd = Inf, param_types = pt_fe)

convention <- list(
  within_der_b1_grp  = unname(conv_grp$der[2]),
  within_der_b1_psu  = unname(conv_psu$der[2]),
  within_flag_grp    = sum(conv_grp$der > C0),
  within_flag_psu    = sum(conv_psu$der > C0),
  within_flag_re_psu = sum(conv_psu$der[idx_re] > C0),
  within_flag_re_grp = sum(conv_grp$der[idx_re] > C0))

# ---- Effective per-state sample sizes under each convention -----------------
w_raw    <- a$weight
w_global <- w_raw / mean(w_raw)     # declared (global unit-mean)
w_within <- d$w                     # alternative (within-state scaling)
eff_global <- tapply(w_global, d$group, sum)
eff_within <- tapply(w_within, d$group, sum)

# ---- Structural protection factors R_k (Theorem 1) --------------------------
# Reproduce the declared-convention DERs from the stored draws (a deterministic
# re-derivation that also serves as an integrity check against the shipped
# refit), then read off the bread-side decomposition.
v2_psu <- der_compute(
  refit$draws, y = d$y, X = d$X, group = d$group, weights = w_global,
  cluster = a$psu_idx, strata = a$stratum_idx, family = "binomial",
  sigma_theta = refit$sigma_theta_hat, normalize = "unit_mean",
  param_types = pt_fe)
v2_grp <- der_compute(
  refit$draws, y = d$y, X = d$X, group = d$group, weights = w_global,
  cluster = d$group, family = "binomial",
  sigma_theta = refit$sigma_theta_hat, normalize = "unit_mean",
  param_types = pt_fe)

reproduce_gap <- c(group = max(abs(v2_grp$der - refit$fit_grp$der)),
                   psu   = max(abs(v2_psu$der - refit$fit_psu$der)))
decomp <- der_decompose(v2_psu)      # R_k is target-independent (bread-side)

supp <- list(
  convention    = convention,
  eff_global    = c(min = unname(min(eff_global)), max = unname(max(eff_global))),
  eff_within    = c(min = unname(min(eff_within)), max = unname(max(eff_within))),
  R_k_within    = unname(decomp$R_k[2]),
  R_k_between   = unname(decomp$R_k[c(1, 3)]),
  reproduce_gap = reproduce_gap)

saveRDS(supp, out_path)

cat(sprintf("Wrote %s\n", out_path))
cat(sprintf("  within-state DER(beta_1): model-group %.2f, design-PSU %.2f\n",
            convention$within_der_b1_grp, convention$within_der_b1_psu))
cat(sprintf("  within-state flagged: model-group %d, design-PSU %d\n",
            convention$within_flag_grp, convention$within_flag_psu))
cat(sprintf("  R_k: poverty (within) %.2f; intercept/tiered (between) %.2f, %.2f\n",
            supp$R_k_within, supp$R_k_between[1], supp$R_k_between[2]))
cat(sprintf("  refit re-derivation gap (should be ~0): group %.2e, psu %.2e\n",
            reproduce_gap["group"], reproduce_gap["psu"]))
