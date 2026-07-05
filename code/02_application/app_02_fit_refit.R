# =============================================================================
# app_02_fit_refit.R -- NSECE 2019 model fit and dual-target DER (Track A)
# Author : JoonHo Lee (jlee296@ua.edu)
# =============================================================================
# Refits the NSECE hierarchical logistic regression under the declared weight
# convention (global unit-mean on the raw design weights). The weights enter the
# pseudo-likelihood, so this is a genuine refit rather than a reweighting of
# stored draws. The script then applies the strict convergence gate and runs the
# dual-target DER analysis (design-PSU and model-group targets) with svyder.
#
# This is a Track A step: it needs the restricted NSECE derivatives and compiles
# the Stan model. The result it produces, nsece_v2_refit.rds, is shipped with
# the package, so readers reproducing the figures and tables (Track B) do not
# have to run it. The restricted derivatives are located via the environment
# variable NSECE_RESTRICTED_DIR (see data/DATA_ACCESS.md); it defaults to the
# git-ignored data/nsece/.
#
# Usage : Rscript code/02_application/app_02_fit_refit.R
# =============================================================================

suppressMessages({
  library(cmdstanr)
  library(posterior)
  library(svyder)
})

t_start   <- Sys.time()
here_root <- tryCatch(here::here(), error = function(e) getwd())
restricted_dir <- Sys.getenv("NSECE_RESTRICTED_DIR",
                             file.path(here_root, "data", "nsece"))
out_dir <- file.path(here_root, "data", "precomputed", "application")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Inputs (restricted NSECE derivatives) ----------------------------------
d <- readRDS(file.path(restricted_dir, "phase3_data.rds"))     # y, X, group (N = 6,785)
a <- readRDS(file.path(restricted_dir, "phase3_analysis.rds")) # raw weight, psu / stratum
stopifnot(length(d$y) == nrow(a), all(d$group == a$state_idx))

# ---- Declared convention: global unit-mean on the raw design weights --------
# The headline v2 fit uses phase3_analysis$weight, not phase3_data$w. The latter
# stores the within-state sensitivity convention used by app_03_der_supplement.R.
w_raw <- a$weight
w_v2  <- w_raw / mean(w_raw)
kish_raw <- length(w_raw) * sum(w_raw^2) / sum(w_raw)^2
cat(sprintf("Raw-weight global Kish DEFF = %.4f (paper: 3.76)\n", kish_raw))

stan_data <- list(N = d$N, J = d$J, p = d$p,
                  y = d$y, X = d$X, group = d$group, w = w_v2)

# ---- Fit --------------------------------------------------------------------
mod <- cmdstan_model(file.path(here_root, "stan", "hlr_weighted.stan"))
fit <- mod$sample(
  data = stan_data, seed = 20260705,
  chains = 4, parallel_chains = 4,
  iter_warmup = 2000, iter_sampling = 2000,
  adapt_delta = 0.95, refresh = 500
)

# ---- Strict convergence gate (all parameters) --------------------------------
sm <- fit$summary(c("beta", "sigma_theta", "theta", "eta_raw"))
diag_sum <- fit$diagnostic_summary()
gate <- list(
  rhat_max     = max(sm$rhat, na.rm = TRUE),
  ess_bulk_min = min(sm$ess_bulk, na.rm = TRUE),
  ess_tail_min = min(sm$ess_tail, na.rm = TRUE),
  n_divergent  = sum(diag_sum$num_divergent),
  pass         = max(sm$rhat, na.rm = TRUE) <= 1.01 &&
                 sum(diag_sum$num_divergent) == 0 &&
                 min(sm$ess_bulk, na.rm = TRUE) >= 100
)
cat(sprintf("STRICT GATE: rhat_max=%.4f | ess_bulk_min=%.0f | divergences=%d | PASS=%s\n",
            gate$rhat_max, gate$ess_bulk_min, gate$n_divergent, gate$pass))

# ---- Draws matrix (beta[1:3], theta[1:51]) ----------------------------------
dm <- as_draws_matrix(fit$draws(c("beta", "theta")))
beta_cols  <- paste0("beta[", 1:d$p, "]")
theta_cols <- paste0("theta[", 1:d$J, "]")
draws <- as.matrix(dm[, c(beta_cols, theta_cols)])
sigma_theta_hat <- mean(as_draws_matrix(fit$draws("sigma_theta")))
cat(sprintf("sigma_theta plug-in = %.4f\n", sigma_theta_hat))

pt_fe <- c("fe_between", "fe_within", "fe_between")

# ---- Dual-target DER (svyder 0.2.0 canonical) --------------------------------
fit_psu <- der_compute(
  draws, y = d$y, X = d$X, group = d$group, weights = w_v2,
  cluster = a$psu_idx, strata = a$stratum_idx,
  family = "binomial", sigma_theta = sigma_theta_hat,
  normalize = "unit_mean", param_types = pt_fe
)
fit_grp <- der_compute(
  draws, y = d$y, X = d$X, group = d$group, weights = w_v2,
  cluster = d$group,
  family = "binomial", sigma_theta = sigma_theta_hat,
  normalize = "unit_mean", param_types = pt_fe
)

cls_psu <- der_classify(fit_psu, tau = 1.2, verbose = FALSE)
cls_grp <- der_classify(fit_grp, tau = 1.2, verbose = FALSE)
cor_psu <- der_correct(cls_psu, method = "block_cholesky")
cor_grp <- der_correct(cls_grp, method = "block_cholesky")

# ---- Headline summary ---------------------------------------------------------
ci <- function(x, p = 0.90) quantile(x, c((1 - p) / 2, 1 - (1 - p) / 2))
b2_naive <- ci(draws[, 2]); b2_psu <- ci(cor_psu$corrected_draws[, 2])
cat("\n===== NSECE dual-target headline =====\n")
cat(sprintf("beta (posterior means): %s\n",
            paste(sprintf("%.3f", colMeans(draws[, 1:3])), collapse = " / ")))
cat(sprintf("PSU target  : beta2 DER = %.3f | flagged = %d/54\n",
            fit_psu$der[2], sum(cls_psu$classification$flagged)))
cat(sprintf("Group target: beta2 DER = %.3f | flagged = %d/54\n",
            fit_grp$der[2], sum(cls_grp$classification$flagged)))
cat(sprintf("beta2 90%% CI naive [%.3f, %.3f] -> PSU-corrected [%.3f, %.3f] (width x%.2f)\n",
            b2_naive[1], b2_naive[2], b2_psu[1], b2_psu[2],
            diff(b2_psu) / diff(b2_naive)))

# ---- Persist -------------------------------------------------------------------
saveRDS(list(
  stan_data_meta = list(N = d$N, J = d$J, p = d$p,
                        convention = "global_unit_mean_raw",
                        kish_raw = kish_raw, seed = 20260705,
                        mcmc = "4x2000/2000, adapt_delta 0.95"),
  gate            = gate,
  draws           = draws,
  sigma_theta_hat = sigma_theta_hat,
  beta_hat        = colMeans(draws[, 1:3]),
  fit_psu = list(der = fit_psu$der, V_sand = fit_psu$V_sand,
                 H_obs = fit_psu$H_obs, J_cluster = fit_psu$J_cluster,
                 target = fit_psu$target,
                 classification = cls_psu$classification,
                 corrected_draws = cor_psu$corrected_draws,
                 deff_j = fit_psu$deff_j, B_j = fit_psu$B_j),
  fit_grp = list(der = fit_grp$der, V_sand = fit_grp$V_sand,
                 target = fit_grp$target,
                 classification = cls_grp$classification,
                 corrected_draws = cor_grp$corrected_draws),
  computed_at = Sys.time()
), file.path(out_dir, "nsece_v2_refit.rds"))

cat(sprintf("\nSaved: %s | elapsed %.1f min\n",
            file.path(out_dir, "nsece_v2_refit.rds"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
