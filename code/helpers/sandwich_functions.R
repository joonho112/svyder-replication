# =============================================================================
# sandwich_functions.R: Sandwich Variance Estimation Functions
# =============================================================================
#
# Purpose : Build the observed information matrix (H_obs), clustered score
#           outer-product matrix (J_cluster), and sandwich variance estimator
#           (V_sand) for hierarchical logistic regression models fitted with
#           survey weights.
# Paper   : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
#           A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
#           arXiv preprint.
# Author  : JoonHo Lee (jlee296@ua.edu)
# License : MIT
#
# Contents:
#   1. build_H_obs_logistic()        -- Observed information (bread)
#   2. build_J_cluster()             -- Clustered meat matrix (CANONICAL:
#                                       stratified, centered, DF-corrected;
#                                       uncentered form reachable via flags)
#   3. build_J_cluster_psu()         -- PSU-level convenience wrapper
#                                       (uncentered defaults)
#   4. compute_sandwich()            -- Full sandwich variance estimator
#   5. compute_mcmc_covariance()     -- Posterior covariance from MCMC draws
#   6. compute_group_deff()          -- Per-group Kish design effect
#   7. compute_shrinkage_factors()   -- Per-group shrinkage factors B_j
#
# This file is the SINGLE canonical home of these functions. It is sourced by
# the simulation pipeline (sim_03_postprocess.R) and the application pipeline
# (app_03_der_analysis.R). Do not re-define these function names elsewhere in
# the package; a single definition of each name avoids the risk of divergent
# copies being sourced.
#
# Key equations (Paper Section 2.2-2.3):
#
#   H_obs = -d^2 ell_w / d phi d phi^T  (observed information, Eq. 5)
#         = H_data + H_prior
#         = sum_i w_i q_i(1-q_i) phi_i phi_i^T + diag(prior precisions)
#
#   J_cluster (default; matches svyder R/sandwich.R .build_J_cluster):
#         J = sum_h [C_h / (C_h - 1)] sum_c (t_hc - tbar_h)(t_hc - tbar_h)^T
#         where t_hc is the score total of cluster c in stratum h, C_h the
#         number of clusters in stratum h, and tbar_h the stratum mean of
#         cluster totals (Binder 1983 / survey::svyrecvar convention).
#
#   J_cluster (uncentered form: strata = NULL, center = FALSE, df_correct = FALSE):
#         J = sum_c t_c t_c^T          (uncentered, single stratum, no DF)
#
#   V_sand = H_obs^{-1} J_cluster H_obs^{-1}   (sandwich, Eq. 7)
#
# Parameter ordering: phi = (beta_1, ..., beta_p, theta_1, ..., theta_J)
# with dimension d = p + J.
#
# =============================================================================

library(Matrix)


###############################################################################
## Section 1 : Observed Information Matrix (H_obs)
###############################################################################

#' Build the observed information matrix for hierarchical logistic regression
#'
#' Constructs H_obs using efficient block structure (Paper Eq. 5):
#'
#'   H_obs = | H_bb   H_bt |   +   | diag(1/sigma_beta^2)   0                |
#'           | H_bt^T H_tt |       | 0                       diag(1/sigma_theta^2) |
#'
#' where:
#'   H_bb[a,b] = sum_i v_i * X[i,a] * X[i,b]       (p x p)
#'   H_bt[a,j] = sum_{i in j} v_i * X[i,a]          (p x J)
#'   H_tt[j,j] = sum_{i in j} v_i                   (J x J diagonal)
#'
#' and v_i = w_i * q_i * (1 - q_i) is the combined survey-working weight.
#'
#' @param X             Numeric matrix (N x p). Design matrix for fixed effects.
#' @param group         Integer vector (N). Group indices (1, ..., J).
#' @param w             Numeric vector (N). Normalized survey weights.
#' @param wt            Numeric vector (N). Working weights: q_i * (1 - q_i).
#' @param p             Integer. Number of fixed effects.
#' @param J             Integer. Number of groups.
#' @param N             Integer. Number of observations.
#' @param sigma_theta   Numeric scalar. Posterior mean of sigma_theta.
#' @param beta_prior_sd Numeric scalar. Prior SD for beta (default 5).
#'   Set to Inf for flat prior.
#'
#' @return Numeric matrix (d x d) where d = p + J.
build_H_obs_logistic <- function(X, group, w, wt, p, J, N,
                                 sigma_theta, beta_prior_sd = 5) {
  d <- p + J
  H <- matrix(0, d, d)

  ## Combined weight: v_i = w_i * q_i * (1 - q_i)
  v <- w * wt

  ## --- Beta-beta block (p x p): H_bb = X^T diag(v) X ---
  ## This is the standard weighted Fisher information for the FE.
  H_bb <- crossprod(X, X * v)
  H[1:p, 1:p] <- H_bb

  ## --- Beta-theta and theta-theta blocks ---
  ## For each group j:
  ##   H_bt[, j] = sum_{i in j} v_i * X[i, ]    (p-vector)
  ##   H_tt[j, j] = sum_{i in j} v_i             (scalar, diagonal)
  ## Cross-group theta entries are zero (groups are independent conditional
  ## on beta and sigma_theta).
  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    v_j   <- v[idx_j]

    ## Beta-theta coupling: reflects the correlation between FE and RE
    ## for group j. Strong coupling indicates the FE parameter absorbs
    ## information from group j (relevant for DER decomposition).
    if (length(idx_j) == 1L) {
      bt_j <- v_j * X[idx_j, ]
    } else {
      bt_j <- colSums(X[idx_j, , drop = FALSE] * v_j)
    }
    H[1:p, p + j] <- bt_j
    H[p + j, 1:p] <- bt_j

    ## Theta-theta diagonal: effective sample size for group j
    H[p + j, p + j] <- sum(v_j)
  }

  ## --- Prior contributions ---

  ## Beta prior: N(0, beta_prior_sd^2) => precision = 1/beta_prior_sd^2
  if (is.finite(beta_prior_sd)) {
    beta_prior_prec <- 1 / beta_prior_sd^2
    for (k in seq_len(p)) {
      H[k, k] <- H[k, k] + beta_prior_prec
    }
  }

  ## Theta prior: N(0, sigma_theta^2) => precision = 1/sigma_theta^2
  ## This is the shrinkage contribution from the hierarchical model.
  theta_prior_prec <- 1 / sigma_theta^2
  for (j in seq_len(J)) {
    H[p + j, p + j] <- H[p + j, p + j] + theta_prior_prec
  }

  ## Verify symmetry (should hold by construction)
  stopifnot(max(abs(H - t(H))) < 1e-10)

  H
}


###############################################################################
## Section 2 : Clustered Meat Matrix (Canonical, Stratified / Centered / DF)
###############################################################################

#' Build the clustered score outer-product matrix (meat)
#'
#' Constructs the meat matrix from cluster-level score totals. This is an
#' exact port of the canonical implementation in the companion svyder
#' package (R/sandwich.R, .build_J_cluster): the general (default) form is
#' the stratified, centered, DF-corrected estimator
#'
#'   J = sum_h [C_h / (C_h - 1)] sum_c (t_hc - tbar_h)(t_hc - tbar_h)^T,
#'
#' where t_hc is the score total of cluster c in stratum h, C_h the number
#' of clusters in stratum h, and tbar_h the stratum mean of cluster totals.
#' Centering matters at the posterior-mean plug-in point because the
#' likelihood score totals do not sum to zero there (the prior gradient
#' does not vanish).
#'
#' With strata = NULL, center = FALSE, df_correct = FALSE this reduces to
#' the uncentered single-stratum form J = sum_c t_c t_c^T (a documented
#' alternative to the default estimator).
#'
#' Strata containing a single cluster cannot be centered; their uncentered
#' contribution is used with a warning.
#'
#' The cluster-level score total t_c has the block-sparse structure
#'   t_c[1:p]   = sum_{i in c} r_i * X[i, ]                (beta block)
#'   t_c[p + j] = sum_{i in c, group = j} r_i              (theta block)
#' so a cluster spanning several model groups (e.g., a design PSU crossing
#' states) loads every theta coordinate it overlaps.
#'
#' @param X       Numeric matrix (N x p). Design matrix.
#' @param group   Integer vector (N). Model group indices (1, ..., J).
#' @param r       Numeric vector (N). Weighted score residuals:
#'   r_i = w_i * (y_i - q_i).
#' @param p       Integer. Number of fixed effects.
#' @param J       Integer. Number of model groups.
#' @param N       Integer. Number of observations.
#' @param cluster Integer vector (N). Sandwich aggregation unit
#'   (design PSU or model group), indices 1, ..., G. Defaults to
#'   \code{group} (model-group clustering, as in the simulation and the
#'   state-level application target).
#' @param strata  Optional integer vector (N). Stratum indicator
#'   (1, ..., H). NULL (default) treats the design as a single stratum.
#'   Each cluster must lie in exactly one stratum.
#' @param center  Logical. Center cluster score totals within strata
#'   (default TRUE; FALSE gives the uncentered form).
#' @param df_correct Logical. Apply the C_h / (C_h - 1) stratum correction
#'   (default TRUE; FALSE omits the correction).
#'
#' @return Numeric matrix (d x d) where d = p + J. Positive semi-definite
#'   by construction.
build_J_cluster <- function(X, group, r, p, J, N,
                            cluster = group, strata = NULL,
                            cluster_strata = NULL,
                            center = TRUE, df_correct = TRUE) {

  ## cluster_strata (optional): stratum of EVERY cluster in the declared
  ## design universe, indexed by cluster id 1..G. Supplying it declares the
  ## universe explicitly, so cluster ids that appear in no observation
  ## (selected-but-empty PSUs under Poisson subsampling) remain in the
  ## meat as zero score-total clusters -- they enter the stratum centering
  ## and the C_h/(C_h-1) correction. Without it, the universe is inferred
  ## from the observed ids (1..max(cluster)), which silently drops empty
  ## clusters.

  d <- p + J
  G <- if (!is.null(cluster_strata)) length(cluster_strata) else max(cluster)

  stopifnot(
    nrow(X) == N,
    length(group)   == N,
    length(r)       == N,
    length(cluster) == N,
    all(group >= 1L & group <= J),
    all(cluster >= 1L & cluster <= G),
    is.null(strata) || length(strata) == N
  )

  ## --- Cluster-to-stratum map (each cluster must lie in one stratum) ---
  if (!is.null(cluster_strata)) {
    if (is.null(strata)) {
      stop("'cluster_strata' requires 'strata'.", call. = FALSE)
    }
    strata_of_cluster <- as.integer(cluster_strata)
    bad <- strata != strata_of_cluster[cluster]
    if (any(bad)) {
      stop(sprintf(
        "'strata' disagrees with 'cluster_strata' at %d observation(s).",
        sum(bad)), call. = FALSE)
    }
  } else if (is.null(strata)) {
    strata_of_cluster <- rep(1L, G)
  } else {
    tab <- unique(data.frame(cluster = cluster, strata = strata))
    if (nrow(tab) != G) {
      stop("Each cluster must belong to exactly one stratum.", call. = FALSE)
    }
    strata_of_cluster <- integer(G)
    strata_of_cluster[tab$cluster] <- tab$strata
  }

  ## --- Cluster-level score totals: T_mat is G x d ---
  T_mat <- matrix(0, G, d)
  for (g in seq_len(G)) {
    idx_g <- which(cluster == g)
    if (length(idx_g) == 0L) next

    s_g <- numeric(d)

    ## Beta score component
    if (length(idx_g) == 1L) {
      s_g[1:p] <- r[idx_g] * X[idx_g, ]
    } else {
      s_g[1:p] <- colSums(X[idx_g, , drop = FALSE] * r[idx_g])
    }

    ## Theta score component: accumulate within each group present in
    ## cluster g (allows non-nested clusters spanning several groups)
    for (j_grp in unique(group[idx_g])) {
      idx_gj <- idx_g[group[idx_g] == j_grp]
      s_g[p + j_grp] <- s_g[p + j_grp] + sum(r[idx_gj])
    }

    T_mat[g, ] <- s_g
  }

  ## --- Stratum-wise (centered, DF-corrected) accumulation ---
  J_cluster <- matrix(0, d, d)
  singleton_strata <- 0L

  for (h in unique(strata_of_cluster)) {
    idx_h <- which(strata_of_cluster == h)
    C_h   <- length(idx_h)
    T_h   <- T_mat[idx_h, , drop = FALSE]

    if (center && C_h > 1L) {
      T_h <- sweep(T_h, 2L, colMeans(T_h))
    } else if (center && C_h == 1L) {
      singleton_strata <- singleton_strata + 1L
    }

    mult <- if (df_correct && C_h > 1L) C_h / (C_h - 1) else 1
    J_cluster <- J_cluster + mult * crossprod(T_h)
  }

  if (singleton_strata > 0L) {
    warning(sprintf(
      "%d stratum/strata contain a single cluster; their contribution is uncentered.",
      singleton_strata
    ), call. = FALSE)
  }

  ## Symmetrize
  J_cluster <- (J_cluster + t(J_cluster)) / 2

  J_cluster
}


###############################################################################
## Section 3 : Clustered Meat Matrix (PSU-Level, Back-Compatibility Wrapper)
###############################################################################

#' Build PSU-level clustered score outer-product matrix (convenience wrapper)
#'
#' Convenience wrapper around the canonical build_J_cluster() for the case
#' where the sandwich aggregation unit is the design PSU rather than the
#' model group (e.g., NSECE: 415 PSUs vs. 51 states). The wrapper's DEFAULTS
#' give the uncentered single-stratum estimator (no DF correction). Pass
#' strata / center / df_correct to obtain the stratified, centered,
#' DF-corrected estimator; new code should call
#' build_J_cluster(cluster = psu, ...) directly.
#'
#' @param X        Numeric matrix (N x p). Design matrix.
#' @param state    Integer vector (N). State (model group) indices
#'   (1, ..., J_states).
#' @param psu      Integer vector (N). PSU indices (1, ..., G).
#' @param r        Numeric vector (N). Weighted score residuals.
#' @param p        Integer. Number of fixed effects.
#' @param J_states Integer. Number of states (random-effect groups).
#' @param N        Integer. Number of observations.
#' @param strata   Optional integer vector (N). Stratum indicator; see
#'   build_J_cluster(). Default NULL (single stratum).
#' @param center   Logical. Default FALSE (uncentered).
#' @param df_correct Logical. Default FALSE (no DF correction).
#'
#' @return Numeric matrix (d x d) where d = p + J_states.
build_J_cluster_psu <- function(X, state, psu, r, p, J_states, N,
                                strata = NULL, center = FALSE,
                                df_correct = FALSE) {
  build_J_cluster(
    X = X, group = state, r = r, p = p, J = J_states, N = N,
    cluster = psu, strata = strata,
    center = center, df_correct = df_correct
  )
}


###############################################################################
## Section 4 : Full Sandwich Variance Estimator
###############################################################################

#' Compute the full sandwich variance estimator
#'
#' Assembles the sandwich variance V_sand = H_obs^{-1} J_cluster H_obs^{-1}
#' (Paper Eq. 7). Includes a nearPD fallback for near-singular H_obs.
#'
#' @param H_obs     Numeric matrix (d x d). Observed information matrix.
#' @param J_cluster Numeric matrix (d x d). Clustered score outer-product.
#'
#' @return A list with components:
#'   \describe{
#'     \item{V_sand}{Numeric matrix (d x d). Sandwich variance estimator.}
#'     \item{H_obs_inv}{Numeric matrix (d x d). Inverse of H_obs.}
#'   }
compute_sandwich <- function(H_obs, J_cluster) {

  stopifnot(
    is.matrix(H_obs),
    is.matrix(J_cluster),
    nrow(H_obs) == ncol(H_obs),
    nrow(J_cluster) == ncol(J_cluster),
    nrow(H_obs) == nrow(J_cluster)
  )

  d <- nrow(H_obs)

  ## Invert H_obs (with nearPD fallback for numerical stability)
  H_obs_inv <- tryCatch(
    solve(H_obs),
    error = function(e) {
      warning("H_obs is singular or near-singular. ",
              "Using nearPD fallback. ",
              "Condition number: ",
              tryCatch(kappa(H_obs, exact = TRUE),
                       error = function(e2) NA_real_))
      H_pd <- as.matrix(Matrix::nearPD(H_obs, keepDiag = TRUE)$mat)
      solve(H_pd)
    }
  )

  ## Sandwich: V_sand = H_inv J H_inv
  V_sand <- H_obs_inv %*% J_cluster %*% H_obs_inv

  ## Enforce symmetry
  V_sand <- (V_sand + t(V_sand)) / 2

  list(
    V_sand    = V_sand,
    H_obs_inv = H_obs_inv
  )
}


###############################################################################
## Section 5 : MCMC Posterior Covariance
###############################################################################

#' Compute posterior covariance matrix from MCMC draws
#'
#' Extracts posterior draws for all parameters (beta and theta), merges
#' across chains, and computes the sample covariance matrix Sigma_MCMC.
#' This serves as the denominator in the DER computation (Paper Eq. 8).
#'
#' @param fit A cmdstanr fit object with a draws() method.
#' @param p   Integer. Number of fixed effects.
#' @param J   Integer. Number of groups (random effects).
#'
#' @return Numeric matrix (d x d) where d = p + J. The posterior covariance
#'   matrix Sigma_MCMC = Cov(phi | y).
compute_mcmc_covariance <- function(fit, p, J) {
  d <- p + J

  ## Parameter names in Stan
  beta_names  <- paste0("beta[", seq_len(p), "]")
  theta_names <- paste0("theta[", seq_len(J), "]")
  all_names   <- c(beta_names, theta_names)

  ## Extract draws as matrix (iterations x parameters), merging all chains
  draws <- fit$draws(variables = all_names, format = "matrix")

  stopifnot(
    is.matrix(draws),
    ncol(draws) == d
  )

  ## Compute sample covariance
  sigma_mcmc <- cov(draws)

  ## Enforce exact symmetry
  sigma_mcmc <- (sigma_mcmc + t(sigma_mcmc)) / 2

  sigma_mcmc
}


###############################################################################
## Section 6 : Per-Group Design Effect (DEFF)
###############################################################################

#' Compute per-group Kish design effect (DEFF)
#'
#' For each group j, computes the Kish (1965) design effect:
#'
#'   DEFF_j = n_j * sum_{i in j} w_i^2 / (sum_{i in j} w_i)^2
#'
#' This measures the variance inflation due to unequal weighting within
#' each group. DEFF = 1 for equal weights; DEFF > 1 for unequal weights.
#' The DEFF enters the DER decomposition (Paper Theorem 2, Eq. 12).
#'
#' @param group Integer vector (N). Group indices (1, ..., J).
#' @param w     Numeric vector (N). Survey weights.
#' @param J     Integer. Number of groups.
#'
#' @return Numeric vector of length J with per-group DEFF values (all >= 1).
compute_group_deff <- function(group, w, J) {
  deff_j <- numeric(J)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)
    w_j   <- w[idx_j]
    n_j   <- length(idx_j)

    ## Kish DEFF: n * sum(w^2) / (sum(w))^2
    deff_j[j] <- n_j * sum(w_j^2) / (sum(w_j))^2
  }

  deff_j
}


###############################################################################
## Section 7 : Shrinkage Factors
###############################################################################

#' Compute per-group shrinkage factors B_j
#'
#' The shrinkage factor B_j measures the data/reliability share for group j's
#' random effect under the hierarchical model (Paper Eq. 9):
#'
#'   B_j = sigma_theta^2 / (sigma_theta^2 + V_tilde_j)
#'
#' where V_tilde_j = 1 / sum_{i in j} w_i * q_i * (1 - q_i) is the
#' effective conditional variance from group j's data.
#'
#' B_j close to 1 means data-dominated / little prior shrinkage.
#' B_j close to 0 means prior-dominated / heavy shrinkage.
#'
#' The shrinkage factor appears in the DER decomposition:
#'   - For RE: DER_j = B_j * DEFF_j * kappa_j (Theorem 2)
#'   - For between-FE: DER_beta = DEFF * (1 - B_bar) (Theorem 1)
#'
#' @param group        Integer vector (N). Group indices (1, ..., J).
#' @param w            Numeric vector (N). Survey weights.
#' @param wt           Numeric vector (N). Working weights: q_i * (1 - q_i).
#' @param sigma_theta  Numeric scalar. Posterior mean of sigma_theta.
#' @param J            Integer. Number of groups.
#'
#' @return Numeric vector of length J with shrinkage factors in (0, 1).
compute_shrinkage_factors <- function(group, w, wt, sigma_theta, J) {
  sigma2 <- sigma_theta^2
  B_j <- numeric(J)

  for (j in seq_len(J)) {
    idx_j <- which(group == j)

    ## Effective information from group j
    info_j <- sum(w[idx_j] * wt[idx_j])

    ## Conditional variance: inverse of effective information
    V_tilde_j <- 1 / info_j

    ## Shrinkage factor (Paper Eq. 9)
    B_j[j] <- sigma2 / (sigma2 + V_tilde_j)
  }

  stopifnot(all(B_j > 0), all(B_j < 1))

  B_j
}


###############################################################################
## Section 8 : Self-Test
###############################################################################

if (interactive()) {
  cat("=== sandwich_functions.R: Self-test ===\n\n")

  set.seed(123)
  N_test <- 100
  J_test <- 5
  p_test <- 3

  ## Synthetic data
  X_test     <- cbind(1, rnorm(N_test), rep(rnorm(J_test), each = N_test / J_test))
  group_test <- rep(1:J_test, each = N_test / J_test)
  w_test     <- runif(N_test, 0.5, 2.0)
  wt_test    <- runif(N_test, 0.1, 0.25)

  ## Test H_obs construction
  H_test <- build_H_obs_logistic(X_test, group_test, w_test, wt_test,
                                 p_test, J_test, N_test,
                                 sigma_theta = 0.5, beta_prior_sd = 5)
  d_test <- p_test + J_test
  stopifnot(
    nrow(H_test) == d_test,
    ncol(H_test) == d_test,
    max(abs(H_test - t(H_test))) < 1e-12
  )
  cat("  build_H_obs_logistic: OK (", d_test, "x", d_test, ")\n")

  ## Test J_cluster construction (default: centered + DF-corrected)
  r_test <- rnorm(N_test)
  Jmat_test <- build_J_cluster(X_test, group_test, r_test,
                               p_test, J_test, N_test)
  stopifnot(
    nrow(Jmat_test) == d_test,
    ncol(Jmat_test) == d_test,
    max(abs(Jmat_test - t(Jmat_test))) < 1e-12
  )
  eigen_vals <- eigen(Jmat_test, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(all(eigen_vals >= -1e-10))
  cat("  build_J_cluster (default): OK (PSD verified)\n")

  ## Uncentered flags reproduce the sum of outer products exactly
  Jmat_v1 <- build_J_cluster(X_test, group_test, r_test,
                             p_test, J_test, N_test,
                             strata = NULL, center = FALSE,
                             df_correct = FALSE)
  T_ref <- matrix(0, J_test, d_test)
  for (j in seq_len(J_test)) {
    idx <- which(group_test == j)
    T_ref[j, 1:p_test] <- colSums(X_test[idx, , drop = FALSE] * r_test[idx])
    T_ref[j, p_test + j] <- sum(r_test[idx])
  }
  J_ref_v1 <- crossprod(T_ref)
  stopifnot(max(abs(Jmat_v1 - J_ref_v1)) < 1e-10)
  cat("  build_J_cluster (legacy flags): OK (matches uncentered form)\n")

  ## Single-stratum default = (G/(G-1)) * centered cross-product
  T_cen <- sweep(T_ref, 2L, colMeans(T_ref))
  J_ref_v2 <- (J_test / (J_test - 1)) * crossprod(T_cen)
  stopifnot(max(abs(Jmat_test - J_ref_v2)) < 1e-10)
  cat("  build_J_cluster (default): OK (matches centered, DF-corrected form)\n")

  ## Stratified accumulation: 2 strata of groups, hand-computed reference
  strata_of_group <- rep(c(1L, 2L), length.out = J_test)
  strata_test <- strata_of_group[group_test]
  J_strat <- build_J_cluster(X_test, group_test, r_test,
                             p_test, J_test, N_test,
                             strata = strata_test,
                             center = TRUE, df_correct = TRUE)
  J_ref_strat <- matrix(0, d_test, d_test)
  for (h in 1:2) {
    idx_h <- which(strata_of_group == h)
    T_h <- sweep(T_ref[idx_h, , drop = FALSE], 2L,
                 colMeans(T_ref[idx_h, , drop = FALSE]))
    C_h <- length(idx_h)
    J_ref_strat <- J_ref_strat + (C_h / (C_h - 1)) * crossprod(T_h)
  }
  stopifnot(max(abs(J_strat - J_ref_strat)) < 1e-10)
  cat("  build_J_cluster (stratified): OK (matches hand computation)\n")

  ## Singleton stratum triggers the uncentered-contribution warning
  strata_singleton <- ifelse(group_test == 1L, 1L, 2L)
  warn_seen <- FALSE
  withCallingHandlers(
    build_J_cluster(X_test, group_test, r_test, p_test, J_test, N_test,
                    strata = strata_singleton,
                    center = TRUE, df_correct = TRUE),
    warning = function(w) {
      warn_seen <<- grepl("single cluster", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  stopifnot(warn_seen)
  cat("  build_J_cluster (singleton stratum): OK (warning emitted)\n")

  ## PSU wrapper with uncentered defaults = uncentered group-level result
  ## when psu == group
  Jmat_psu <- build_J_cluster_psu(X_test, group_test, group_test, r_test,
                                  p_test, J_test, N_test)
  stopifnot(max(abs(Jmat_psu - J_ref_v1)) < 1e-10)
  cat("  build_J_cluster_psu (legacy wrapper): OK\n")

  ## Test compute_sandwich
  sand_test <- compute_sandwich(H_test, Jmat_test)
  stopifnot(
    nrow(sand_test$V_sand) == d_test,
    ncol(sand_test$V_sand) == d_test
  )
  cat("  compute_sandwich: OK\n")

  ## Test compute_group_deff
  deff_test <- compute_group_deff(group_test, w_test, J_test)
  stopifnot(
    length(deff_test) == J_test,
    all(deff_test >= 1)
  )
  cat("  compute_group_deff: OK (all >= 1)\n")

  ## Test compute_shrinkage_factors
  B_test <- compute_shrinkage_factors(group_test, w_test, wt_test,
                                      sigma_theta = 0.5, J_test)
  stopifnot(
    length(B_test) == J_test,
    all(B_test > 0),
    all(B_test < 1)
  )
  cat("  compute_shrinkage_factors: OK (all in (0,1))\n")

  cat("\nAll self-tests passed.\n")
}
