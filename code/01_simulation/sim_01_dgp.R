# =============================================================================
# sim_01_dgp.R: Data-Generating Mechanism with a genuine two-stage design
# Author : JoonHo Lee (jlee296@ua.edu)
# =============================================================================
#
# Data-generating mechanism: a genuine two-stage informative sampling design.
# Sourced after sim_00_config.R, it provides generate_survey_data(), which the
# fitting and post-processing stages consume. Key properties:
#
#   1. A finite population is generated first.
#   2. Samples are drawn by an explicit two-stage mechanism whose inclusion
#      probabilities depend on the random effects (stage 1) and the response
#      (stage 2) when informativeness is on.
#   3. Survey weights are EXACT inverse first-order inclusion probabilities
#      of the design (capping and certainty units are reflected in the
#      probabilities actually used; empty Poisson samples are retained,
#      never redrawn -- see poisson_sample() and the replicate guard).
#   4. Normalization is GLOBAL unit-mean (the declared convention), applied
#      after the weights exist -- never within-cluster rescaling, which
#      cancels cluster-level informative components exactly.
#
# Two structural arms (config: sim_00_config.R):
#   * sampled_groups   : stage 1 = PPS(exp(gamma_grp*theta_g + e1)) sampling
#                        of J from M population groups; stage 2 = Poisson
#                        (exp(gamma_unit*y + e2)) within groups.
#                        Design PSU = model group.
#   * psu_within_groups: all J groups observed; stage 1 = PPS sampling of
#                        k of K PSUs within each group (sizes depend on a
#                        PSU-level outcome shock upsilon_c that the analysis
#                        model omits); stage 2 = Poisson within PSUs.
#                        Design PSU != model group; PSU stage stratified by
#                        group.
#
# All functions are base R. Sourcing this file after sim_00_config.R
# provides generate_survey_data(); dgp_acceptance_test() implements the
# acceptance checks on the analysis weights.
# =============================================================================


# =============================================================================
# Section 1: Sampling primitives with exact first-order probabilities
# =============================================================================

#' Systematic PPS sampling with certainty-unit handling
#'
#' Selects n units with probability proportional to size. Units whose
#' implied probability exceeds 1 are taken with certainty and the remainder
#' re-scaled (standard iterative procedure). Selection is systematic on a
#' randomly permuted list, so first-order probabilities are exactly
#' pi_i = n' * s_i / S' on the non-certainty set.
#'
#' @param sizes Positive size measures (length M).
#' @param n Number of units to select.
#' @return list(idx = selected indices, pi = exact first-order inclusion
#'   probabilities for ALL M units)
pps_systematic <- function(sizes, n) {
  stopifnot(all(sizes > 0), n >= 1, n <= length(sizes))
  M  <- length(sizes)
  pi <- rep(NA_real_, M)

  certain <- logical(M)
  repeat {
    rem   <- which(!certain)
    n_rem <- n - sum(certain)
    if (n_rem == 0L) break
    p_rem <- n_rem * sizes[rem] / sum(sizes[rem])
    over  <- p_rem >= 1
    if (!any(over)) {
      pi[rem] <- p_rem
      break
    }
    certain[rem[over]] <- TRUE
  }
  pi[certain] <- 1

  # Systematic selection on the non-certainty set, random order
  sel <- which(certain)
  rem <- which(!certain)
  n_rem <- n - length(sel)
  if (n_rem > 0L) {
    ord   <- sample(rem)
    cum   <- cumsum(pi[ord])          # total = n_rem
    start <- runif(1, 0, 1)
    pts   <- start + seq_len(n_rem) - 1
    take  <- findInterval(pts, c(0, cum), rightmost.closed = TRUE)
    sel   <- c(sel, ord[unique(take)])
    # Numerical guard: systematic PPS yields exactly n_rem distinct hits;
    # duplicated interval indices cannot occur when all pi < 1.
    stopifnot(length(sel) == n)
  }

  list(idx = sort(sel), pi = pi)
}


#' Poisson sampling with capped probabilities (exact first-order design)
#'
#' Each unit is included independently with pi_i = min(cap,
#' n_target * s_i / sum(s)). The capped values are the probabilities
#' actually used, so weights 1/pi are EXACT first-order inverse
#' inclusion probabilities. Empty samples are a legitimate outcome of
#' Poisson sampling and are RETURNED AS EMPTY -- no redraw, no forced
#' inclusion (a redraw-until-nonempty rule would condition the design
#' and silently shift the effective first-order probabilities).
#' Degenerate whole-group outcomes are handled at the replicate level
#' (see generate_survey_data).
#'
#' @return list(idx (possibly empty), pi (all units))
poisson_sample <- function(sizes, n_target, cap = 0.95) {
  stopifnot(all(sizes > 0), n_target > 0)
  pi  <- pmin(cap, n_target * sizes / sum(sizes))
  idx <- which(runif(length(sizes)) < pi)
  list(idx = idx, pi = pi)
}


# =============================================================================
# Section 2: Population generators
# =============================================================================

#' Population for the sampled_groups arm
#'
#' @return list of population vectors; units stacked by group.
generate_population_A <- function(M, N_g, sigma_theta, beta) {
  theta_g <- rnorm(M, 0, sigma_theta)
  z_g     <- rnorm(M)                       # group-level covariate
  n_tot   <- M * N_g
  group   <- rep(seq_len(M), each = N_g)
  x       <- rnorm(n_tot)                   # unit-level covariate
  eta     <- beta[1] + beta[2] * x + beta[3] * z_g[group] + theta_g[group]
  y       <- rbinom(n_tot, 1L, plogis(eta))
  list(M = M, N_g = N_g, group = group, x = x, y = y,
       theta_g = theta_g, z_g = z_g)
}


#' Population for the psu_within_groups arm
#'
#' All J groups will be observed; K PSUs per group with a PSU-level shock
#' upsilon that the analysis model omits.
generate_population_B <- function(J, K, N_c, sigma_theta, sigma_psu, beta) {
  theta_g <- rnorm(J, 0, sigma_theta)
  z_g     <- rnorm(J)
  n_psu   <- J * K
  ups_c   <- rnorm(n_psu, 0, sigma_psu)
  n_tot   <- n_psu * N_c
  psu     <- rep(seq_len(n_psu), each = N_c)         # global PSU id
  group   <- rep(rep(seq_len(J), each = K), each = N_c)
  x       <- rnorm(n_tot)
  eta     <- beta[1] + beta[2] * x + beta[3] * z_g[group] +
             theta_g[group] + ups_c[psu]
  y       <- rbinom(n_tot, 1L, plogis(eta))
  list(J = J, K = K, N_c = N_c, group = group, psu = psu,
       x = x, y = y, theta_g = theta_g, z_g = z_g, ups_c = ups_c)
}


# =============================================================================
# Section 3: Two-stage samplers (exact inverse-probability weights)
# =============================================================================

#' Two-stage sample for the sampled_groups arm
sample_two_stage_A <- function(pop, J, nbar, gamma_grp, gamma_unit,
                               tau1, tau2) {

  # ---- Stage 1: PPS on groups, sizes depend on theta_g -----------------------
  s_grp  <- exp(gamma_grp * pop$theta_g + rnorm(pop$M, 0, tau1))
  st1    <- pps_systematic(s_grp, J)
  g_sel  <- st1$idx

  # ---- Stage 2: Poisson within each sampled group, sizes depend on y ---------
  keep_idx <- integer(0)
  pi2      <- numeric(0)

  for (g in g_sel) {
    u_g   <- which(pop$group == g)
    s_uni <- exp(gamma_unit * pop$y[u_g] + rnorm(length(u_g), 0, tau2))
    st2   <- poisson_sample(s_uni, n_target = nbar)
    keep_idx <- c(keep_idx, u_g[st2$idx])
    pi2      <- c(pi2, st2$pi[st2$idx])
  }

  # ---- Assemble the sample frame ---------------------------------------------
  g_old  <- pop$group[keep_idx]
  g_new  <- match(g_old, g_sel)              # relabel sampled groups 1..J
  w_raw  <- 1 / (st1$pi[g_old] * pi2)

  list(
    y       = pop$y[keep_idx],
    x       = pop$x[keep_idx],
    z       = pop$z_g[g_sel][g_new],
    group   = g_new,
    psu     = g_new,                          # design PSU = model group
    strata  = NULL,                           # single stratum
    w_raw   = w_raw,
    truth   = list(theta = pop$theta_g[g_sel], z_g = pop$z_g[g_sel]),
    diag    = list(pi1 = st1$pi[g_sel], f1 = J / pop$M)
  )
}


#' Two-stage sample for the psu_within_groups arm
#'
#' Groups are all observed (pi_g = 1); k of K PSUs sampled per group with
#' sizes depending on the omitted PSU shock; Poisson within PSUs.
sample_two_stage_B <- function(pop, k, nbar, gamma_grp, gamma_unit,
                               tau1, tau2) {

  keep_idx <- integer(0)
  pi_all   <- numeric(0)
  psu_sel_global <- integer(0)

  for (g in seq_len(pop$J)) {
    psu_g <- ((g - 1) * pop$K + 1):(g * pop$K)      # global PSU ids in group g

    # ---- Stage 1: PPS on PSUs within group g (stratified by group) ----------
    s_psu <- exp(gamma_grp * pop$ups_c[psu_g] + rnorm(pop$K, 0, tau1))
    st1   <- pps_systematic(s_psu, k)
    c_sel <- psu_g[st1$idx]
    psu_sel_global <- c(psu_sel_global, c_sel)

    # ---- Stage 2: Poisson within each sampled PSU ----------------------------
    for (ci in seq_along(c_sel)) {
      u_c   <- which(pop$psu == c_sel[ci])
      s_uni <- exp(gamma_unit * pop$y[u_c] + rnorm(length(u_c), 0, tau2))
      st2   <- poisson_sample(s_uni, n_target = nbar / k)
      keep_idx <- c(keep_idx, u_c[st2$idx])
      pi_all   <- c(pi_all, st1$pi[st1$idx[ci]] * st2$pi[st2$idx])
    }
  }

  # Design-PSU cluster UNIVERSE = ALL selected PSUs (J*k), including any
  # whose second-stage Poisson sample came back empty. A selected-but-empty
  # PSU is a legitimate design outcome whose estimated score total is
  # exactly zero; it must remain in the cluster universe so that the
  # stratified meat centers over, DF-corrects with, and sums over all J*k
  # selected clusters. Dropping empties would silently change the declared
  # design-PSU target.
  psu_map <- sort(psu_sel_global)              # universe id j <-> global id psu_map[j]
  psu_new <- match(pop$psu[keep_idx], psu_map)
  psu_strata_map <- (psu_map - 1L) %/% pop$K + 1L  # stratum (= group) of EVERY selected PSU

  list(
    y       = pop$y[keep_idx],
    x       = pop$x[keep_idx],
    z       = pop$z_g[pop$group[keep_idx]],
    group   = pop$group[keep_idx],            # all J groups present
    psu     = psu_new,                        # design PSU id in the selected universe (1..J*k)
    strata  = pop$group[keep_idx],            # PSU stage stratified by group
    psu_strata_map = psu_strata_map,          # stratum of each universe PSU, incl. empty ones
    w_raw   = 1 / pi_all,
    truth   = list(theta = pop$theta_g, z_g = pop$z_g,
                   ups_c = pop$ups_c, psu_map = psu_map),
    diag    = list(f1 = k / pop$K)
  )
}


# =============================================================================
# Section 4: Orchestrator
# =============================================================================

#' Generate one replication data set
#'
#' @param cell One row of build_scenario_grid().
#' @param seed Integer seed for this replication.
#' @return list(y, X, group, psu, strata, weights, w_raw, truth, diagnostics)
#'
#' Degenerate-replicate guard: under Poisson stage-2 sampling, a model
#' group can (rarely) end up with ZERO sampled units -- e.g., in the
#' smallest psu_within_groups cells the per-group probability is about
#' (0.072)^4 = 2.7e-5. A group with no data breaks the fitted model's
#' indexing, so such replicates are REGENERATED with a deterministic
#' seed offset (+1,000,000 per attempt). This is a replicate-level
#' rejection rule, applied identically across arms; the conditioning it
#' induces on the design is of order 1e-5 and is disclosed in the
#' manuscript. Within accepted replicates the weights are the exact
#' unconditional first-order inverse inclusion probabilities.
generate_survey_data <- function(cell, seed) {

  P <- SIM_PARAMS
  attempt   <- 0L
  seed_used <- as.integer(seed)

  repeat {
    set.seed(seed_used)

    if (cell$structure == "sampled_groups") {
      pop <- generate_population_A(
        M = P$M_groups, N_g = P$N_g_factor * cell$nbar_j,
        sigma_theta = cell$sigma_theta, beta = P$beta_true
      )
      smp <- sample_two_stage_A(
        pop, J = cell$J, nbar = cell$nbar_j,
        gamma_grp = cell$gamma_grp, gamma_unit = cell$gamma_unit,
        tau1 = P$tau1, tau2 = P$tau2
      )
    } else {
      pop <- generate_population_B(
        J = cell$J, K = P$K_psu,
        N_c = P$N_c_factor * cell$nbar_j / P$k_psu,
        sigma_theta = cell$sigma_theta, sigma_psu = cell$sigma_psu,
        beta = P$beta_true
      )
      smp <- sample_two_stage_B(
        pop, k = P$k_psu, nbar = cell$nbar_j,
        gamma_grp = cell$gamma_grp, gamma_unit = cell$gamma_unit,
        tau1 = P$tau1, tau2 = P$tau2
      )
    }

    # Accept iff every model group contributed at least one observation
    if (length(unique(smp$group)) == cell$J) break
    attempt   <- attempt + 1L
    if (attempt > 5L) stop("degenerate-replicate guard exceeded 5 attempts")
    seed_used <- seed_used + 1000000L
  }

  # ---- Declared weight convention: global unit-mean --------------------------
  w <- smp$w_raw / mean(smp$w_raw)

  # ---- Within-group centering of the unit covariate --------------------------
  x_wc <- smp$x - ave(smp$x, smp$group)

  X <- cbind(intercept = 1, x_wc = x_wc, z_bc = smp$z)

  # ---- Realized-design diagnostics (measured, never asserted) ----------------
  n_j <- as.integer(table(factor(smp$group, levels = sort(unique(smp$group)))))
  deff_global <- length(w) * sum(w^2) / sum(w)^2

  list(
    y = smp$y, X = X, group = smp$group, psu = smp$psu,
    strata = smp$strata, psu_strata_map = smp$psu_strata_map,
    weights = w, w_raw = smp$w_raw,
    truth = c(smp$truth,
              list(beta = P$beta_true, sigma_theta = cell$sigma_theta)),
    diagnostics = list(
      N = length(w), n_j = n_j,
      deff_kish_global = deff_global,
      cv_w = sd(w) / mean(w),
      n_reseed = attempt,
      f1 = smp$diag$f1,
      n_psu = length(unique(smp$psu)),
      n_psu_universe = if (!is.null(smp$psu_strata_map))
        length(smp$psu_strata_map) else length(unique(smp$psu)),
      n_psu_empty = if (!is.null(smp$psu_strata_map))
        length(smp$psu_strata_map) - length(unique(smp$psu)) else 0L,
      seed = seed, seed_used = seed_used
    )
  )
}


# =============================================================================
# Section 5: Acceptance tests on the analysis weights
# =============================================================================
#
# Within-cluster normalization would cancel the informative component, so the
# declared convention is global unit-mean. These tests verify, on the ANALYSIS
# weights (post-normalization):
#   (T1) informative arms carry a real association between the weights and
#        the STAGE-1 SELECTION SHOCK -- which is theta_g in the
#        sampled_groups arm but the omitted PSU shock upsilon_c in the
#        psu_within_groups arm (there, all groups are observed with
#        pi_g = 1, so weights are theta-free BY DESIGN); the none arms
#        carry no such association;
#   (T2) response-dependent selection is visible: E[w | y=1] < E[w | y=0]
#        in informative arms;
#   (T3) sampling primitives deliver their exact first-order probabilities
#        (Monte Carlo check);
#   (T4) realized DEFF separates the gamma levels.
# =============================================================================

dgp_acceptance_test <- function(n_rep = 40L, verbose = TRUE) {

  grid <- build_scenario_grid()
  out  <- list()

  # ---- T3: primitive probability checks --------------------------------------
  set.seed(1)
  sizes <- exp(rnorm(40, 0, 0.8))
  pi_an <- pps_systematic(sizes, 10)$pi     # analytic (selection discarded)
  hits  <- numeric(40)
  for (r in 1:4000) hits <- hits + (seq_len(40) %in% pps_systematic(sizes, 10)$idx)
  t3_pps <- max(abs(hits / 4000 - pi_an))

  set.seed(2)
  sizes2 <- exp(rnorm(200, 0, 0.8))
  pi_an2 <- poisson_sample(sizes2, 20)$pi
  hits2  <- numeric(200)
  for (r in 1:4000) hits2 <- hits2 + (seq_len(200) %in% poisson_sample(sizes2, 20)$idx)
  t3_poi <- max(abs(hits2 / 4000 - pi_an2))

  # Small-cell Poisson check at the actual smallest design values
  # (N_c = 25, n_target = 2.5): empty samples must be retained so that
  # empirical inclusion frequencies match the analytic pi exactly.
  set.seed(3)
  sizes3 <- exp(rnorm(25, 0, 0.8))
  pi_an3 <- poisson_sample(sizes3, 2.5)$pi
  hits3  <- numeric(25); n_empty <- 0L
  for (r in 1:8000) {
    d3 <- poisson_sample(sizes3, 2.5)
    hits3 <- hits3 + (seq_len(25) %in% d3$idx)
    n_empty <- n_empty + (length(d3$idx) == 0L)
  }
  t3_poi_small <- max(abs(hits3 / 8000 - pi_an3))
  empty_rate   <- n_empty / 8000
  empty_theory <- prod(1 - pi_an3)

  # ---- T1/T2/T4 across gamma levels, both structures --------------------------
  for (struct in c("sampled_groups", "psu_within_groups")) {
    for (gam in c("none", "moderate", "strong")) {
      cell <- grid[grid$J == 20 & grid$nbar_j == 50 &
                   grid$informativeness == gam &
                   grid$structure == struct, ]
      cor_tw <- dw <- deff <- numeric(n_rep)
      for (r in seq_len(n_rep)) {
        d <- generate_survey_data(cell, seed = 55000 + r)
        if (struct == "sampled_groups") {
          # stage-1 shock = theta_g of the sampled groups
          lw_g <- tapply(log(d$weights), d$group, mean)
          sh   <- d$truth$theta[as.integer(names(lw_g))]
        } else {
          # stage-1 shock = omitted PSU shock upsilon_c of the sampled PSUs
          lw_g <- tapply(log(d$weights), d$psu, mean)
          sh   <- d$truth$ups_c[d$truth$psu_map[as.integer(names(lw_g))]]
        }
        cor_tw[r] <- suppressWarnings(cor(lw_g, sh))
        dw[r]     <- mean(d$weights[d$y == 1]) - mean(d$weights[d$y == 0])
        deff[r]   <- d$diagnostics$deff_kish_global
      }
      out[[paste(struct, gam, sep = ".")]] <- c(
        cor_shock_logw = mean(cor_tw, na.rm = TRUE),
        dw_y1_minus_y0 = mean(dw),
        deff_realized  = mean(deff)
      )
    }
  }

  # ---- T5: selected-but-empty PSUs stay in the design-PSU universe -----------
  # In the smallest PW cells the per-PSU empty-sample probability is about
  # 7%, so selected PSUs with zero sampled units are routine. They must
  # (i) remain in the returned universe
  # with their stratum and (ii) enter the stratified meat as zero score-total
  # clusters. (ii) is checked against a hand-rolled meat when the canonical
  # builder is available.
  if (!exists("build_J_cluster") &&
      file.exists(file.path("code", "helpers", "sandwich_functions.R"))) {
    source(file.path("code", "helpers", "sandwich_functions.R"))
  }
  t5_meat_checked <- exists("build_J_cluster")
  cellT5 <- grid[grid$J == 20 & grid$nbar_j == 10 &
                 grid$informativeness == "none" &
                 grid$structure == "psu_within_groups", ]
  t5_universe_ok <- TRUE; t5_meat_maxerr <- 0; t5_empty_reps <- 0L
  for (r in 1:8) {
    d5 <- generate_survey_data(cellT5, seed = 77000 + r)
    Gu <- length(d5$psu_strata_map)
    t5_universe_ok <- t5_universe_ok &&
      identical(Gu, 20L * 4L) &&
      max(d5$psu) <= Gu &&
      all(d5$strata == d5$psu_strata_map[d5$psu])
    if (d5$diagnostics$n_psu_empty > 0L) t5_empty_reps <- t5_empty_reps + 1L

    if (t5_meat_checked) {
      p5 <- ncol(d5$X); J5 <- max(d5$group); dd <- p5 + J5
      rr <- d5$weights * (d5$y - 0.5)      # dummy residuals: algebra check only
      Tm <- matrix(0, Gu, dd)
      for (i in seq_along(rr)) {
        g <- d5$psu[i]
        Tm[g, 1:p5] <- Tm[g, 1:p5] + rr[i] * d5$X[i, ]
        Tm[g, p5 + d5$group[i]] <- Tm[g, p5 + d5$group[i]] + rr[i]
      }
      Jman <- matrix(0, dd, dd)
      for (h in unique(d5$psu_strata_map)) {
        idx <- which(d5$psu_strata_map == h); C_h <- length(idx)
        Th  <- sweep(Tm[idx, , drop = FALSE], 2,
                     colMeans(Tm[idx, , drop = FALSE]))
        Jman <- Jman + C_h / (C_h - 1) * crossprod(Th)
      }
      Jman <- (Jman + t(Jman)) / 2
      Jbld <- build_J_cluster(
        X = d5$X, group = d5$group, r = rr, p = p5, J = J5,
        N = length(rr), cluster = d5$psu, strata = d5$strata,
        cluster_strata = d5$psu_strata_map,
        center = TRUE, df_correct = TRUE)
      t5_meat_maxerr <- max(t5_meat_maxerr, max(abs(Jman - Jbld)))
    }
  }

  res <- list(t3_pps_maxerr = t3_pps, t3_poisson_maxerr = t3_poi,
              t3_poisson_smallcell_maxerr = t3_poi_small,
              empty_rate_observed = empty_rate,
              empty_rate_theory   = empty_theory,
              t5_universe_ok = t5_universe_ok,
              t5_empty_reps = t5_empty_reps,
              t5_meat_checked = t5_meat_checked,
              t5_meat_maxerr = t5_meat_maxerr,
              arms = do.call(rbind, out))

  # ---- Verdicts ---------------------------------------------------------------
  a <- res$arms
  informative_rows <- grep("moderate|strong", rownames(a))
  none_rows        <- grep("none", rownames(a))
  res$pass <- list(
    T1_signal_in_informative = all(abs(a[informative_rows, "cor_shock_logw"]) > 0.15),
    T1_no_signal_in_none     = all(abs(a[none_rows, "cor_shock_logw"]) < 0.10),
    T2_response_selection    = all(a[informative_rows, "dw_y1_minus_y0"] < 0),
    T3_exact_probabilities   = (t3_pps < 0.025 && t3_poi < 0.025 &&
                                t3_poi_small < 0.025 &&
                                abs(empty_rate - empty_theory) < 0.015),
    T4_deff_separation       = all(a[grep("strong", rownames(a)), "deff_realized"] >
                                   a[grep("none", rownames(a)), "deff_realized"] + 0.3),
    T5_empty_psu_universe    = (t5_universe_ok && t5_empty_reps > 0L &&
                                (!t5_meat_checked || t5_meat_maxerr < 1e-10))
  )
  res$all_pass <- all(unlist(res$pass))

  if (verbose) {
    cat("== DGP acceptance test ==\n")
    cat(sprintf("T3 pps max |pi_hat - pi|     : %.4f\n", t3_pps))
    cat(sprintf("T3 poisson max |pi_hat - pi| : %.4f\n", t3_poi))
    cat(sprintf("T3 small-cell (N=25, n=2.5)  : %.4f | empty rate %.4f (theory %.4f)\n",
                t3_poi_small, empty_rate, empty_theory))
    cat(sprintf("T5 empty-PSU universe        : universe ok %s | reps with empty PSU %d/8 | meat check %s (maxerr %.1e)\n",
                t5_universe_ok, t5_empty_reps,
                ifelse(t5_meat_checked, "run", "skipped"), t5_meat_maxerr))
    print(round(res$arms, 3))
    cat("Pass flags:\n"); print(unlist(res$pass))
    cat(sprintf("ALL PASS: %s\n", res$all_pass))
  }

  invisible(res)
}
