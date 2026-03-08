// ============================================================================
// hlr_weighted.stan
// Weighted Hierarchical Logistic Regression (NCP Parameterization)
// ============================================================================
//
// Paper    : Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models:
//            A Diagnostic Framework for Identifying Survey-Sensitive Parameters.
//            arXiv preprint.
// Author   : JoonHo Lee (jlee296@ua.edu)
// License  : MIT
//
// Description:
//   Hierarchical logistic regression with survey weights, using a
//   non-centered parameterization (NCP) for the group-level random effects.
//   Implements pseudo-maximum likelihood by scaling the log-likelihood
//   contributions by normalized survey weights.
//
//   This model corresponds to Section 2 of the paper. The key elements:
//
//   (1) Model specification (Paper Eq. 1-3):
//       y_ij | p_ij ~ Bernoulli(p_ij)
//       logit(p_ij) = X_i * beta + theta_{group[i]}
//       theta_j = sigma_theta * eta_raw_j    (NCP)
//       eta_raw_j ~ N(0, 1)
//
//   (2) Weighted pseudo-likelihood (Paper Section 2.1):
//       The standard log-likelihood is replaced by:
//         target += w_i * bernoulli_logit_lpmf(y_i | linpred_i)
//       where w_i are normalized survey weights (sum(w) = N).
//       This implements the pseudo-maximum likelihood approach of
//       Savitsky & Toth (2016) for survey-weighted Bayesian models.
//
//   (3) Score residuals (Paper Section 2.3, Algorithm 1):
//       The generated quantities block computes score residuals:
//         s_i = w_i * (y_i - q_i)
//       These are the building blocks for the sandwich variance estimator
//       and subsequent DER computation:
//         J_cluster = sum_j (sum_{i in j} s_i * phi_i^T)
//                         * (sum_{i in j} s_i * phi_i^T)^T
//
//   Design matrix X is expected to have p columns, typically:
//     Column 1: Intercept (1's)
//     Column 2: Within-group covariate (group-mean centered)
//     Column 3: Between-group covariate
//
// Priors:
//   beta ~ N(0, 5^2)               [weakly informative]
//   sigma_theta ~ half-N(0, 2^2)   [weakly informative, via <lower=0>]
//   eta_raw ~ N(0, 1)              [NCP standard normal]
//
// ============================================================================

data {
  int<lower=1> N;                        // total number of observations
  int<lower=1> J;                        // number of groups (clusters)
  int<lower=1> p;                        // number of fixed-effect predictors
  array[N] int<lower=0, upper=1> y;      // binary outcome
  matrix[N, p] X;                        // design matrix
  array[N] int<lower=1, upper=J> group;  // group membership indicator
  vector<lower=0>[N] w;                  // normalized survey weights
}

parameters {
  vector[p] beta;                        // fixed effects (length p)
  real<lower=0> sigma_theta;             // random-effect SD
  vector[J] eta_raw;                     // NCP: standardized random effects
}

transformed parameters {
  // Non-centered parameterization: theta_j = sigma_theta * eta_raw_j
  // This improves sampling efficiency when groups have few observations
  // or when sigma_theta is small (Paper Section 2.1, footnote 3).
  vector[J] theta = sigma_theta * eta_raw;
}

model {
  // --- Priors (Paper Section 2.1) ---
  beta ~ normal(0, 5);                   // weakly informative
  sigma_theta ~ normal(0, 2);            // half-normal via <lower=0>
  eta_raw ~ std_normal();                // NCP prior

  // --- Weighted pseudo-likelihood (Paper Section 2.1, Eq. 4) ---
  // Each observation's log-likelihood contribution is scaled by its
  // survey weight w_i. This implements:
  //   ell_w(phi) = sum_i w_i * log Bernoulli(y_i | logit^{-1}(eta_i))
  // where eta_i = X[i,] * beta + theta[group[i]].
  {
    vector[N] linpred;
    for (i in 1:N) {
      linpred[i] = X[i] * beta + theta[group[i]];
    }
    for (i in 1:N) {
      target += w[i] * bernoulli_logit_lpmf(y[i] | linpred[i]);
    }
  }
}

generated quantities {
  // Score residuals for sandwich variance estimation (Paper Section 2.3)
  //
  // For each observation i:
  //   score_residual[i] = w_i * (y_i - q_i)
  // where q_i = inv_logit(X[i,] * beta + theta[group[i]])
  //
  // These are aggregated post-hoc to form the clustered meat matrix:
  //   J_cluster = sum_j s_j s_j^T
  // where s_j = sum_{i in cluster j} score_residual[i] * phi_i
  // and phi_i is the d-vector of partial derivatives d(eta_i)/d(phi).
  vector[N] score_residual;
  {
    real linpred_i;
    real q_i;
    for (i in 1:N) {
      linpred_i = X[i] * beta + theta[group[i]];
      q_i = inv_logit(linpred_i);
      score_residual[i] = w[i] * (y[i] - q_i);
    }
  }
}
