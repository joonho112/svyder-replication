# Replication package: Design Effect Ratios for Bayesian Survey Models

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R >= 4.5](https://img.shields.io/badge/R-%3E%3D%204.5-1f65b7.svg)](https://www.r-project.org/)

This repository is the replication package for

> **Design Effect Ratios for Bayesian Survey Models: A Diagnostic Framework for Identifying Survey-Sensitive Parameters.**
> JoonHo Lee, Matthew R. Williams, and Terrance D. Savitsky (2026). [arXiv:2603.07791 [stat.ME]](https://arxiv.org/abs/2603.07791).

It contains the code, the Stan model, and the shipped precomputed results needed to reproduce every figure and table in the paper, together with the scripts that regenerate those results from raw data. A companion R package, [**svyder**](https://github.com/joonho112/svyder), implements the diagnostic; this package uses it and reproduces the paper's two studies.

A rendered walkthrough of the package is available as the [replication guide](https://joonho112.github.io/svyder-replication/).

---

## Overview

When a Bayesian hierarchical model is fit to complex survey data through a weighted pseudo-posterior, not every parameter requires a design correction. The **Design Effect Ratio (DER)** is a per-parameter diagnostic: for a declared variance target, it is the ratio of a parameter's design-based sandwich variance to its model-based posterior variance. A DER near 1 means the posterior is already calibrated and should be left alone. A DER above 1 means the posterior understates design uncertainty and should be widened. A DER below 1 means the posterior is already wider than the target, so correcting it would wrongly narrow the interval and remove the shrinkage protection that a hierarchical prior provides.

The DER drives a **Compute–Classify–Correct (CCC)** workflow: compute the DER for every parameter, flag the parameters whose DER exceeds the threshold `c0 = 1.2`, and apply a selective block-Cholesky correction to only the flagged block. Parameters sort into four tiers — I-a within-group fixed effects (exposed), I-b between-group fixed effects (shielded), II random effects (protected by shrinkage), and III hyperparameters (excluded). A *declared variance target* is fixed by two choices: the aggregation unit (design PSU versus model group) and the weight-scaling convention. The same posterior draws can therefore be evaluated against more than one target, and the paper shows that the choice of target materially changes which parameters are flagged.

---

## What is reproduced

The paper reports two studies, both of which this package reproduces exactly.

**Simulation study.** A `2 × 2 × 3 × 2 = 24`-cell factorial design crossing the number of groups `J ∈ {20, 50}`, the expected group size `n̄_j ∈ {10, 50}`, the informativeness of the sampling design `∈ {none, moderate, strong}`, and a structural arm `∈ {sampled groups (design PSU = model group), PSUs within groups (design PSU ≠ model group)}`, with 200 replications per cell, giving **4,800 fitted models**. The data-generating mechanism is a genuine two-stage informative sampling design over finite populations, with inclusion probabilities that depend on the random effects and the outcome and exact inverse-probability weights. The outcome model is a hierarchical logistic regression,

```
y | theta_j ~ Bernoulli( expit( beta_0 + beta_1 * x_within + beta_2 * x_between + theta_j ) ),
```

with true `(beta_0, beta_1, beta_2) = (-0.5, 0.5, 0.3)` and the latent-scale intraclass correlation fixed at 0.15. A strict convergence gate (split-Rhat `<= 1.01` on all parameters and zero divergent transitions) is applied; 99.3% of replications pass, and 33 of the 4,800 are excluded. Selective correction restores near-nominal coverage of the exposed within-group coefficient while leaving the protected parameters untouched, whereas blanket correction collapses random-effect coverage. The within-group DER exceeds every non-target DER by a factor of at least 7.1, and at `c0 = 1.2` the per-replication false-positive rate is 2.1%.

**Empirical application.** The 2019 National Survey of Early Care and Education (NSECE), restricted-use data (ICPSR Study 38445), with `N = 6,785` center-based providers, `J = 51` states (including the District of Columbia), and 415 design PSUs within 30 strata. The outcome is participation in a state Infant-Toddler quality-improvement system, modeled with a survey-weighted hierarchical logistic regression with 54 parameters (3 fixed effects and 51 state random intercepts). The DER is evaluated under two declared targets. Under the model-group target, exactly 1 of the 54 parameters is flagged — the within-state poverty coefficient, with DER 4.59. Under the design-PSU target, 28 of the 54 are flagged: the poverty coefficient (DER 3.55) together with 27 of the 51 state effects. Selective correction widens the poverty coefficient's 90% interval by a factor of 1.88, whereas blanket correction would instead have narrowed 53 of the 54 intervals. After MCMC, the full DER pipeline runs in well under a second.

---

## Two-track replication

Both tracks reproduce the paper's numbers exactly. Choose the one that matches what you need.

| | Track A — full rerun | Track B — figures and tables |
|---|---|---|
| **Starts from** | raw data | shipped precomputed results |
| **Simulation** | regenerates all 4,800 replications on synthetic data (HPC recommended) | reads per-cell summaries |
| **Application** | fits the NSECE model; needs the restricted-use data and Stan | reads the shipped model fit |
| **Requires** | R, CmdStan, `svyder`, and (for the application) ICPSR 38445 | R and CRAN packages only |
| **Approximate time** | many hours for the simulation; a couple of hours for the NSECE fit | minutes on a laptop |

Track B is the default entry point and needs no restricted data and no Stan.

---

## Quick start (Track B)

From the package root, or after opening `svyder-replication.Rproj` in RStudio:

```r
source("code/00_setup.R")                        # one-time Track B setup
source("code/03_results/results_01_theory.R")    # Figure 1
source("code/03_results/results_02_aggregate.R") # Figures 2-4 and Tables 2-5
```

Outputs are written to `output/figures/` (PDF) and `output/tables/` (LaTeX and CSV). These scripts read only the shipped precomputed results; they require no restricted data and no Stan.

---

## Requirements

- **R ≥ 4.5.** Track B needs only CRAN packages used by the reporting scripts and companion book. Track A additionally needs **CmdStan (≥ 2.35)** through `cmdstanr`, the **svyder** package, and restricted NSECE data for the empirical application.
- **Track B setup** — from the package root, run:

  ```r
  source("code/00_setup.R")
  ```

  This is the default setup path and skips `cmdstanr`, CmdStan, and `svyder`.
- **Track A setup** — opt in explicitly before sourcing the setup script:

  ```r
  Sys.setenv(SVYDER_TRACK = "A")
  source("code/00_setup.R")
  ```

  Track A installs/checks `cmdstanr`, CmdStan, and `svyder`.
- **svyder** — for Track A or your own survey models, install with

  ```r
  pak::pak("joonho112/svyder")
  # or: remotes::install_github("joonho112/svyder")
  ```

  Documentation is at <https://joonho112.github.io/svyder/>.

`code/00_setup.R` is track-aware. It defaults to Track B; set `SVYDER_TRACK = "A"` for the full Track A toolchain.

---

## Repository structure

```
svyder-replication/
├── README.md
├── LICENSE                       MIT
├── CITATION.cff
├── svyder-replication.Rproj
├── code/
│   ├── 00_setup.R                track-aware setup (Track B default; Track A with SVYDER_TRACK=A)
│   ├── 01_simulation/            build the 4,800 simulation replications (Track A)
│   │   ├── sim_00_config.R        the 24-cell factorial design, seeds
│   │   ├── sim_01_dgp.R           two-stage informative-sampling data generator
│   │   ├── sim_02_fit.R           fit the weighted HLR in Stan
│   │   ├── sim_03_postprocess.R   sandwich variance, DER, coverage per replication
│   │   ├── sim_04_run_single.R    one replication end to end
│   │   └── sim_05_run_batch.R     parallel driver over all cells and replications
│   ├── 02_application/           the NSECE 2019 analysis (Track A; restricted data)
│   │   ├── app_01_prepare_data.R  build the model-ready arrays from the NSECE extract
│   │   ├── app_02_fit_refit.R     fit the model, dual-target DER  -> nsece_v2_refit.rds
│   │   └── app_03_der_supplement.R weight-convention sensitivity and R_k -> nsece_supp.rds
│   ├── 03_results/               regenerate the paper's figures and tables (Track B)
│   │   ├── results_01_theory.R    Figure 1 (closed-form DER curves)
│   │   └── results_02_aggregate.R Figures 2-4 and Tables 2-5
│   └── helpers/                  self-contained DER / sandwich / plotting utilities
├── stan/
│   └── hlr_weighted.stan          weighted hierarchical logistic regression
├── data/
│   ├── DATA_ACCESS.md             how to obtain the restricted NSECE data
│   └── precomputed/               shipped, non-identifiable results for Track B
│       ├── simulation/            per-cell simulation summary
│       └── application/           NSECE model fit and supplement
├── output/
│   ├── simulation/                the 4,800 simulation replication files
│   ├── figures/                   generated figures (PDF)
│   └── tables/                    generated tables (LaTeX + CSV)
└── book/                          source of the companion guide (rendered to GitHub Pages)
```

---

## Data availability

The simulation study runs entirely on synthetic data generated within the package, so it needs no external inputs.

The empirical application uses restricted-use microdata from the 2019 NSECE, distributed through ICPSR as **Study 38445**. These data cannot be redistributed and are required only for Track A of the application. Instructions for requesting access are in [`data/DATA_ACCESS.md`](data/DATA_ACCESS.md). Track B reproduces every application figure and table from the shipped, non-identifiable precomputed results in `data/precomputed/application/`, so no restricted data are needed to reproduce the paper's reported numbers.

---

## Reproducibility notes

All seeds are set in `code/01_simulation/sim_00_config.R`, so a full Track A rerun reproduces the same replications. The precomputed results shipped in `data/precomputed/` reproduce the paper's figures and tables exactly, and Track B depends on nothing else. The NSECE restricted-use data are required only for Track A of the application; every reported number can otherwise be regenerated from the shipped results.

---

## Citation

If you use this package, please cite both the paper and the software.

**Paper**

```bibtex
@article{lee2026der,
  author        = {Lee, JoonHo and Williams, Matthew R. and Savitsky, Terrance D.},
  title         = {Design Effect Ratios for Bayesian Survey Models: A Diagnostic
                   Framework for Identifying Survey-Sensitive Parameters},
  year          = {2026},
  journal       = {arXiv preprint arXiv:2603.07791},
  eprint        = {2603.07791},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME},
  doi           = {10.48550/arXiv.2603.07791}
}
```

**Software (svyder)**

```bibtex
@manual{svyder,
  author = {Lee, JoonHo and Williams, Matthew R. and Savitsky, Terrance D.},
  title  = {svyder: Design Effect Ratios for Bayesian Survey Models},
  year   = {2026},
  note   = {R package},
  url    = {https://github.com/joonho112/svyder}
}
```

---

## Authors

- **JoonHo Lee**, The University of Alabama — corresponding author and package maintainer ([jlee296@ua.edu](mailto:jlee296@ua.edu), GitHub [@joonho112](https://github.com/joonho112), ORCID [0009-0006-4019-8703](https://orcid.org/0009-0006-4019-8703))
- **Matthew R. Williams**, U.S. Bureau of Labor Statistics
- **Terrance D. Savitsky**, U.S. Bureau of Labor Statistics

> Any opinions expressed in this package are those of the authors and do not constitute policy of the U.S. Bureau of Labor Statistics.

---

## Links

- Paper: [arXiv:2603.07791 [stat.ME]](https://arxiv.org/abs/2603.07791) — <https://doi.org/10.48550/arXiv.2603.07791>
- svyder package: <https://github.com/joonho112/svyder> — documentation at <https://joonho112.github.io/svyder/>
- Replication guide (companion website): <https://joonho112.github.io/svyder-replication/>
- NSECE restricted data: ICPSR Study 38445 (see [`data/DATA_ACCESS.md`](data/DATA_ACCESS.md))

---

## License

Released under the MIT License. Copyright (c) 2026 JoonHo Lee, Matthew R. Williams, and Terrance D. Savitsky. See [`LICENSE`](LICENSE).
