# Replication Package: Design Effect Ratios for Bayesian Survey Models

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.5-276DC3.svg)](https://www.r-project.org/)
[![Stan](https://img.shields.io/badge/CmdStan-%E2%89%A5%202.35-B2001E.svg)](https://mc-stan.org/users/interfaces/cmdstan)

**Paper:** "Design Effect Ratios for Bayesian Survey Models: A Diagnostic Framework for Identifying Survey-Sensitive Parameters"

**Author:** JoonHo Lee, The University of Alabama ([jlee296@ua.edu](mailto:jlee296@ua.edu))

**Preprint:** Lee, J. (2026). Design Effect Ratios for Bayesian Survey Models: A Diagnostic Framework for Identifying Survey-Sensitive Parameters. arXiv preprint.

**Companion R packages:**
- [`svyder`](https://github.com/joonho112/svyder) -- compute, classify, and correct Design Effect Ratios

---

## Paper Overview

When you fit a Bayesian hierarchical model to complex survey data, not every parameter needs the same treatment. Some parameters -- those identified from within-cluster variation -- are sensitive to the survey design and need variance correction. Others -- those protected by hierarchical shrinkage -- are already well calibrated by the model structure. The trouble is knowing which is which.

This paper introduces the **Design Effect Ratio (DER)**, a per-parameter diagnostic that answers that question. The DER is the ratio of the design-corrected posterior variance to the model-based posterior variance. A DER near 1.0 means the parameter is already fine; a DER well above 1.0 means the sampling design is inflating the posterior variance beyond what the model accounts for, and correction is needed.

The key insight is a closed-form decomposition linking the DER to two quantities practitioners already reason about: the classical Kish design effect and the hierarchical shrinkage factor. This decomposition reveals a conservation law -- hierarchical shrinkage does not eliminate design sensitivity, it redistributes it. Random effects absorb protection from the prior, while fixed effects identified from between-cluster variation inherit the displaced sensitivity.

These theoretical results lead to a practical **compute-classify-correct (CCC)** workflow:

1. **Compute** the DER for every parameter (< 0.03 seconds after MCMC)
2. **Classify** each parameter using a threshold of 1.2 into protected vs. survey-sensitive
3. **Correct** only the flagged parameters via targeted sandwich variance adjustment

The framework is validated across **54 simulation scenarios** (10,800 replications of hierarchical logistic regression) and applied to the **2019 National Survey of Early Care and Education** (6,785 childcare providers, 51 U.S. states), where the DER flags exactly 1 of 54 parameters for correction.

---

## Quick Start (Track B -- Figures and Tables from Pre-Computed Results)

Want to see the results without running any models? Track B reproduces every publication figure and table in just a few minutes on an ordinary laptop, using pre-computed posterior summaries and simulation output.

```r
# 1. Install required R packages
source("code/00_setup.R")

# 2. Generate simulation figures (Figures 2-3)
source("code/03_figures/fig_01_simulation.R")

# 3. Generate application figures (Figure 4)
source("code/03_figures/fig_02_application.R")

# 4. Generate application tables
source("code/02_application/app_04_tables.R")

# 5. Outputs appear in output/tables/ and output/figures/
```

That is it. No restricted data, no Stan compilation, no HPC cluster. Everything you need is already in `data/precomputed/`.

---

## Two-Track Replication Design

This package supports two replication tracks to accommodate different levels of data access and computational resources.

| | Track A: Full Replication | Track B: Partial Replication |
|:--|:--------------------------|:----------------------------|
| **What it does** | Reruns the entire pipeline from raw data through final outputs | Reproduces all figures and tables from pre-computed results |
| **Data required** | NSECE 2019 restricted-use data (application only) | None (pre-computed results included) |
| **Compute required** | HPC cluster recommended (10,800 simulation datasets) | Standard laptop (a few minutes) |
| **Time estimate** | ~12--24 hours for simulation; ~2--4 hours for application | ~5 minutes |

Both tracks produce identical final outputs. Choose Track B if you want to inspect the results and verify the figures; choose Track A if you want to confirm the entire computational pipeline end to end.

---

## Requirements

### Software

| Software | Version | Purpose |
|:---------|:--------|:--------|
| R | >= 4.5 | Statistical computing |
| CmdStan | >= 2.35 | Stan model compilation and sampling (Track A only) |
| LaTeX | Any modern distribution | Manuscript compilation (optional) |

### R Packages

Install all dependencies at once with `source("code/00_setup.R")`, or install manually:

**Core packages:**

| Package | Purpose |
|:--------|:--------|
| `svyder` | DER computation, classification, and selective correction |
| `cmdstanr` | Interface to CmdStan for Bayesian model fitting |
| `posterior` | Posterior draws manipulation and summary |
| `survey` | Design-based survey inference and sandwich variance |

**Data manipulation:**

| Package | Purpose |
|:--------|:--------|
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |
| `purrr` | Functional programming utilities |
| `forcats` | Factor manipulation |
| `Matrix` | Sparse matrix computations |
| `mvtnorm` | Multivariate normal distribution |

**Visualization:**

| Package | Purpose |
|:--------|:--------|
| `ggplot2` | Publication figures |
| `patchwork` | Multi-panel figure composition |
| `scales` | Axis formatting |
| `viridis` | Colorblind-friendly color scales |
| `ggrepel` | Non-overlapping text labels |

**Table generation:**

| Package | Purpose |
|:--------|:--------|
| `xtable` | LaTeX table generation |
| `knitr` | Dynamic report generation |

### Installing the Companion Packages

```r
# Install svyder (DER computation, classification, and correction)
# install.packages("pak")
pak::pak("joonho112/svyder")

# Install csSampling (pseudo-posterior and Cholesky correction)
pak::pak("joonho112/csSampling")

# Or using remotes
# install.packages("remotes")
remotes::install_github("joonho112/svyder")
```

Full package documentation and vignettes are available at [https://joonho112.github.io/svyder/](https://joonho112.github.io/svyder/).

---

## Repository Structure

```
svyder-replication/
|
+-- README.md                            # This file
+-- LICENSE                              # MIT License
+-- CITATION.cff                         # Machine-readable citation metadata
+-- .gitignore                           # Git ignore rules
+-- svyder-replication.Rproj             # RStudio project file
|
+-- code/
|   +-- 00_setup.R                       # Environment check + package installation
|   +-- 01_simulation/                   # Simulation study (Section 4)
|   |   +-- sim_00_config.R              # Factorial design: 54 scenarios
|   |   +-- sim_01_dgp.R                 # Data generating process
|   |   +-- sim_02_fit.R                 # Stan model fitting
|   |   +-- sim_03_postprocess.R         # Sandwich variance and DER computation
|   |   +-- sim_04_run.R                 # Parallel execution wrapper
|   |   +-- sim_05_analyze.R             # Aggregation and summary tables
|   +-- 02_application/                  # NSECE application (Section 5)
|   |   +-- app_01_data_prep.R           # Data preparation (requires restricted data)
|   |   +-- app_02_model_fit.R           # Bayesian model fitting
|   |   +-- app_03_der_analysis.R        # DER computation and classification
|   |   +-- app_04_tables.R              # Application tables
|   +-- 03_figures/                      # Publication figures
|   |   +-- fig_01_simulation.R          # Figures 2-3 (simulation results)
|   |   +-- fig_02_application.R         # Figure 4 (NSECE results)
|   +-- helpers/                         # Utility functions
|       +-- der_functions.R              # Core DER computation
|       +-- sandwich_functions.R         # Sandwich variance matrices
|       +-- plotting_helpers.R           # ggplot2 theme and utilities
|
+-- stan/
|   +-- hlr_weighted.stan                # Weighted hierarchical logistic regression
|
+-- data/
|   +-- DATA_ACCESS.md                   # How to obtain restricted NSECE data
|   +-- precomputed/                     # Pre-computed results for Track B
|       +-- simulation/                  # Aggregated simulation summaries
|       +-- application/                 # NSECE model outputs (non-identifiable)
|
+-- output/
    +-- tables/                          # Generated tables (CSV + TEX)
    +-- figures/                         # Generated figures (PDF + PNG)
```

---

## Track B: Reproduce Tables and Figures (No Restricted Data)

**Time estimate:** ~5 minutes on a standard laptop.

Track B uses pre-computed posterior summaries and simulation results stored in `data/precomputed/`. It reproduces every publication figure and table without requiring the restricted-use NSECE data or refitting any models.

### Step 1. Verify your R environment

```r
source("code/00_setup.R")
```

This script checks your R version, installs any missing packages (including `svyder`), and confirms that all dependencies are available. CmdStan is not needed for Track B.

### Step 2. Generate simulation figures

```r
source("code/03_figures/fig_01_simulation.R")
```

This reads from `data/precomputed/simulation/` and produces Figures 2--3 from the paper (DER separation across 54 scenarios, coverage comparison between selective and blanket correction).

### Step 3. Generate application figures

```r
source("code/03_figures/fig_02_application.R")
```

This reads from `data/precomputed/application/` and produces Figure 4 from the paper (NSECE DER profile with threshold at 1.2, selective correction results).

### Step 4. Generate tables

```r
source("code/02_application/app_04_tables.R")
```

Produces all publication tables in both CSV and LaTeX formats.

### Step 5. Inspect outputs

Browse `output/figures/` for PDF and PNG figures, and `output/tables/` for CSV and TEX tables. Each output file is named to match its corresponding table or figure number in the paper.

---

## Track A: Full Replication (Restricted Data + HPC Required)

**Time estimate:** ~12--24 hours for the simulation study on a multi-core machine; ~2--4 hours for the NSECE application.

Track A reproduces the entire analysis pipeline from raw data through final outputs. The simulation study and NSECE application are independent; the simulation requires only synthetic data, while the NSECE application requires restricted-use data.

### Part 1: Simulation Study

**Step 1.** Configure the simulation:

```r
source("code/01_simulation/sim_00_config.R")
```

Sets up the 54-scenario factorial grid (3 J × 3 CV_w × 3 ICC × 2 Informativeness) with 200 replications per scenario.

**Step 2.** Run the full simulation pipeline:

```r
# Option A: Run everything in parallel (recommended)
source("code/01_simulation/sim_04_run.R")

# Option B: Run step by step
source("code/01_simulation/sim_01_dgp.R")      # Generate synthetic data
source("code/01_simulation/sim_02_fit.R")       # Fit Stan models
source("code/01_simulation/sim_03_postprocess.R") # DER computation + coverage
```

`sim_04_run.R` wraps Steps 2--4 in a parallelized execution pipeline for multi-core and HPC environments.

**Step 3.** Aggregate results:

```r
source("code/01_simulation/sim_05_analyze.R")
```

Computes coverage, bias, RMSE, and DER calibration summaries across all replications and scenarios.

### Part 2: NSECE Application

**Step 4.** Place the NSECE restricted-use data:

```
data/restricted/nsece_2019.rds
```

This directory is git-ignored and must never be committed to version control. See `data/DATA_ACCESS.md` for detailed instructions on obtaining these data through OPRE/ICPSR.

**Step 5.** Prepare the analysis data:

```r
source("code/02_application/app_01_data_prep.R")
```

Produces the Stan-ready data objects: 6,785 center-based childcare providers across 51 U.S. jurisdictions (50 states + DC).

**Step 6.** Fit the Bayesian model:

```r
source("code/02_application/app_02_model_fit.R")
```

Fits the survey-weighted hierarchical logistic regression via CmdStan (4 chains, 2,000 warmup + 2,000 sampling iterations).

**Step 7.** Run the DER diagnostic:

```r
source("code/02_application/app_03_der_analysis.R")
```

This is the heart of the replication. It runs the compute-classify-correct pipeline:

1. Computes the DER for all 54 parameters
2. Classifies each parameter using the threshold of 1.2
3. Flags exactly 1 parameter (the within-state poverty coefficient, DER = 2.643) for correction
4. Applies selective sandwich variance correction to the flagged parameter

The entire DER pipeline completes in under 0.03 seconds after MCMC.

**Step 8.** Generate all outputs:

```r
source("code/03_figures/fig_01_simulation.R")
source("code/03_figures/fig_02_application.R")
source("code/02_application/app_04_tables.R")
```

---

## Simulation Study

The simulation evaluates the DER diagnostic framework across a factorial design of survey conditions using hierarchical logistic regression as the working model.

### Data-Generating Process

The DGP is a two-level hierarchical logistic regression:

```
y_ij | theta_j, beta ~ Bernoulli(expit(beta_0 + beta_1 * x_ij^(W) + beta_2 * x_j^(B) + theta_j))
theta_j ~ Normal(0, sigma_theta^2)
```

where `x_ij^(W)` is a within-cluster covariate and `x_j^(B)` is a between-cluster covariate. True parameter values are (beta_0, beta_1, beta_2) = (-0.5, 0.5, 0.3) with n_j = 50 observations per cluster.

This design provides a clean test of the DER classification:
- **beta_1** (within-cluster): should be Tier I-a (exposed to design effects)
- **beta_2** (between-cluster): should be Tier I-b (shielded by hierarchical structure)
- **theta_j** (random effects): should be Tier II (protected by shrinkage)

### Factorial Grid

| Factor | Levels | Values |
|:-------|:-------|:-------|
| Number of clusters (J) | 3 | 20, 50, 100 |
| Weight CV (CV_w) | 3 | 0.3, 1.0, 2.0 |
| Intraclass correlation (ICC) | 3 | 0.05, 0.15, 0.30 |
| Informativeness | 2 | Non-informative, Informative |

Total: 3 × 3 × 3 × 2 = **54 scenarios** × 200 replications = **10,800 datasets**.

### Key Findings

- Selective correction achieves **87--88% coverage** for survey-sensitive parameters, matching blanket correction
- Blanket correction **collapses coverage to 20--21%** for protected parameters; selective correction preserves near-nominal coverage
- Threshold τ = 1.2 produces **zero false positives** with a separation ratio exceeding **4:1**
- The DER classification is robust across all 54 scenarios

---

## Stan Model

| File | Description |
|:-----|:------------|
| `stan/hlr_weighted.stan` | Hierarchical logistic regression with pseudo-posterior survey weighting |

The Stan model implements pseudo-posterior inference by exponentiating individual log-likelihood contributions by the survey weights, following the approach of Savitsky and Toth (2016). The posterior draws are then post-processed using the DER pipeline for selective sandwich correction.

---

## NSECE Application

The empirical application analyzes the 2019 National Survey of Early Care and Education (NSECE), a nationally representative survey of center-based childcare providers.

| Characteristic | Value |
|:---------------|:------|
| Sample size | 6,785 center-based childcare providers |
| Clusters | 51 (U.S. states + District of Columbia) |
| Outcome | Binary indicator of infant/toddler enrollment |
| Model | Survey-weighted hierarchical logistic regression |
| Parameters | 54 (3 fixed effects + 51 state-level random intercepts) |
| Parameters flagged | 1 (within-state poverty coefficient, DER = 2.643) |
| DER pipeline time | < 0.03 seconds |

The DER identified exactly 1 of 54 parameters as survey-sensitive: the within-state poverty coefficient, which draws its identifying variation from within-cluster (within-state) contrasts and therefore inherits the full design effect. All 53 remaining parameters -- including all state-level random intercepts -- were classified as protected, with DER values close to 1.0.

Had blanket correction been applied, the worst-affected interval would have been narrowed to 4.3% of its original width -- effectively destroying the inference for that parameter. Selective correction preserved these intervals while appropriately widening the single flagged parameter.

---

## The `svyder` Package

The [`svyder`](https://github.com/joonho112/svyder) R package ("**s**ur**v**e**y** **d**esign **e**ffect **r**atios") implements the compute-classify-correct workflow described in Algorithm 1 of the paper.

### Key Functions

| Function | Purpose |
|:---------|:--------|
| `der_compute()` | Compute DER values from posterior draws and survey design |
| `der_classify()` | Classify parameters using threshold τ |
| `der_correct()` | Apply selective Cholesky correction to flagged parameters |
| `der_diagnose()` | All-in-one convenience wrapper for the full CCC pipeline |
| `der_decompose()` | Decompose DER into design effect, shrinkage, and protection factors |
| `der_theorem_check()` | Compare closed-form predictions against numerical DER values |
| `plot()` | Diagnostic visualizations (profile, decomposition, correction plots) |
| `summary()` | Tabular summary of DER values, tiers, and flagging decisions |

### Backend Support

`svyder` accepts posterior draws from any Bayesian backend that produces an S × d draws matrix. Native S3 methods are provided for:

- `brms` (brmsfit objects)
- `cmdstanr` (CmdStanMCMC objects)
- `rstanarm` (stanreg objects)

For other samplers, supply the draws matrix along with the response vector, design matrix, group indicators, and survey design variables.

### Minimal Example

```r
library(svyder)
data(nsece_demo)

## All-in-one DER diagnostic
result <- der_diagnose(
  nsece_demo$draws,
  y = nsece_demo$y, X = nsece_demo$X,
  group = nsece_demo$group,
  weights = nsece_demo$weights,
  psu = nsece_demo$psu,
  family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)

summary(result)
plot(result, type = "profile")

## Or step-by-step for more control
result <- der_compute(
  nsece_demo$draws,
  y = nsece_demo$y, X = nsece_demo$X,
  group = nsece_demo$group,
  weights = nsece_demo$weights,
  psu = nsece_demo$psu,
  family = "binomial",
  sigma_theta = nsece_demo$sigma_theta,
  param_types = nsece_demo$param_types
)
result <- der_classify(result, tau = 1.2)
result <- der_correct(result)
corrected_draws <- as.matrix(result)
```

Full documentation and vignettes are available at [https://joonho112.github.io/svyder/](https://joonho112.github.io/svyder/).

---

## File Descriptions

### Code

| File | Track | Description |
|:-----|:------|:------------|
| `code/00_setup.R` | Both | Checks R version, installs missing packages, verifies CmdStan |
| `code/01_simulation/sim_00_config.R` | A | Defines the 54-scenario factorial design |
| `code/01_simulation/sim_01_dgp.R` | A | Data-generating process for hierarchical logistic regression |
| `code/01_simulation/sim_02_fit.R` | A | Fits weighted HLR model to each simulated dataset |
| `code/01_simulation/sim_03_postprocess.R` | A | Computes sandwich variance, DER, classification, coverage |
| `code/01_simulation/sim_04_run.R` | A | Master launcher orchestrating the full simulation pipeline |
| `code/01_simulation/sim_05_analyze.R` | A/B | Aggregates coverage, bias, RMSE; generates summary tables |
| `code/02_application/app_01_data_prep.R` | A | Processes raw NSECE data into Stan-ready format |
| `code/02_application/app_02_model_fit.R` | A | Fits survey-weighted HLR to NSECE data via CmdStan |
| `code/02_application/app_03_der_analysis.R` | A | Runs the full CCC pipeline on NSECE posterior draws |
| `code/02_application/app_04_tables.R` | Both | Generates publication tables from computed or pre-computed results |
| `code/03_figures/fig_01_simulation.R` | Both | Generates simulation study figures (Figures 2--3) |
| `code/03_figures/fig_02_application.R` | Both | Generates NSECE application figures (Figure 4) |
| `code/helpers/der_functions.R` | Both | DER computation, decomposition, and classification helpers |
| `code/helpers/sandwich_functions.R` | Both | Sandwich variance estimation following Binder (1983) |
| `code/helpers/plotting_helpers.R` | Both | ggplot2 publication theme and formatting utilities |

### Stan

| File | Description |
|:-----|:------------|
| `stan/hlr_weighted.stan` | Survey-weighted hierarchical logistic regression with pseudo-likelihood |

### Data

| File | Description |
|:-----|:------------|
| `data/DATA_ACCESS.md` | Instructions for obtaining NSECE 2019 restricted-use data from OPRE/ICPSR |
| `data/precomputed/simulation/` | Aggregated simulation results across 54 scenarios and 10,800 replications |
| `data/precomputed/application/` | NSECE posterior summaries, DER values, and classification results |

### Output

| Directory | Contents |
|:----------|:---------|
| `output/tables/` | Publication-ready tables in CSV and LaTeX formats |
| `output/figures/` | Publication-ready figures in PDF (vector) and PNG (300 dpi) formats |

---

## Table and Figure Index

### Main Body Figures

| Paper | Description | Script |
|:------|:------------|:-------|
| Figure 1 | Conceptual diagram (DER decomposition) | External (TikZ) |
| Figure 2 | DER separation across 54 scenarios | `code/03_figures/fig_01_simulation.R` |
| Figure 3 | Coverage comparison: selective vs. blanket correction | `code/03_figures/fig_01_simulation.R` |
| Figure 4 | NSECE DER profile and selective correction | `code/03_figures/fig_02_application.R` |

---

## Helper Functions

The `code/helpers/` directory contains utility functions used across the replication pipeline.

| File | Description |
|:-----|:------------|
| `der_functions.R` | Core DER computation: sandwich variance, DER values, threshold classification |
| `sandwich_functions.R` | Sandwich variance matrix construction following Binder (1983) |
| `plotting_helpers.R` | ggplot2 publication theme, tier-based color palette, shared plot utilities |

These helper functions mirror the core logic of the `svyder` package. The package provides the production-quality implementation; the helper files here are self-contained versions for replication transparency.

---

## Computational Notes

### Runtime Estimates

| Task | Hardware | Approximate Time |
|:-----|:---------|:-----------------|
| Track B (figures + tables from pre-computed) | Standard laptop | ~5 minutes |
| Simulation: 1 scenario (200 replications) | 8-core machine | ~15 minutes |
| Simulation: all 54 scenarios | 8-core machine | ~12--18 hours |
| NSECE model fitting (4 chains) | 8-core machine | ~2--4 hours |
| DER pipeline (54 parameters) | Any machine | < 0.03 seconds |

### Parallelization

The simulation wrapper `sim_04_run.R` supports parallelization via the `future` and `furrr` packages. The number of parallel workers is set automatically based on available cores. For HPC environments, the script can be adapted for job array submission (SLURM, PBS, etc.).

### Reproducibility

- All random seeds are set explicitly in configuration files
- Stan model compilation is deterministic given the same CmdStan version
- Pre-computed results in `data/precomputed/` were generated on R 4.5.1 with CmdStan 2.35

---

## Restricted Data Access

The NSECE 2019 restricted-use data are required for Track A (application component only). The simulation study runs entirely on synthetic data and does not require any restricted data access.

To obtain the restricted data:
1. Visit the NSECE page at [ICPSR Study 38445](https://www.icpsr.umich.edu/web/ICPSR/studies/38445)
2. Apply for restricted-use data access through OPRE (Office of Planning, Research, and Evaluation)
3. Follow the data use agreement procedures

Detailed instructions are provided in `data/DATA_ACCESS.md`.

The pre-computed results in `data/precomputed/application/` contain only aggregated, non-identifiable summaries (posterior means, standard deviations, DER values) that cannot be used to reconstruct individual-level data.

---

## Citation

If you use this replication package or the associated methodology, please cite:

```bibtex
@article{Lee2026der,
  author  = {Lee, JoonHo},
  title   = {Design Effect Ratios for {B}ayesian Survey Models: A Diagnostic
             Framework for Identifying Survey-Sensitive Parameters},
  journal = {arXiv preprint},
  year    = {2026}
}
```

If you use the `svyder` R package, please also cite:

```bibtex
@software{Lee2026svyder,
  author  = {Lee, JoonHo},
  title   = {{svyder}: Survey Design Effect Ratios for {B}ayesian Models},
  year    = {2026},
  url     = {https://github.com/joonho112/svyder}
}
```

A machine-readable citation file is available in [`CITATION.cff`](CITATION.cff).

---

## Related Work

This paper builds on and extends the following:

- **Pseudo-posterior inference:** Savitsky and Toth (2016); Williams and Savitsky (2020)
- **Cholesky sandwich correction:** Williams and Savitsky (2021); Williams, McGuire, and Savitsky (2025)
- **Classical design effects:** Kish (1965); Skinner (1989)
- **Calibrated Bayes:** Little (2004, 2012)
- **Parent project:** Lee (2026), "A Bayesian Hierarchical Hurdle Beta-Binomial Model for Survey-Weighted Bounded Counts," arXiv preprint. Replication package: [`hurdlebb-replication`](https://github.com/joonho112/hurdlebb-replication)

---

## License

This replication package is released under the [MIT License](LICENSE).

Copyright (c) 2026 JoonHo Lee.

---

## Contact

JoonHo Lee
Department of Educational Studies in Psychology, Research Methodology, and Counseling
The University of Alabama
[jlee296@ua.edu](mailto:jlee296@ua.edu)
