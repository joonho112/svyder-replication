# Pre-computed Results

This directory contains pre-computed intermediate and final results that enable
**Track B** replication (figures and tables only) without re-running the
computationally intensive simulation study or requiring restricted-use data.

## Directory Structure

```
precomputed/
├── simulation/
│   ├── summary/              # Per-scenario summaries (54 files)
│   │   ├── J020_CV030_ICC005_IN_summary.rds
│   │   ├── J020_CV030_ICC005_NI_summary.rds
│   │   └── ... (54 files total, one per scenario)
│   └── analysis/             # Aggregated analysis objects
│       ├── analysis_all_j.rds          # Master analysis across all J values
│       ├── der_master_j20.rds          # DER values for J=20 scenarios
│       ├── coverage_master_j20.rds     # Coverage results for J=20
│       ├── coverage_by_strategy_j20.rds
│       ├── coverage_deviation_j20.rds
│       ├── der_by_scenario_j20.rds
│       ├── correction_intensity_j20.rds
│       ├── convergence_j20.rds
│       ├── scenario_info_j20.rds
│       └── timing_master_j20.rds
└── application/
    ├── phase3_summary.rds              # Posterior summary statistics
    ├── phase3_sandwich.rds             # Sandwich variance components
    ├── phase3_classification.rds       # DER three-tier classification
    ├── phase3_corrected_draws.rds      # Selectively corrected draws
    └── phase3_theorem_verification.rds # Theorem verification results
```

## Simulation Results

Each scenario summary file (`*_summary.rds`) contains, for all 200
replications within that scenario:

- DER values (exact and Laplace-approximated) for all parameters
- Coverage rates under four correction strategies
- Convergence diagnostics (Rhat, ESS, divergences)
- Posterior summaries (means, SDs, quantiles)
- Sandwich variance components (H_obs, J_cluster, V_sand)

The `analysis_all_j.rds` file aggregates results across all 54 scenarios
and provides the data underlying Tables 1--3 and Figures 2--3 in the paper.

### Naming Convention

Scenario IDs follow the pattern `J{nnn}_CV{nnn}_ICC{nnn}_{NI|IN}`:

| Code | Meaning | Example |
|:-----|:--------|:--------|
| `J020` | 20 clusters | `J020_CV030_ICC005_NI` |
| `CV030` | CV_w = 0.3 | DEFF ≈ 1.09 |
| `CV100` | CV_w = 1.0 | DEFF ≈ 2.00 |
| `CV200` | CV_w = 2.0 | DEFF ≈ 5.00 |
| `ICC005` | ICC = 0.05 | Weak clustering |
| `ICC015` | ICC = 0.15 | Moderate clustering |
| `ICC030` | ICC = 0.30 | Strong clustering |
| `NI` | Non-informative weights | |
| `IN` | Informative weights | |

## Application Results (NSECE 2019)

Application pre-computed results are derived from the restricted-use
NSECE 2019 data (ICPSR Study 38445). These files contain only model
outputs (posterior summaries, variance components, classifications) and
do **not** contain any individual-level survey responses.

**Excluded from this directory** (restricted data or oversized):
- `phase3_fit.rds` — Full CmdStan fit object (387 MB, exceeds GitHub limit)
- `phase3_data.rds` — Individual-level survey data (restricted)
- `phase3_analysis.rds` — Analysis object with identifiable variables

## Usage

To use these pre-computed results with Track B scripts:

```r
source("code/00_setup.R")

# Simulation analysis (Track B)
sim_results <- readRDS(file.path(PATHS$precomputed, "simulation",
                                  "analysis", "analysis_all_j.rds"))

# Application analysis (Track B)
classification <- readRDS(file.path(PATHS$precomputed, "application",
                                     "phase3_classification.rds"))
```

## Provenance

All results were generated using:
- R 4.5.1 with cmdstanr 0.8.1 and CmdStan 2.35.0
- MCMC settings: 4 chains × 1,500 post-warmup draws (6,000 total)
- Random seed: 20260303
- Platform: macOS (Apple Silicon)
- Total computation time: ~72 hours for 10,800 simulation replications
