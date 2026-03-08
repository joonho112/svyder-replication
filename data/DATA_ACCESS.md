# Data Access Instructions

## Overview

The empirical application in this paper uses restricted-use data from the
2019 National Survey of Early Care and Education (NSECE). Due to data use
agreements, the raw microdata cannot be included in this repository.

This replication package provides **two tracks**:

- **Track A (Full replication):** Requires obtaining the restricted-use NSECE
  data and running all analyses from scratch.
- **Track B (Pre-computed results):** Uses pre-computed summary statistics and
  posterior draws stored in `precomputed/` to reproduce all figures and tables
  without access to the restricted data.

## Track A: Full Replication

### Data Source

**NSECE 2019 Restricted-Use Data**

- **Repository:** ICPSR (Inter-university Consortium for Political and Social
  Research)
- **Study Number:** 38445
- **URL:** <https://www.icpsr.umich.edu/web/ICPSR/studies/38445>
- **Title:** National Survey of Early Care and Education (NSECE), 2019:
  [United States] (ICPSR 38445)

### Application Process

1. Visit the ICPSR study page linked above.
2. Click "Get Data" and follow the restricted-use data application process.
3. You will need to submit a data use agreement (DUA) signed by your
   institution's authorized representative.
4. **Expected timeline:** 4-8 weeks for approval.

### Required Variables

The analysis uses the following variables from the NSECE 2019 household
survey (HH file):

| Variable         | Description                                    |
|------------------|------------------------------------------------|
| `FIPS_STATE`     | State FIPS code (grouping variable, J = 51)    |
| `PSU`            | Primary sampling unit identifier               |
| `HHWEIGHT`       | Household survey weight                        |
| `POVER_RATIO`    | Poverty-to-income ratio                        |
| `TIERED_REIM`    | Whether state has tiered reimbursement (0/1)   |
| `ENROLLMENT`     | Binary: child enrolled in ECE program (0/1)    |

### Data Preparation

After obtaining the restricted-use data:

1. Place the raw data files in `data/nsece/` (this directory is gitignored).
2. Run `code/01_data_prep.R` to clean and prepare the analysis dataset.
3. Run `code/02_fit_model.R` to fit the Stan model and compute DER values.
4. Run `code/03_figures/` scripts to generate all figures.

## Track B: Pre-Computed Results

All summary statistics, DER values, and posterior summaries needed to
reproduce the paper's figures and tables are stored in `precomputed/`:

| File                          | Contents                                   |
|-------------------------------|-------------------------------------------|
| `sim_summary_all_j.rds`      | Simulation study aggregated results        |
| `nsece_der_results.rds`      | NSECE application DER values and tiers     |
| `nsece_posterior_summary.rds` | Posterior means and credible intervals     |

To reproduce figures using pre-computed results:

```r
source("code/00_setup.R")
source("code/03_figures/fig_01_simulation.R")
source("code/03_figures/fig_02_application.R")
```

No restricted data access is required for Track B.

## Data Citation

National Center for Education Statistics. (2022). National Survey of Early
Care and Education (NSECE), 2019: [United States]. ICPSR 38445-v1. Ann
Arbor, MI: Inter-university Consortium for Political and Social Research
[distributor]. https://doi.org/10.3886/ICPSR38445.v1
