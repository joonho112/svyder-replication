# Data Access

This note accompanies the replication package for

> Lee, J., Williams, M. R., and Savitsky, T. D. (2026).
> *Design Effect Ratios for Bayesian Survey Models.*

It explains which data are needed, how to obtain the one restricted-use
dataset the package relies on, and — most importantly — the exact structure
of the data objects that the empirical application scripts expect, so that a
reader who has obtained the raw survey files can rebuild those objects and
reproduce the application from scratch.

---

## 1. Who needs this, and who does not

The package is organized into two tracks:

- **Track B (reporting).** Regenerate every figure and table in the paper from
  the precomputed results shipped in this repository. Track B fits no models
  and reads no restricted data. If your goal is to reproduce the paper's
  figures and tables, **you can ignore the rest of this document** and run
  `code/03_results/`.

- **Track A (from-scratch application).** Refit the empirical model and
  recompute its Design Effect Ratios from the underlying survey microdata.
  This is the only part of the package that requires restricted-use data.

The **simulation study** (`code/01_simulation/`) uses only synthetic data
generated internally by the package. It needs nothing described here and no
special access of any kind.

In short, the restricted **National Survey of Early Care and Education (NSECE)
2019** data are required **only for Track A of the application**
(`code/02_application/`). Everything else in the package runs without them.

---

## 2. Obtaining the NSECE 2019 restricted-use data

The application uses restricted-use files from the **2019 National Survey of
Early Care and Education (NSECE)**. The NSECE was conducted by **NORC at the
University of Chicago** under the **Office of Planning, Research, and
Evaluation (OPRE)**, Administration for Children and Families, U.S. Department
of Health and Human Services.

The restricted-use files are distributed through the Inter-university
Consortium for Political and Social Research (ICPSR):

- **Study:** *National Survey of Early Care and Education (NSECE), 2019
  \[United States\]: Restricted-Use Files*
- **ICPSR study number:** **38445**
- **Landing page:** <https://www.icpsr.umich.edu/web/ICPSR/studies/38445>

### Access procedure

The restricted-use files are **not** available for direct download and cannot
be redistributed. Access is granted through ICPSR's standard restricted-data
process. At a general level, this involves:

1. Registering for an ICPSR account and locating study **38445**.
2. Submitting a **restricted-data application** describing the research use,
   the data-security environment, and the personnel who will access the data.
3. Executing a **Restricted Data Use Agreement** (and any accompanying
   institutional or IRB documentation ICPSR requires).
4. Receiving the files through ICPSR's approved secure delivery method once
   the application is approved.

Requirements and delivery mechanisms are set by ICPSR and OPRE and may change
over time; the ICPSR study page above is the authoritative source for the
current terms. We describe the process only in general terms and make no
representation about approval or timelines.

### Suggested data citation

> National Opinion Research Center (NORC) at the University of Chicago, and
> Office of Planning, Research, and Evaluation, Administration for Children and
> Families, U.S. Department of Health and Human Services. *National Survey of
> Early Care and Education (NSECE), 2019 \[United States\]: Restricted-Use
> Files.* Inter-university Consortium for Political and Social Research
> \[distributor\], ICPSR 38445. <https://www.icpsr.umich.edu/web/ICPSR/studies/38445>

Please cite the data in accordance with the citation and acknowledgment terms
specified by ICPSR and OPRE for study 38445.

---

## 3. Where to place the data

The Track A scripts locate the restricted derivatives through the environment
variable **`NSECE_RESTRICTED_DIR`**. If that variable is unset, they fall back
to the git-ignored directory **`data/nsece/`** inside this package. You may
either:

- point the variable at wherever you keep the derivatives, e.g.

  ```r
  Sys.setenv(NSECE_RESTRICTED_DIR = "/path/to/your/nsece/derivatives")
  ```

- or place the derivative files directly in `data/nsece/`.

> **Do not commit any restricted data.** The directories `data/nsece/` and
> `data/restricted/`, along with common microdata file extensions
> (`.csv`, `.dta`, `.sav`, `.sas7bdat`), are listed in `.gitignore` precisely
> so that respondent-level files cannot be added to version control by
> accident. Keep the restricted files inside these ignored locations (or
> outside the repository entirely) and never add them to a commit, a fork, or
> any public mirror. The terms of the Restricted Data Use Agreement govern how
> these files may be stored and shared.

---

## 4. Data contract: the objects the scripts expect

This is the part that lets you rebuild the inputs from your own NSECE extract.
The Track A application scripts read two RDS objects from
`NSECE_RESTRICTED_DIR`. Their required structure is documented below. If you
construct objects with these names, types, and dimensions, the scripts will
run against your data without modification.

The analysis sample is **N = 6,785** provider observations nested in **J = 51**
groups (the 50 states plus the District of Columbia). The outcome is
participation in the state Infant–Toddler quality-improvement system. The
sample comes from a **stratified, multistage cluster design: 415 primary
sampling units (PSUs) within 30 strata.**

### 4.1 `phase3_data.rds` — model arrays and sensitivity weights

A **list** carrying the arrays consumed by the pseudo-likelihood model:

| Element | Type | Description |
|---|---|---|
| `N` | integer, scalar | Number of observations. **6785**. |
| `J` | integer, scalar | Number of groups (states). **51**. |
| `p` | integer, scalar | Number of fixed-effect columns. **3**. |
| `y` | integer, length `N` | Binary outcome (participation in the state Infant–Toddler quality-improvement system), coded `0/1`. |
| `X` | numeric matrix, `N × 3` | Fixed-effects design matrix; columns `[intercept, poverty_cwc, tiered_reim]` (see below). |
| `group` | integer, length `N` | Group (state) index in `1..J`. |
| `w` | numeric, length `N` | Within-state unit-mean weights used by the weight-convention sensitivity analysis. |

Important convention note: `phase3_data$w` is **not** the headline v2 fit
weight. The headline fit in `app_02_fit_refit.R` constructs
`w_v2 <- phase3_analysis$weight / mean(phase3_analysis$weight)` from the raw
weight column. `phase3_data$w` is the within-state alternative used by
`app_03_der_supplement.R`.

The three columns of `X` are:

1. **`intercept`** — a column of ones.
2. **`poverty_cwc`** — the community poverty rate, **group-mean-centered within
   each state** (centering within context). After centering it carries only
   within-state variation, so it is a **purely within-state** covariate; its
   state means are zero to numerical precision.
3. **`tiered_reim`** — a **between-state (state-level)** indicator for a tiered
   reimbursement policy. It is constant within each state.

### 4.2 `phase3_analysis.rds` — enriched analysis data frame

A **data frame with `N` rows** (one per observation, aligned row-for-row with
the arrays in `phase3_data.rds`). The scripts require at least the following
columns:

| Column | Type | Description |
|---|---|---|
| `state_idx` | integer | Group (state) index. **Must equal `group`** in `phase3_data.rds`. |
| `weight` | numeric | The **raw NSECE design weight**, before any normalization. |
| `psu_idx` | integer | Design PSU index, in `1..415`. |
| `stratum_idx` | integer | Design stratum index, in `1..30`. |
| `poverty_cwc` | numeric | Within-state centered poverty rate (as in `X`, column 2). |
| `tiered_reim` | numeric | Between-state tiered-reimbursement indicator (as in `X`, column 3). |

The scripts enforce the alignment explicitly:

```r
stopifnot(length(d$y) == nrow(a), all(d$group == a$state_idx))
```

where `d` is `phase3_data.rds` and `a` is `phase3_analysis.rds`.

### 4.3 The two declared weight conventions

The paper reports results under two weight conventions, and **both are derived
by the scripts from the single raw column `weight`** in `phase3_analysis.rds`.
You supply the raw weight; the scripts form the conventions.

- **Declared convention — global unit-mean.** Divide each raw weight by the
  overall mean of the raw weights, so that the weights average to one across
  the full sample:

  ```r
  w_v2 <- phase3_analysis$weight / mean(phase3_analysis$weight)
  ```

  This is the primary convention used for the headline model fit in
  `app_02_fit_refit.R`. Its Kish design effect on the raw weights is reported
  as **DEFF ≈ 3.76**.

- **Alternative convention — within-state scaling.** Rescale the weights so
  that, within each state, they sum to that state's sample size `n_j` (a
  "sum-to-n" normalization applied group by group):

  ```r
  w_within <- weight * n_j[state] / tapply(weight, state, sum)[state]
  ```

  This convention is stored as `phase3_data$w` by `app_01_prepare_data.R` and
  is used for the weight-convention sensitivity analysis in
  `app_03_der_supplement.R`.

Because the weights enter the pseudo-likelihood itself, switching conventions
yields a **different pseudo-posterior**, not a reweighting of the same draws.
That is why the sensitivity analysis is a genuine re-evaluation rather than a
post-hoc rescaling.

### 4.4 Supplementary derivatives (optional)

`app_03_der_supplement.R` additionally reads two small derivative files if you
wish to reproduce the supplementary table exactly:

- `phase3_corrected_draws.rds` — a list whose element `draws_naive` holds the
  stored posterior draws of the **within-state–convention** fit.
- `phase3_sandwich.rds` — a list whose element `sigma_theta_hat` holds the
  plug-in group-scale estimate.

These carry only posterior summaries and scalars, never respondent-level data.
If you refit the model yourself you can regenerate them from your own run.

---

## 5. What is, and is not, shipped in this package

This package ships **non-identifiable precomputed results only.** It contains
**no respondent-level microdata of any kind.**

Shipped (under `data/precomputed/application/`):

- `nsece_v2_refit.rds` — posterior draws and Design Effect Ratio results of the
  fitted model under the declared convention (both design-PSU and model-group
  targets).
- `nsece_supp.rds` — the weight-convention sensitivity results and the
  structural protection factors.

**Never shipped, and never to be committed:**

- `phase3_data.rds`, `phase3_analysis.rds`, and any raw NSECE extract.
- Any file containing respondent-level records.

Track B reads the two shipped `.rds` files above and needs nothing further.

---

## 6. Running Track A once you have the data

With the restricted derivatives in place, point the scripts at them and run the
two Track A steps from the package root:

```r
# 1. Tell the scripts where your restricted derivatives live
Sys.setenv(NSECE_RESTRICTED_DIR = "/path/to/your/nsece/derivatives")

# 2. Refit the model and compute the dual-target DERs -> nsece_v2_refit.rds
source("code/02_application/app_02_fit_refit.R")

# 3. Weight-convention sensitivity and structural factors R_k -> nsece_supp.rds
source("code/02_application/app_03_der_supplement.R")
```

The first step compiles the Stan model in `stan/hlr_weighted.stan` and refits
under the declared global unit-mean convention, then applies the convergence
gate and the dual-target DER analysis. The second step computes the
supplementary quantities. Both write their outputs to
`data/precomputed/application/`, overwriting the shipped copies with your
freshly produced versions.

Afterwards, the Track B reporting script uses the files you just produced:

```r
source("code/03_results/results_02_aggregate.R")
```

so the figures and tables you regenerate reflect your own refit rather than the
shipped results.

---

## 7. Questions

For questions about the replication package, contact the corresponding author,
JoonHo Lee (<jlee296@ua.edu>). Questions about eligibility for, or access to,
the NSECE 2019 restricted-use files should be directed to ICPSR (study 38445)
and OPRE, which administer the data and set the terms of access.
