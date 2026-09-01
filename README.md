<div align="center">

# SSDALS : Semi-supervised Discriminant Analysis under Label Shift

Demonstrates the superiority of the proposed semiparametric efficient QDA estimator over existing label-shift correction methods (BBSE, RLLS, Naive) via Monte Carlo simulation under label shift $P(Y)$.
</div>

---

## Table of Contents

- [Overview](#overview)
- [Methods](#methods)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Simulation Settings](#simulation-settings)
- [Core Functions](#core-functions)
- [Output](#output)
- [Results Interpretation](#results-interpretation)
- [References](#references)
- [Research Status](#research-status)
- [Acknowledgement](#acknowledgement)

---

## Overview

Label shift is a distribution-shift setting where P(Y) differs between the source and target domains while P(X | Y) is assumed to remain the same. Applying a QDA classifier trained on the source domain directly to an unlabeled target domain leads to degraded performance because of this prior-probability mismatch.

This project implements:

- Four parameter-correction methods that jointly use labeled source data and unlabeled target data
- An `optim()`-based M-estimation (estimating-equation / Z-estimator) framework that drives each method's estimating equation toward zero
- Monte Carlo simulation (1,000 replicates by default, run in parallel) comparing Bias / SE / MSE / RMSE
- Classification performance (Accuracy, MCC) on target data, including Oracle and Source baselines

## Methods

| Model | Description | Key Function |
|---|---|---|
| **Efficient** | Estimation based on the semiparametric efficient influence function | `eff.f`, `E_star.f`, `a.f` |
| **Naive** | Simple correction using only an importance weight (ρ) | `naive.f`, `rho.f` |
| **BBSE** | Black Box Shift Estimation — label ratio estimated from the confusion matrix | `bbse.f`, `CM.f` |
| **RLLS** | Regularized Learning under Label Shift — label ratio estimated via regularized optimization | `rlls.f`, `rho.rlls.f` |

## Repository Structure

```
.
├── QDA_function.R           # Core functions: score function, precompute, estimating functions
├── QDA_simulation.R         # Simulation: data generation, parallel estimation, performance eval
├── all_summary.csv          # (output) Parameter estimation summary
├── performance_summary.csv  # (output) Classification performance summary
└── README.md
```

## Requirements

- R >= 4.0
- Packages: `MASS`, `mvtnorm`, `snowfall`, `dplyr`, `tidyr`

## Installation

```r
install.packages(c("MASS", "mvtnorm", "snowfall", "dplyr", "tidyr"))
```

```bash
git clone https://github.com/<username>/<repo-name>.git
cd <repo-name>
```

## Usage

```r
# 1. Load core functions
source("QDA_function.R")

# 2. Run simulation (parallelized via snowfall)
source("QDA_simulation.R")
```

The number of Monte Carlo replicates and parallel cores can be adjusted at the top of `QDA_simulation.R`:

```r
n_replicates <- 1000   # number of Monte Carlo replicates
n_cpus <- 62            # number of parallel cores
```

## Simulation Settings

| Parameter | Value |
|---|---|
| Total sample size (`n`) | 1,000 |
| Source sample size (`n1`) | 400 |
| P(Y=1) in Source (`mu.yp`) | 0.2 |
| P(Y=1) in Target (`theta`) | 0.9 |
| μ0 | (0, 0) |
| μ1 | (2, 2) |
| Σ0 | diag(0.8, 0.8) |
| Σ1 | diag(0.5, 0.5) |

If an outlier occurs (divergence or a parameter falling outside the allowed range), the run is retried with a new seed. The outlier criteria are defined by `ALPHA_MIN/MAX`, `MU_MIN/MAX`, and `SIGMA_MIN/MAX`.

## Core Functions

<details>
<summary>Click to expand</summary>

| Function | Description |
|---|---|
| `classifier_P()` | Fits a QDA classifier on the source data (P) |
| `U.f()` | Computes the score function with respect to the parameters (alpha, mu0, mu1, Sigma0, Sigma1) |
| `precompute.f()` | Computes, once before the optimization loop, all quantities (pyx, rho, w, rho_bbse, rho_rlls, etc.) that do not depend on the parameters (thetas) |
| `theta.f()` | Extracts initial parameters from the fitted QDA model |
| `pyx.f()` | Predicts the posterior probability P(Y=1\|X) from the model |
| `make_sigma()` / `sigma_to_param()` | Converts between a covariance matrix and its (log-scale) optimization parameters via Cholesky decomposition (diagonal covariance model) |
| `CM.f()` / `CM_joint.f()` | Computes the confusion matrix (used by BBSE / RLLS) |

</details>

## Output

| File | Description |
|---|---|
| `all_summary.csv` | True Value, Est_Mean, Bias, SE, MSE, RMSE for each model × parameter |
| `performance_summary.csv` | Target classification performance (TP, TN, FP, FN, Accuracy, MCC) for Oracle / Source / the 4 models |

> The examples below were produced with a reduced run (`n_replicates = 30`); actual values will vary with the number of replicates and the random seed. For numbers intended for a paper or report, run the full simulation (`n_replicates = 1000`) yourself.

<details>
<summary><code>all_summary.csv</code> example (Efficient model, 9 parameters)</summary>

| Parameter | True_Value | Est_Mean | Bias | SE | MSE | RMSE |
|---|---|---|---|---|---|---|
| Alpha | 0.9 | 0.9017 | 0.0017 | 0.0036 | 0.0000 | 0.0039 |
| Mu0_1 | 0.0 | -0.0764 | -0.0764 | 0.0878 | 0.0133 | 0.1152 |
| Mu0_2 | 0.0 | -0.0732 | -0.0732 | 0.0838 | 0.0121 | 0.1102 |
| Mu1_1 | 2.0 | 1.9843 | -0.0157 | 0.0846 | 0.0072 | 0.0847 |
| Mu1_2 | 2.0 | 1.9840 | -0.0160 | 0.0992 | 0.0098 | 0.0988 |
| Sigma0_11 | 0.8 | 0.8306 | 0.0306 | 0.2616 | 0.0671 | 0.2590 |
| Sigma0_22 | 0.8 | 0.8175 | 0.0175 | 0.2483 | 0.0599 | 0.2448 |
| Sigma1_11 | 0.5 | 0.4613 | -0.0387 | 0.2246 | 0.0503 | 0.2242 |
| Sigma1_22 | 0.5 | 0.4717 | -0.0283 | 0.1692 | 0.0285 | 0.1688 |

The full CSV repeats this table for the `Naive`, `BBSE`, and `RLLS` models in the same format.

</details>

**`performance_summary.csv` example** (n_replicates = 30, n_test = 10,000)

| Model | Pos | Neg | TP | TN | FP | FN | Accuracy | MCC |
|---|---|---|---|---|---|---|---|---|
| Oracle | 9024 | 976 | 8964.0 | 866.0 | 110.0 | 60.0 | 0.9830 | 0.9016 |
| Source | 9024 | 976 | 8251.0 | 961.0 | 15.0 | 773.0 | 0.9212 | 0.7047 |
| Efficient | 9024 | 976 | 8943.0 | 871.1 | 104.9 | 81.0 | 0.9814 | 0.8939 |
| Naive | 9024 | 976 | 8855.1 | 905.0 | 71.0 | 168.9 | 0.9760 | 0.8728 |
| BBSE | 9024 | 976 | 8885.6 | 887.1 | 88.9 | 138.4 | 0.9773 | 0.8765 |
| RLLS | 9024 | 976 | 8885.7 | 887.8 | 88.2 | 138.3 | 0.9773 | 0.8768 |

Between Oracle (using the true parameters) and Source (no correction), all four estimated models improve Accuracy and MCC over Source and move closer to Oracle, with Efficient coming closest to Oracle.

## Results Interpretation

- In `all_summary.csv`, smaller **Bias** and **RMSE** indicate more accurate estimation.
- In `performance_summary.csv`, **Oracle** represents the performance ceiling when the true parameters are known, and **Source** represents the baseline with no correction. The closer an estimated model's Accuracy/MCC is to Oracle — and the more it exceeds Source — the more effective the label-shift correction.

## References

- Lipton, Z., Wang, Y. X., & Smola, A. (2018). *Detecting and Correcting for Label Shift with Black Box Predictors*. ICML.
- Azizzadenesheli, K., Liu, A., Yang, F., & Anandkumar, A. (2019). *Regularized Learning for Domain Adaptation under Label Shifts*. ICLR.

## Research Status

This repository contains **ongoing and unpublished research**.

The code is provided to document the implementation and experimental development of the project. Full theoretical derivations, detailed algorithms, simulation settings, and complete numerical results are intentionally omitted at this stage.

## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/
