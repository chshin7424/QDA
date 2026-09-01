<div align="center">

# QDA Label Shift Estimation

**Label Shift 하에서의 준모수 효율적 QDA 파라미터 추정**

[![R](https://img.shields.io/badge/R-%3E%3D4.0-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-research-orange.svg)](#)

Source 도메인과 라벨 비율(P(Y))이 다른 Target 도메인에서 QDA 분류기 파라미터를 보정 추정하고, 4가지 방법론(Efficient / Naive / BBSE / RLLS)의 성능을 시뮬레이션으로 비교합니다.

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
- [License](#license)

---

## Overview

Label shift(사전확률 이동)는 P(Y)는 Source와 Target 도메인 간 다르지만 P(X|Y)는 동일하다고 가정하는 분포 이동(distribution shift) 상황입니다. Source에서 학습된 QDA 분류기를 라벨이 없는 Target 도메인에 그대로 적용하면 사전확률 차이로 인해 성능이 저하됩니다.

이 프로젝트는 다음을 구현합니다.

- Source 데이터(라벨 有)와 Target 데이터(라벨 無)를 함께 사용하는 4가지 파라미터 보정 추정 방법
- `optim()` 기반 M-estimation 프레임워크로 각 방법의 추정방정식(estimating equation) 최소화
- Monte Carlo 시뮬레이션(기본 1,000회 반복, 병렬 처리)을 통한 Bias / SE / MSE / RMSE 비교
- Target 데이터에 대한 분류 성능(Accuracy, MCC) 비교 (Oracle / Source baseline 포함)

## Methods

| Model | Description | Key Function |
|---|---|---|
| **Efficient** | 준모수 효율적 영향함수(semiparametric efficient influence function) 기반 추정 | `eff.f`, `E_star.f`, `a.f` |
| **Naive** | Importance weight(ρ)만 이용한 단순 보정 | `naive.f`, `rho.f` |
| **BBSE** | Black Box Shift Estimation — confusion matrix 기반 라벨 비율 추정 | `bbse.f`, `CM.f` |
| **RLLS** | Regularized Learning under Label Shift — 정규화 최적화 기반 라벨 비율 추정 | `rlls.f`, `rho.rlls.f` |

## Repository Structure

```
.
├── 010_QDA_func_최종.R      # Core functions: score function, precompute, estimating functions
├── 010_QDA_simul_최종.r     # Simulation: data generation, parallel estimation, performance eval
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
source("010_QDA_func_최종.R")

# 2. Run simulation (parallelized via snowfall)
source("010_QDA_simul_최종.r")
```

시뮬레이션 반복 횟수와 병렬 코어 수는 `010_QDA_simul_최종.r` 상단에서 조정할 수 있습니다.

```r
n_replicates <- 1000   # Monte Carlo 반복 횟수
n_cpus <- 62            # 병렬 처리 코어 수
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

이상치(발산, 파라미터 범위 이탈)가 발생하면 시드를 변경해 재시도하며, 판정 기준은 `ALPHA_MIN/MAX`, `MU_MIN/MAX`, `SIGMA_MIN/MAX`로 정의됩니다.

## Core Functions

<details>
<summary>펼쳐서 보기</summary>

| Function | Description |
|---|---|
| `classifier_P()` | Source 데이터(P)로 QDA 분류기 학습 |
| `U.f()` | 파라미터(alpha, mu0, mu1, Sigma0, Sigma1)에 대한 score function 계산 |
| `precompute.f()` | 파라미터(thetas)와 무관한 값(pyx, rho, w, rho_bbse, rho_rlls 등)을 최적화 루프 진입 전 1회 계산 |
| `theta.f()` | QDA 적합 결과로부터 초기 파라미터 추출 |
| `pyx.f()` | 모델 기반 P(Y=1\|X) 사후확률 예측 |
| `make_sigma()` / `sigma_to_param()` | Cholesky 분해 기반 공분산 ↔ 최적화 파라미터(log-scale) 변환 (대각 공분산 모델) |
| `CM.f()` / `CM_joint.f()` | Confusion matrix 계산 (BBSE / RLLS용) |

</details>

## Output

| File | Description |
|---|---|
| `all_summary.csv` | 모델 × 파라미터별 True Value, Est_Mean, Bias, SE, MSE, RMSE |
| `performance_summary.csv` | Oracle / Source / 4개 모델의 Target 분류 성능 (TP, TN, FP, FN, Accuracy, MCC) |

## Results Interpretation

- `all_summary.csv`의 **Bias**, **RMSE**가 작을수록 추정이 정확합니다.
- `performance_summary.csv`에서 **Oracle**은 참값을 알 때의 성능 상한, **Source**는 보정 없이 사용할 때의 baseline입니다. 추정 모델의 Accuracy/MCC가 Oracle에 가깝고 Source보다 높을수록 label shift 보정이 효과적인 것입니다.

## References

- Lipton, Z., Wang, Y. X., & Smola, A. (2018). *Detecting and Correcting for Label Shift with Black Box Predictors*. ICML.
- Azizzadenesheli, K., Liu, A., Yang, F., & Anandkumar, A. (2019). *Regularized Learning for Domain Adaptation under Label Shifts*. ICLR.

## License

이 프로젝트는 [MIT License](LICENSE)를 따릅니다.
