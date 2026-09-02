library(snowfall)
library(MASS)
library(mvtnorm)
library(dplyr)
library(tidyr)

# 시뮬레이션 설정
n_replicates <- 1000
n_cpus <- 62

# 이상치 기준 설정
ALPHA_MIN <- 0.0001
ALPHA_MAX <- 2
MU_MIN <- -2
MU_MAX <- 8
SIGMA_MIN <- 0.0000001
SIGMA_MAX <- 5

# ============================================================================
# [전역] 데이터 생성 참값(true parameter) 설정
# ============================================================================
p <- 2

mu.yp <- 0.2       # P(Y=1): source 도메인 prior
theta <- 0.9        # Q(Y=1): target 도메인 prior
alpha_true <- theta

mu0_true <- rep(0, p)
mu1_true <- rep(2, p)

Sigma0_true <- matrix(c(0.8, 0,
                        0, 0.8), nrow = 2, byrow = TRUE)
Sigma1_true <- matrix(c(0.5, 0,
                        0, 0.5), nrow = 2, byrow = TRUE)

stopifnot(all(eigen(Sigma0_true, symmetric = TRUE)$values > 0))
stopifnot(all(eigen(Sigma1_true, symmetric = TRUE)$values > 0))

run_simulation_worker <- function(iter) {
  ecount <- 0
  eflag  <- 1
  
  while (eflag > 0) {
    
    seed_val <- (as.integer(iter) + as.integer(ecount) * 1000L) %% .Machine$integer.max
    set.seed(seed_val)
    
    # sample size
    n <- 1000
    n1 <- 400
    pi1 <- n1 / n
    
    # 데이터 생성 --------------------------------------------------------------
    yp <- rbinom(n1, 1, mu.yp)
    xp <- matrix(NA, nrow = n1, ncol = p)
    
    for (i in 1:n1) {
      if (yp[i] == 1) {
        xp[i, ] <- mvrnorm(1, mu1_true, Sigma1_true)
      } else {
        xp[i, ] <- mvrnorm(1, mu0_true, Sigma0_true)
      }
    }
    
    P <- data.frame(yp = factor(yp), xp)
    
    yq_true <- rbinom(n - n1, 1, alpha_true)
    yq_masked <- factor(yq_true * 0)
    
    xq <- matrix(NA, nrow = n - n1, ncol = p)
    for (i in 1:(n - n1)) {
      if (yq_true[i] == 1) {
        xq[i, ] <- mvrnorm(1, mu1_true, Sigma1_true)
      } else {
        xq[i, ] <- mvrnorm(1, mu0_true, Sigma0_true)
      }
    }
    
    Q <- data.frame(yp = yq_masked, xq)
    colnames(Q) <- colnames(P)
    
    # 초기 모델 적합 -----------------------------------------------------------
    fit.qda <- classifier_P(P)
    
    # 결과 저장 객체
    results <- list()
    
    # optim 시작값 ------------------------------------------------------------
    start_p1 <- log(alpha_true / (1 - alpha_true))
    start_p2 <- mu0_true
    start_p3 <- mu1_true
    start_p4 <- sigma_to_param(Sigma0_true)
    start_p5 <- sigma_to_param(Sigma1_true)
    
    start_values_transformed <- c(start_p1, start_p2, start_p3, start_p4, start_p5)

    # =========================================================================
    # thetas에 무관한 공통 값을 optim 루프 바깥에서 1회만 계산
    # =========================================================================
    precomp <- precompute.f(model = fit.qda, P = P, Q = Q, pi1 = pi1)

    
    #=========================================================================
    # Efficient 모델
    #=========================================================================
    results$Efficient <- tryCatch({
      estimating_function_eff <- function(params_transformed) {
        p1 <- params_transformed[1]
        p2 <- params_transformed[2:3]
        p3 <- params_transformed[4:5]
        p4 <- params_transformed[6:7]
        p5 <- params_transformed[8:9]
        
        current_alpha  <- 1 / (1 + exp(-p1))
        current_mu0    <- p2
        current_mu1    <- p3
        current_Sigma0 <- make_sigma(p4)
        current_Sigma1 <- make_sigma(p5)
        
        current_thetas <- list(
          alpha  = current_alpha,
          mu0    = current_mu0,
          mu1    = current_mu1,
          Sigma0 = current_Sigma0,
          Sigma1 = current_Sigma1
        )
        
        phi_values <- eff.f(
          P       = P,
          pi1     = pi1,
          precomp = precomp,
          thetas  = current_thetas
        )
        
        colMeans(phi_values)
      }
      
      sol_eff <- optim(
        start_values_transformed,
        fn = function(param) sum(estimating_function_eff(param)^2),
        method = "BFGS",
        control = list(reltol = 1e-5)
      )
      
      final_p <- sol_eff$par
      
      final_alpha  <- 1 / (1 + exp(-final_p[1]))
      final_mu0    <- final_p[2:3]
      final_mu1    <- final_p[4:5]
      final_Sigma0 <- make_sigma(final_p[6:7])
      final_Sigma1 <- make_sigma(final_p[8:9])
      
      c(
        final_alpha,
        final_mu0[1], final_mu0[2],
        final_mu1[1], final_mu1[2],
        final_Sigma0[1,1], final_Sigma0[2,2],
        final_Sigma1[1,1], final_Sigma1[2,2]
      )
      
    }, error = function(e) {
      rep(NA, 9)
    })
    
    #=========================================================================
    # Naive 모델
    #=========================================================================
    results$Naive <- tryCatch({
      estimating_function_naive <- function(params_transformed) {
        p1 <- params_transformed[1]
        p2 <- params_transformed[2:3]
        p3 <- params_transformed[4:5]
        p4 <- params_transformed[6:7]
        p5 <- params_transformed[8:9]
        
        current_alpha  <- 1 / (1 + exp(-p1))
        current_mu0    <- p2
        current_mu1    <- p3
        current_Sigma0 <- make_sigma(p4)
        current_Sigma1 <- make_sigma(p5)
        
        current_thetas <- list(
          alpha  = current_alpha,
          mu0    = current_mu0,
          mu1    = current_mu1,
          Sigma0 = current_Sigma0,
          Sigma1 = current_Sigma1
        )
        
        phi_values <- naive.f(
          P       = P,
          pi1     = pi1,
          precomp = precomp,
          thetas  = current_thetas
        )
        
        colMeans(phi_values)
      }
      
      sol_naive <- optim(
        start_values_transformed,
        fn = function(param) sum(estimating_function_naive(param)^2),
        method = "BFGS",
        control = list(reltol = 1e-5)
      )
      
      final_p <- sol_naive$par
      
      final_alpha  <- 1 / (1 + exp(-final_p[1]))
      final_mu0    <- final_p[2:3]
      final_mu1    <- final_p[4:5]
      final_Sigma0 <- make_sigma(final_p[6:7])
      final_Sigma1 <- make_sigma(final_p[8:9])
      
      c(
        final_alpha,
        final_mu0[1], final_mu0[2],
        final_mu1[1], final_mu1[2],
        final_Sigma0[1,1], final_Sigma0[2,2],
        final_Sigma1[1,1], final_Sigma1[2,2]
      )
      
    }, error = function(e) {
      rep(NA, 9)
    })
    
    #=========================================================================
    # BBSE 모델
    #=========================================================================
    results$BBSE <- tryCatch({
      estimating_function_bbse <- function(params_transformed) {
        p1 <- params_transformed[1]
        p2 <- params_transformed[2:3]
        p3 <- params_transformed[4:5]
        p4 <- params_transformed[6:7]
        p5 <- params_transformed[8:9]
        
        current_alpha  <- 1 / (1 + exp(-p1))
        current_mu0    <- p2
        current_mu1    <- p3
        current_Sigma0 <- make_sigma(p4)
        current_Sigma1 <- make_sigma(p5)
        
        current_thetas <- list(
          alpha  = current_alpha,
          mu0    = current_mu0,
          mu1    = current_mu1,
          Sigma0 = current_Sigma0,
          Sigma1 = current_Sigma1
        )
        
        phi_values <- bbse.f(
          P       = P,
          pi1     = pi1,
          precomp = precomp,
          thetas  = current_thetas
        )
        
        colMeans(phi_values)
      }
      
      sol_bbse <- optim(
        start_values_transformed,
        fn = function(param) sum(estimating_function_bbse(param)^2),
        method = "BFGS",
        control = list(reltol = 1e-5)
      )
      
      final_p <- sol_bbse$par
      
      final_alpha  <- 1 / (1 + exp(-final_p[1]))
      final_mu0    <- final_p[2:3]
      final_mu1    <- final_p[4:5]
      final_Sigma0 <- make_sigma(final_p[6:7])
      final_Sigma1 <- make_sigma(final_p[8:9])
      
      c(
        final_alpha,
        final_mu0[1], final_mu0[2],
        final_mu1[1], final_mu1[2],
        final_Sigma0[1,1], final_Sigma0[2,2],
        final_Sigma1[1,1], final_Sigma1[2,2]
      )
      
    }, error = function(e) {
      rep(NA, 9)
    })
    
    #=========================================================================
    # RLLS 모델
    #=========================================================================
    results$RLLS <- tryCatch({
      estimating_function_rlls <- function(params_transformed) {
        p1 <- params_transformed[1]
        p2 <- params_transformed[2:3]
        p3 <- params_transformed[4:5]
        p4 <- params_transformed[6:7]
        p5 <- params_transformed[8:9]
        
        current_alpha  <- 1 / (1 + exp(-p1))
        current_mu0    <- p2
        current_mu1    <- p3
        current_Sigma0 <- make_sigma(p4)
        current_Sigma1 <- make_sigma(p5)
        
        current_thetas <- list(
          alpha  = current_alpha,
          mu0    = current_mu0,
          mu1    = current_mu1,
          Sigma0 = current_Sigma0,
          Sigma1 = current_Sigma1
        )
        
        phi_values <- rlls.f(
          P       = P,
          pi1     = pi1,
          precomp = precomp,
          thetas  = current_thetas
        )
        
        colMeans(phi_values)
      }
      
      sol_rlls <- optim(
        start_values_transformed,
        fn = function(param) sum(estimating_function_rlls(param)^2),
        method = "BFGS",
        control = list(reltol = 1e-5)
      )
      
      final_p <- sol_rlls$par
      
      final_alpha  <- 1 / (1 + exp(-final_p[1]))
      final_mu0    <- final_p[2:3]
      final_mu1    <- final_p[4:5]
      final_Sigma0 <- make_sigma(final_p[6:7])
      final_Sigma1 <- make_sigma(final_p[8:9])
      
      c(
        final_alpha,
        final_mu0[1], final_mu0[2],
        final_mu1[1], final_mu1[2],
        final_Sigma0[1,1], final_Sigma0[2,2],
        final_Sigma1[1,1], final_Sigma1[2,2]
      )
      
    }, error = function(e) {
      rep(NA, 9)
    })
    
    #=========================================================================
    # 결과 정리 및 이상치 판정
    #=========================================================================
    result_vec <- c(results$Efficient, results$Naive, results$BBSE, results$RLLS)
    is_outlier <- FALSE
    
    if (sum(is.na(result_vec)) > 0) {
      is_outlier <- TRUE
    } else {
      n_models <- 4
      
      for (i in 1:n_models) {
        model_result <- result_vec[(1 + (i - 1) * 9):(9 + (i - 1) * 9)]
        
        if (model_result[1] < ALPHA_MIN || model_result[1] > ALPHA_MAX) {
          is_outlier <- TRUE
          break
        }
        
        if (any(model_result[2:5] < MU_MIN) || any(model_result[2:5] > MU_MAX)) {
          is_outlier <- TRUE
          break
        }
        
        if (any(model_result[c(6:9)] < SIGMA_MIN) ||
            any(model_result[c(6:9)] > SIGMA_MAX)) {
          is_outlier <- TRUE
          break
        }
      }
    }
    
    if (is_outlier) {
      ecount <- ecount + 1
    } else {
      eflag <- 0
    }
  }
  
  return(c(result_vec, ecount))
}


# 병렬 실행 ===============================================

sfInit(parallel = TRUE, cpus = n_cpus)

sfLibrary(MASS)
sfLibrary(mvtnorm)
sfLibrary(dplyr)
sfLibrary(tidyr)

sfExport(
  "run_simulation_worker",
  "classifier_P",
  "make_sigma",
  "sigma_to_param",
  "theta.f",
  "pyx.f",
  "w.f",
  "E.f",
  "a.f",
  "E_star.f",
  "U.f",
  "rho.f",
  "naive.f",
  "eff.f",
  "CM.f",
  "rho.bbse.f",
  "bbse.f",
  "pred_dist.f",
  "CM_joint.f",
  "norm.f",
  "rho.rlls.f",
  "rlls.f",
  "precompute.f",
  "ALPHA_MIN", "ALPHA_MAX",
  "MU_MIN", "MU_MAX",
  "SIGMA_MIN", "SIGMA_MAX",
  "p", "mu.yp", "theta", "alpha_true",
  "mu0_true", "mu1_true",
  "Sigma0_true", "Sigma1_true"
)

results_list <- sfLapply(1:n_replicates, run_simulation_worker)
sfStop()

# 데이터 프레임으로 변환 ===================================

results_df <- do.call(rbind, results_list) %>% as.data.frame()

param_names <- c(
  "Alpha",
  "Mu0_1", "Mu0_2",
  "Mu1_1", "Mu1_2",
  "Sigma0_11", "Sigma0_22",
  "Sigma1_11", "Sigma1_22"
)

model_names <- c("Efficient", "Naive", "BBSE", "RLLS")

col_names <- c()
for (model in model_names) {
  col_names <- c(col_names, paste0(model, "_", param_names))
}
col_names <- c(col_names, "ecount")

colnames(results_df) <- col_names

results_df_analysis <- results_df[, 1:(length(model_names) * length(param_names))]

# 결과 요약 ===============================================

true_values <- data.frame(
  Parameter = param_names,
  True_Value = c(
    theta,
    mu0_true[1], mu0_true[2],
    mu1_true[1], mu1_true[2],
    Sigma0_true[1, 1], Sigma0_true[2, 2],
    Sigma1_true[1, 1], Sigma1_true[2, 2]
  )
)

summary_list <- list()

for (model in model_names) {
  model_cols <- paste0(model, "_", param_names)
  
  model_data <- results_df_analysis[, model_cols]
  colnames(model_data) <- param_names
  
  model_long <- model_data %>%
    mutate(Iter = row_number()) %>%
    pivot_longer(cols = -Iter, names_to = "Parameter", values_to = "Estimate")
  
  model_summary <- model_long %>%
    left_join(true_values, by = "Parameter") %>%
    group_by(Parameter) %>%
    summarise(
      Model = model,
      True_Value = mean(True_Value),
      Est_Mean = mean(Estimate, na.rm = TRUE),
      Bias = mean(Estimate, na.rm = TRUE) - mean(True_Value),
      SE = sd(Estimate, na.rm = TRUE),
      MSE = mean((Estimate - True_Value)^2, na.rm = TRUE),
      RMSE = sqrt(mean((Estimate - True_Value)^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ round(., 4)))
  
  summary_list[[model]] <- model_summary
}

all_summary <- bind_rows(summary_list)
print(all_summary, n = Inf)

ecount_summary <- sum(results_df$ecount)
cat("\n=== ecount ===\n")
ecount_summary

#결과 엑셀로 저장
write.csv(all_summary, "all_summary.csv", row.names = FALSE)
getwd()


# ============================================================================
# [4] Classification Accuracy Evaluation (QDA 버전) - [수정] 속도 최적화
# ============================================================================

# 대각공분산 전용 log-density: mvtnorm::dmvnorm(..., sigma = diag(sigma_diag))와 수치적으로 동일
log_dnorm_diag <- function(x, mu, sigma_diag) {
  p_dim <- length(mu)
  quad <- rowSums(sweep(x, 2, mu, "-")^2 /
                    matrix(sigma_diag, nrow(x), p_dim, byrow = TRUE))
  logdet <- sum(log(sigma_diag))
  -0.5 * (p_dim * log(2 * pi) + logdet + quad)
}

eval_performance <- function(xq_data, yq_actual, est_alpha,
                              est_mu1, est_mu0,
                              est_sigma1_diag, est_sigma0_diag, p) {
  # [수정] mvtnorm::dmvnorm(..., sigma = diag(...)) 대신 대각전용 계산 사용
  log_f1_hat <- log_dnorm_diag(xq_data, est_mu1, est_sigma1_diag)
  log_f0_hat <- log_dnorm_diag(xq_data, est_mu0, est_sigma0_diag)

  log_prob1_hat <- log(est_alpha) + log_f1_hat
  log_prob0_hat <- log(1 - est_alpha) + log_f0_hat

  y_pred <- ifelse(log_prob1_hat > log_prob0_hat, 1, 0)

  TP <- as.numeric(sum(yq_actual == 1 & y_pred == 1))
  TN <- as.numeric(sum(yq_actual == 0 & y_pred == 0))
  FP <- as.numeric(sum(yq_actual == 0 & y_pred == 1))
  FN <- as.numeric(sum(yq_actual == 1 & y_pred == 0))

  Pos <- TP + FN
  Neg <- TN + FP
  accuracy <- (TP + TN) / (TP + TN + FP + FN)

  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(denom == 0, 0, (TP * TN - FP * FN) / denom)

  return(c(Pos = Pos, Neg = Neg, TP = TP, TN = TN, FP = FP, FN = FN,
           Accuracy = accuracy, MCC = mcc))
}

# Q_test 데이터셋 재생성 (전역 참값을 그대로 재사용 → run_simulation_worker와 항상 동일) --
set.seed(999)
n_test <- 10000
p_dim  <- p

yq_test <- rbinom(n_test, 1, theta)
xq_test <- matrix(NA, nrow = n_test, ncol = p_dim)
for (i in 1:n_test) {
  if (yq_test[i] == 1) {
    xq_test[i, ] <- mvrnorm(1, mu1_true, Sigma1_true)
  } else {
    xq_test[i, ] <- mvrnorm(1, mu0_true, Sigma0_true)
  }
}

# base line -------------------------------------------------------------

oracle_thetas <- list(
  alpha = theta, mu1 = mu1_true, mu0 = mu0_true,
  sigma1_diag = c(Sigma1_true[1, 1], Sigma1_true[2, 2]),
  sigma0_diag = c(Sigma0_true[1, 1], Sigma0_true[2, 2])
)

# Source: P 도메인의 (보정 안 된) prior인 mu.yp를 그대로 사용,
# 조건부분포(mu0, mu1, Sigma0, Sigma1)는 covariate shift가 없으므로 참값과 동일
source_thetas <- list(
  alpha = mu.yp, mu1 = mu1_true, mu0 = mu0_true,
  sigma1_diag = c(Sigma1_true[1, 1], Sigma1_true[2, 2]),
  sigma0_diag = c(Sigma0_true[1, 1], Sigma0_true[2, 2])
)

# evaluation --------------------------------------------------------------
performance_list <- list()

res_oracle <- eval_performance(
  xq_test, yq_test,
  oracle_thetas$alpha, oracle_thetas$mu1, oracle_thetas$mu0,
  oracle_thetas$sigma1_diag, oracle_thetas$sigma0_diag, p_dim
)
performance_list[["Oracle"]] <- data.frame(Model = "Oracle", t(round(res_oracle, 4)))

res_source <- eval_performance(
  xq_test, yq_test,
  source_thetas$alpha, source_thetas$mu1, source_thetas$mu0,
  source_thetas$sigma1_diag, source_thetas$sigma0_diag, p_dim
)
performance_list[["Source"]] <- data.frame(Model = "Source", t(round(res_source, 4)))

# 시뮬레이션으로 추정된 4가지 모델 성능 평가 --------------------------------
# [수정] for + list 누적 대신 vapply로 반복 계산 (재할당 오버헤드 제거)
for (model in model_names) {
  col_alpha  <- paste0(model, "_Alpha")
  col_mu0    <- paste0(model, "_Mu0_", 1:p_dim)
  col_mu1    <- paste0(model, "_Mu1_", 1:p_dim)
  col_sigma0 <- paste0(model, "_Sigma0_", c("11", "22"))
  col_sigma1 <- paste0(model, "_Sigma1_", c("11", "22"))

  valid_idx <- which(!is.na(results_df_analysis[[col_alpha]]))

  metrics_mat <- vapply(valid_idx, function(i) {
    eval_performance(
      xq_data         = xq_test,
      yq_actual       = yq_test,
      est_alpha       = results_df_analysis[i, col_alpha],
      est_mu1         = as.numeric(results_df_analysis[i, col_mu1]),
      est_mu0         = as.numeric(results_df_analysis[i, col_mu0]),
      est_sigma1_diag = as.numeric(results_df_analysis[i, col_sigma1]),
      est_sigma0_diag = as.numeric(results_df_analysis[i, col_sigma0]),
      p               = p_dim
    )
  }, FUN.VALUE = numeric(8))

  # metrics_mat: 8(지표) x length(valid_idx) 행렬 -> 행별 평균이 지표별 평균
  avg_metrics <- rowMeans(metrics_mat, na.rm = TRUE)
  performance_list[[model]] <- data.frame(Model = model, t(round(avg_metrics, 4)))
}

performance_df <- do.call(rbind, performance_list)
rownames(performance_df) <- NULL

cat("\n=== [4] 분류 예측 성능 (시뮬레이션", n_replicates, "회 평균) ===\n")
performance_df

# 분류 성능(정확도, MCC 등) 결과를 csv로 저장
write.csv(performance_df, "performance_summary.csv", row.names = FALSE)
getwd()
``
