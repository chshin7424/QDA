# 비대각 원소 추정 X (Sigma0_11, Sigma0_22, Sigma1_11, Sigma1_22 각각 독립 추정하는 대각 2x2 모델)
# [수정] precompute.f 도입: model/P/Q/pi1에만 의존하고 thetas에는 무관한 값들을
#        optim 루프 진입 전 1회만 계산 -> eff.f/naive.f/bbse.f/rlls.f 호출 시그니처가
#        (model, P, Q, pi1, thetas) -> (P, pi1, precomp, thetas) 로 변경됨
# [수정] U.f의 sigma_score를 for-loop 없이 벡터화 (invSigma 대칭성 이용)

library(MASS)
library(mvtnorm)

# working classifier ------------------------------------------------------
classifier_P <- function(P) {
  fit.qda <- qda(yp ~ ., data = P)
  return(fit.qda)
}


# score function -------------------------------------------------------------
U.f <- function(y, x, thetas) {
  x <- as.matrix(x)
  
  alpha  <- as.numeric(thetas$alpha)
  mu0    <- as.numeric(thetas$mu0)
  mu1    <- as.numeric(thetas$mu1)
  Sigma0 <- as.matrix(thetas$Sigma0)
  Sigma1 <- as.matrix(thetas$Sigma1)
  
  invSigma0 <- solve(Sigma0)
  invSigma1 <- solve(Sigma1)
  
  f0 <- mvtnorm::dmvnorm(x, mean = mu0, sigma = Sigma0) + 1e-100
  f1 <- mvtnorm::dmvnorm(x, mean = mu1, sigma = Sigma1) + 1e-100
  q1 <- alpha * f1 / (alpha * f1 + (1 - alpha) * f0)
  
  u1 <- (y - q1) / (alpha * (1 - alpha))
  
  x_minus_mu0 <- t(t(x) - mu0)
  u2 <- -(y - q1) * (x_minus_mu0 %*% invSigma0)
  
  x_minus_mu1 <- t(t(x) - mu1)
  u3 <- (y - q1) * (x_minus_mu1 %*% invSigma1)
  
  # [수정] for-loop 제거: invSigma가 항상 대칭행렬이라는 성질을 이용해 벡터화.
  # v_i = invSigma %*% x_centered_i 라고 하면
  # invSigma %*% (xi %*% t(xi)) %*% invSigma 의 (j,j) 대각원소는 v_i[j]^2 와 같다.
  sigma_score <- function(x_centered, invSigma, Sigma) {
    V <- x_centered %*% invSigma   # n x 2, V[i, j] = v_i[j] (invSigma 대칭성 이용)
    
    diag1 <- 0.5 * (V[, 1]^2 - invSigma[1, 1])
    diag2 <- 0.5 * (V[, 2]^2 - invSigma[2, 2])
    
    out <- cbind(diag1 * 2 * Sigma[1, 1],
                 diag2 * 2 * Sigma[2, 2])
    return(out)
  }
  
  u4_base <- sigma_score(x_minus_mu0, invSigma0, Sigma0)
  u4 <- -(y - q1) * u4_base
  
  u5_base <- sigma_score(x_minus_mu1, invSigma1, Sigma1)
  u5 <- (y - q1) * u5_base
  
  result <- cbind(u1, u2, u3, u4, u5)
  return(result)
}


# ============================================================================
# 사전 계산 함수 (precompute)
# ============================================================================
# thetas에 무관한 값들을 optim 루프 바깥에서 1회만 계산하기 위한 함수.
# 반환값(precomp 리스트)을 각 모델 함수에 전달하면 내부 재계산이 발생하지 않는다.

precompute.f <- function(model, P, Q, pi1) {
  nP    <- nrow(P)
  nQ    <- nrow(Q)

  # P와 Q를 합친 전체 데이터프레임 (pyx 계산용)
  X_all <- data.frame(rbind(P[, -1], Q[, -1]))
  A     <- data.frame(y_dummy = NA, X_all)

  # pyx: 모델과 데이터에만 의존 → thetas 무관
  pyx_all <- pyx.f(model, A)   # nP+nQ 행
  pyx_P   <- pyx_all[1:nP, ]
  pyx_Q   <- pyx_all[(nP + 1):(nP + nQ), ]

  # rho (naive/eff 공통): model, P, Q에만 의존 → thetas 무관
  p_hat   <- mean(P$yp == "1")
  q_hat   <- mean(pyx_Q[, 2])
  rho_vec <- c(`0` = (1 - q_hat) / (1 - p_hat),
               `1` =  q_hat      /  p_hat)

  # w: rho와 pyx에만 의존 → thetas 무관
  ratio <- pi1 / (1 - pi1)
  w0    <- (rho_vec[1]^2 + ratio * rho_vec[1]) * pyx_all[, 1]
  w1    <- (rho_vec[2]^2 + ratio * rho_vec[2]) * pyx_all[, 2]
  w_all <- 1 / (w0 + w1 + 1e-16)

  # rho_bbse: model, P, Q에만 의존 → thetas 무관
  C_hat        <- as.matrix(CM.f(model, P))
  pred_Q_prob  <- colMeans(pyx_Q)
  theta_hat_b  <- qr.solve(C_hat + diag(1e-6, nrow(C_hat)), pred_Q_prob)
  theta_hat_b  <- as.numeric(theta_hat_b)
  theta_hat_b[theta_hat_b < 0] <- 0
  theta_hat_b  <- theta_hat_b / sum(theta_hat_b)
  q_hat_b      <- theta_hat_b[2]
  rho_bbse_vec <- c(`0` = (1 - q_hat_b) / (1 - p_hat),
                    `1` =  q_hat_b       /  p_hat)

  # rho_rlls: model, P, Q에만 의존 → thetas 무관
  rho_rlls_vec <- rho.rlls.f(model, P, Q)

  # yp_numeric: 루프마다 반복 변환되던 것을 1회 계산
  yp_numeric <- as.numeric(as.character(P$yp))

  list(
    nP           = nP,
    nQ           = nQ,
    X_all        = X_all,
    pyx_all      = pyx_all,
    pyx_P        = pyx_P,
    pyx_Q        = pyx_Q,
    rho_vec      = rho_vec,
    w_all        = w_all,
    rho_bbse_vec = rho_bbse_vec,
    rho_rlls_vec = rho_rlls_vec,
    yp_numeric   = yp_numeric
  )
}


# ============================================================================
# eff model
# ============================================================================

# theta.f 
theta.f <- function(model, P) {
  alpha_hat <- model$prior[2]
  mu0_hat <- model$means["0", ]
  mu1_hat <- model$means["1", ]
  
  x_data <- P[, -1, drop = FALSE]
  y_data <- P$yp
  
  x_group0 <- x_data[y_data == "0", , drop = FALSE]
  x_group1 <- x_data[y_data == "1", , drop = FALSE]
  
  Sigma0_hat <- cov(x_group0)
  Sigma1_hat <- cov(x_group1)
  
  thetas_list <- list(
    alpha  = alpha_hat,
    mu0    = mu0_hat,
    mu1    = mu1_hat,
    Sigma0 = Sigma0_hat,
    Sigma1 = Sigma1_hat
  )
  
  return(thetas_list)
}

# pyx.f --------------------------------------------------------------------
pyx.f <- function(model, Q) {
  new_data <- Q[, -1, drop = FALSE]
  pred_p <- predict(model, newdata = new_data)$posterior
  p1 <- pred_p[, "1"]
  return(data.frame("0" = 1 - p1, "1" = p1))
}

# working importance weight  -------------------------------------------------
rho.f <- function(model, P, Q) {
  p_hat <- mean(P$yp == "1")
  q_pred_prob <- pyx.f(model, Q)[, 2]
  q_hat <- mean(q_pred_prob)
  return(c(`0` = (1 - q_hat) / (1 - p_hat), `1` = q_hat / p_hat))
}

# conditional expectation term -----------------------------------------------
w.f <- function(precomp) {
  return(precomp$w_all)
}

E.f <- function(precomp, U_y0, U_y1) {
  rho_0 <- precomp$rho_vec[1]
  rho_1 <- precomp$rho_vec[2]
  
  p_y0.x <- precomp$pyx_all[, 1]
  p_y1.x <- precomp$pyx_all[, 2]
  
  Term0 <- (rho_0^2 * U_y0) * p_y0.x
  Term1 <- (rho_1^2 * U_y1) * p_y1.x
  
  E_values <- Term0 + Term1
  return(E_values)
}

a.f <- function(P, pi1, precomp, U_y0, U_y1, U_P) {
  yp_numeric <- precomp$yp_numeric
  nP         <- precomp$nP

  w_all   <- precomp$w_all
  pxy_all <- precomp$pyx_all
  rho_vec <- precomp$rho_vec

  E_all <- E.f(precomp, U_y0, U_y1)

  u_P       <- U_P
  w_P       <- w_all[1:nP]
  pxy_all_P <- pxy_all[1:nP, ]
  E_P       <- E_all[1:nP, ]
  
  M_00 <- mean(w_P[yp_numeric == 0] * rho_vec[1] * pxy_all_P[yp_numeric == 0, 1])
  M_01 <- mean(w_P[yp_numeric == 0] * rho_vec[2] * pxy_all_P[yp_numeric == 0, 2])
  M_10 <- mean(w_P[yp_numeric == 1] * rho_vec[1] * pxy_all_P[yp_numeric == 1, 1])
  M_11 <- mean(w_P[yp_numeric == 1] * rho_vec[2] * pxy_all_P[yp_numeric == 1, 2])
  M <- matrix(c(M_00, M_01, M_10, M_11), nrow = 2, byrow = TRUE)
  
  R_inside <- u_P - (w_P * E_P)
  R_row1 <- colMeans(R_inside[yp_numeric == 0, , drop = FALSE], na.rm = TRUE)
  R_row2 <- colMeans(R_inside[yp_numeric == 1, , drop = FALSE], na.rm = TRUE)
  R <- rbind(R_row1, R_row2)
  
  A_mat <- tryCatch({
    solve(M, R)
  }, error = function(e) {
    qr.solve(M + diag(1e-6, nrow(M)), R)
  })
  
  return(A_mat)
}

E_star.f <- function(P, pi1, precomp, U_y0, U_y1, U_P) {
  rho_0 <- precomp$rho_vec[1]
  rho_1 <- precomp$rho_vec[2]
  
  p_y0.x <- precomp$pyx_all[, 1]
  p_y1.x <- precomp$pyx_all[, 2]
  
  a   <- a.f(P, pi1, precomp, U_y0, U_y1, U_P)
  a_0 <- a[1, ]
  a_1 <- a[2, ]
  
  Term0 <- (rho_0^2 * U_y0) + rho_0 * matrix(rep(a_0, nrow(U_y0)), nrow = nrow(U_y0), byrow = TRUE)
  Term0 <- Term0 * p_y0.x
  
  Term1 <- (rho_1^2 * U_y1) + rho_1 * matrix(rep(a_1, nrow(U_y1)), nrow = nrow(U_y1), byrow = TRUE)
  Term1 <- Term1 * p_y1.x
  
  E_values <- Term0 + Term1
  return(E_values)
}

# working estimating function : influence function ---------------------------
eff.f <- function(P, pi1, precomp, thetas) {
  nP     <- precomp$nP
  nTotal <- precomp$nP + precomp$nQ
  
  U_y0 <- U.f(y = 0, x = precomp$X_all, thetas)
  U_y1 <- U.f(y = 1, x = precomp$X_all, thetas)
  yp_num <- precomp$yp_numeric
  U_P  <- as.matrix(U.f(yp_num, P[, -1], thetas))
  
  w_star_values <- precomp$w_all
  E_star_values <- E_star.f(P, pi1, precomp, U_y0, U_y1, U_P)
  rho_vec       <- precomp$rho_vec
  
  w_P    <- as.vector(w_star_values[1:nP])
  E_P    <- as.matrix(E_star_values[1:nP, ])
  rho_P  <- ifelse(P$yp == "0", rho_vec[1], rho_vec[2])
  
  term_P_in <- U_P - (w_P * E_P)
  phi_P     <- (1 / pi1) * (rho_P * term_P_in)
  
  w_Q   <- as.vector(w_star_values[(nP + 1):nTotal])
  E_Q   <- as.matrix(E_star_values[(nP + 1):nTotal, ])
  phi_Q <- (1 / (1 - pi1)) * (w_Q * E_Q)
  
  phi_all <- rbind(phi_P, phi_Q)
  return(phi_all)
}

# 공분산 파라미터화 ----------------------------------------------------------
# [수정] 2개의 파라미터(v[1], v[2])를 받아 Sigma_11, Sigma_22를 독립적으로 반환
# (비대각 원소는 여전히 0으로 고정, 두 대각 원소만 서로 다르게 추정)
make_sigma <- function(v) {
  L <- matrix(c(exp(v[1]), 0,
                0, exp(v[2])), nrow = 2, byrow = TRUE)
  Sigma <- L %*% t(L)
  return(Sigma)
}

sigma_to_param <- function(Sigma) {
  L <- t(chol(Sigma))
  return(c(log(L[1, 1]), log(L[2, 2])))
}


# ============================================================================
# naive model
# ============================================================================

naive.f <- function(P, pi1, precomp, thetas) {
  rho_vec <- precomp$rho_vec
  rho_P   <- ifelse(P$yp == "0", rho_vec[1], rho_vec[2])
  
  yp  <- precomp$yp_numeric
  U_P <- as.matrix(U.f(yp, P[, -1], thetas))
  
  naive_values <- (1 / pi1) * (rho_P * U_P)
  return(naive_values)
}


# ============================================================================
# BBSE model
# ============================================================================

CM.f <- function(model, P) {
  pred_class <- predict(model, newdata = P[, -1, drop = FALSE])$class
  class_matrix <- table(
    factor(pred_class, levels = c("0","1")),
    factor(P$yp,       levels = c("0","1")))
  C_hat <- prop.table(class_matrix, margin = 2)
  return(C_hat)
}

rho.bbse.f <- function(model, P, Q) {
  C_hat        <- as.matrix(CM.f(model, P))
  pred_Q_prob  <- colMeans(pyx.f(model, Q))
  theta_hat    <- qr.solve(C_hat + diag(1e-6, nrow(C_hat)), pred_Q_prob)
  theta_hat    <- as.numeric(theta_hat)
  theta_hat[theta_hat < 0] <- 0
  theta_hat    <- theta_hat / sum(theta_hat)
  q_hat        <- theta_hat[2]
  p_hat        <- mean(P$yp == "1")
  return(c(`0` = (1 - q_hat) / (1 - p_hat),
           `1` =  q_hat       /  p_hat))
}

bbse.f <- function(P, pi1, precomp, thetas) {
  rho_bbse_vec <- precomp$rho_bbse_vec
  rho_bbse     <- ifelse(P$yp == "0", rho_bbse_vec[1], rho_bbse_vec[2])
  
  yp  <- precomp$yp_numeric
  U_P <- as.matrix(U.f(yp, P[, -1], thetas))
  
  bbse_values <- (1 / pi1) * (rho_bbse * U_P)
  return(bbse_values)
}


# ============================================================================
# RLLS model
# ============================================================================

pred_dist.f <- function(model, P) {
  post  <- predict(model, newdata = P[, -1, drop = FALSE])$posterior
  hat1  <- mean(post[, "1"])
  hat0  <- 1 - hat1
  return(c(hat0, hat1))
}

CM_joint.f <- function(model, P) {
  pred_class <- predict(model, newdata = P[, -1, drop = FALSE])$class
  M <- table(
    factor(pred_class, levels = c("0","1")),
    factor(P$yp,       levels = c("0","1"))
  )
  return(as.matrix(M) / nrow(P))
}

norm.f <- function(x) sqrt(sum(x^2))

rho.rlls.f <- function(model, P, Q, delta = 0.05, alpha_rlls = 0.01) {
  np       <- nrow(P)
  K        <- 2
  CM_joint <- CM_joint.f(model, P)
  mu_p_hat <- colMeans(pyx.f(model, P))
  mu_q_hat <- colMeans(pyx.f(model, Q))
  b        <- mu_q_hat - mu_p_hat
  rho      <- 3 * (2 * log(2 * K / delta) / (3 * np) + sqrt(2 * log(2 * K / delta) / np))
  lambda   <- alpha_rlls * rho
  opt <- optim(
    par    = rep(0, K),
    fn     = function(theta) norm.f(CM_joint %*% theta - b) + lambda * norm.f(theta),
    method = "L-BFGS-B",
    lower  = rep(-1, K)
  )
  theta_hat        <- opt$par
  rho_rlls_vec     <- 1 + theta_hat
  rho_rlls_vec[rho_rlls_vec < 0] <- 0
  names(rho_rlls_vec) <- c("0", "1")
  return(rho_rlls_vec)
}

rlls.f <- function(P, pi1, precomp, thetas) {
  rho_rlls_vec <- precomp$rho_rlls_vec
  rho_rlls     <- ifelse(P$yp == "0", rho_rlls_vec[1], rho_rlls_vec[2])
  
  yp  <- precomp$yp_numeric
  U_P <- as.matrix(U.f(yp, P[, -1, drop = FALSE], thetas))
  
  rlls_values <- (1 / pi1) * (rho_rlls * U_P)
  return(rlls_values)
}
