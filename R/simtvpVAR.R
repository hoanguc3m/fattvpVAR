#' Simulate data from TVP VAR model with SV
#'
#' This function simulate data from a TVP VAR model with SV.
#' \deqn{A_t y_t = B_t x_t + SQRT(w_t) H_t eps_t}
#' @param dist The variable specifies the VAR error distribution. It should be one of
#' c("Gaussian","Student").
#' @param K The number of variables in VAR model.
#' @param p The number of lags in VAR model.
#' @param t_max The number of observations
#' @param b0 The coefficients of B matrix.
#' @param a0 The coefficients of A matrix.
#' @param h The initial log diag(Sigma_t) matrix. sigma_h is set at diag(0.01, K)
#' @param nu The degree of freedom.
#' @param seednum The default seed.
#' @param burn_in The discarded observations.
#' @return A list of simulated data with its specification.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.tvpVAR.SV(dist="Gaussian")
#' y <- datagen$y
#' }
sim.tvpVAR.SV <- function(dist, K = 3, p = 2, t_max = 1000,
                       b0 = 0.5, a0 = 0.5, h = -3, sigma_ab = NULL, sigma_h = NULL,
                       y0 = matrix(0, ncol = K, nrow = p),
                       nu = 6, is_tv = NULL,
                       seednum = 0, burn_in = 0){
  if (!(dist %in% c("Gaussian","Student") ))
    stop("dist is not implemented.")

  if (dist == "Gaussian") datagen <- sim.tvpVAR.Gaussian.SV(K, p, t_max, b0, a0, h, sigma_ab, sigma_h, y0, is_tv, seednum, burn_in)
  if (dist == "Student") datagen <- sim.tvpVAR.Student.SV(K, p, t_max, b0, a0, h, sigma_ab, sigma_h, y0, nu, is_tv, seednum, burn_in)
  return(datagen)
}

#' @export
sim.tvpVAR.Gaussian.SV <- function(K = 3, p = 2, t_max = 1000,
                                b0 = 0.5, a0 = 0.5, h = -3, sigma_ab = NULL, sigma_h = NULL,
                                y0 = matrix(0, ncol = K, nrow = p),
                                is_tv = NULL, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- B0_mat(b0, K, p)

  # Sample matrix corr A0
  A0 <- A0_mat(a0, K)
  tA0 <- t(A0)
  a0_vec <- tA0[upper.tri(tA0)]

  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # No tail
  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(1e-2, 1e-2, length.out = K)
  } else {
    Vh <-  rep(sigma_h,K)
  }

  if (is.null(sigma_ab)){
    V_ab <- (2e-3)
  } else {
    V_ab <- sigma_ab
  }

  if (is.null(is_tv)){
    is_tv <- rep(1,K)
  }

  # Identify the tvp element
  k_alp <- K*(K-1)/2 # dimension of the impact matrix
  k_beta <- K^2*p + K # number of VAR coefficients
  k_beta_div_K <- K*p + 1
  n_tv <- sum(is_tv)# # of time-varying equations
  idx_b_tv <- (kronecker(is_tv,matrix(1,nrow = K*p+1,ncol = 1))==1)   # index for time-varying betas
  idx_a_tv <- matrix(FALSE, nrow = k_alp, ncol = 1)            # construct index for time-varying alphas
  count <- 0
  for (j in 2:K){
    if (is_tv[j] == 1){
      idx_a_tv[ (count+1):(count+j-1)] <- TRUE
    }
    count <- count + j-1
  }


  B_t <- reprow( as.numeric(t(B0)), t_max) + reprow(idx_b_tv, t_max) * apply(matrix(rnorm(t_max * k_beta) * V_ab, nrow = t_max), MARGIN = 2, FUN = cumsum)
  A_t <- reprow( as.numeric(t(a0_vec)), t_max) + reprow(idx_a_tv, t_max) * apply(matrix(rnorm(t_max * k_alp) * V_ab, nrow = t_max), MARGIN = 2, FUN = cumsum)
  H_t <- reprow( as.numeric(h), t_max) + apply(matrix(rnorm(t_max * K), nrow = t_max) * reprow(Vh, t_max), MARGIN = 2, FUN = cumsum)

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  eps <- matrix(rnorm(t_max*K), ncol = K)



  for (i in c(1:t_max)){
    inv_A0 <- solve(A0_mat(A_t[i,], K))
    Sigma_t <-  inv_A0 %*% diag(exp(0.5*H_t[i,]), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    xt <- rbind(1, matrixcalc::vec( t(ystar[(p+i-1):i,])))
    Bt_mat <- B0_mat(B_t[i,], K, p)
    y_mean[i,] <- inv_A0 %*% (Bt_mat %*% xt)
    ysim <-  y_mean[i,] + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       Vh = Vh, B_t = B_t, A_t = A_t, H_t = H_t,
       dist = "Gaussian", SV = TRUE)
}

#' @export
sim.tvpVAR.Student.SV <- function(K = 3, p = 2, t_max = 1000,
                                  b0 = 0.5, a0 = 0.5, h = -3, sigma_ab = NULL, sigma_h = NULL,
                                  y0 = matrix(0, ncol = K, nrow = p),
                                  nu = 6, is_tv = NULL, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- B0_mat(b0, K, p)

  # Sample matrix corr A0
  A0 <- A0_mat(a0, K)
  tA0 <- t(A0)
  a0_vec <- tA0[upper.tri(tA0)]

  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # Tail of student
  w_t <- matrix(rinvgamma(t_max*K, shape = nu/2, rate = nu/2), ncol = K)
  w_sqrt_t <- sqrt(w_t)

  # Volatility volatility
  if (is.null(sigma_h)){
    Vh <- seq(1e-2, 1e-2, length.out = K)
  } else {
    Vh <-  rep(sigma_h,K)
  }

  if (is.null(sigma_ab)){
    V_ab <- (2e-3)
  } else {
    V_ab <- sigma_ab
  }

  if (is.null(is_tv)){
    is_tv <- rep(1,K)
  }

  # Identify the tvp element
  k_alp <- K*(K-1)/2 # dimension of the impact matrix
  k_beta <- K^2*p + K # number of VAR coefficients
  k_beta_div_K <- K*p + 1
  n_tv <- sum(is_tv)# # of time-varying equations
  idx_b_tv <- (kronecker(is_tv,matrix(1,nrow = K*p+1,ncol = 1))==1)   # index for time-varying betas
  idx_a_tv <- matrix(FALSE, nrow = k_alp, ncol = 1)            # construct index for time-varying alphas
  count <- 0
  for (j in 2:K){
    if (is_tv[j] == 1){
      idx_a_tv[ (count+1):(count+j-1)] <- TRUE
    }
    count <- count + j-1
  }


  B_t <- reprow( as.numeric(t(B0)), t_max) + reprow(idx_b_tv, t_max) * apply(matrix(rnorm(t_max * k_beta) * V_ab, nrow = t_max), MARGIN = 2, FUN = cumsum)
  A_t <- reprow( as.numeric(t(a0_vec)), t_max) + reprow(idx_a_tv, t_max) * apply(matrix(rnorm(t_max * k_alp) * V_ab, nrow = t_max), MARGIN = 2, FUN = cumsum)
  H_t <- reprow( as.numeric(h), t_max) + apply(matrix(rnorm(t_max * K), nrow = t_max) * reprow(Vh, t_max), MARGIN = 2, FUN = cumsum)

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)



  for (i in c(1:t_max)){
    inv_A0 <- solve(A0_mat(A_t[i,], K))
    Sigma_t <-  inv_A0 %*% diag(w_sqrt_t[i,] * exp(0.5*H_t[i,]), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    xt <- rbind(1, matrixcalc::vec( t(ystar[(p+i-1):i,])))
    Bt_mat <- B0_mat(B_t[i,], K, p)
    y_mean[i,] <- inv_A0 %*% (Bt_mat %*% xt)
    ysim <-  y_mean[i,] + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max)],
       Vh = Vh, B_t = B_t, A_t = A_t, H_t = H_t,
       dist = "Student", SV = TRUE)
}


#' Simulate data from TVP VAR model without SV
#'
#' This function simulate data from a TVP VAR model without SV.
#' \deqn{A_t y_t = B_t x_t + SQRT(w_t) H eps_t}
#' @param dist The variable specifies the VAR error distribution. It should be one of
#' c("Gaussian","Student").
#' @param K The number of variables in VAR model.
#' @param p The number of lags in VAR model.
#' @param t_max The number of observations
#' @param b0 The coefficients of B matrix.
#' @param a0 The coefficients of A matrix.
#' @param h The initial log diag(Sigma_t) matrix. sigma_h is set at diag(0.01, K)
#' @param nu The degree of freedom.
#' @param seednum The default seed.
#' @param burn_in The discarded observations.
#' @return A list of simulated data with its specification.
#' @export
#' @examples
#' \dontrun{
#' datagen <- sim.tvpVAR.SV(dist="Gaussian")
#' y <- datagen$y
#' }
sim.tvpVAR.nonSV <- function(dist, K = 3, p = 2, t_max = 1000,
                          b0 = 0.5, a0 = 0.5, h = -3, sigma_ab = NULL,
                          y0 = matrix(0, ncol = K, nrow = p),
                          nu = 6, is_tv = NULL,
                          seednum = 0, burn_in = 0){
  if (!(dist %in% c("Gaussian","Student") ))
    stop("dist is not implemented.")

  if (dist == "Gaussian") datagen <- sim.tvpVAR.Gaussian.nonSV(K, p, t_max, b0, a0, h, sigma_ab, y0, is_tv, seednum, burn_in)
  if (dist == "Student") datagen <- sim.tvpVAR.Student.nonSV(K, p, t_max, b0, a0, h, sigma_ab, y0, nu, is_tv, seednum, burn_in)
  return(datagen)
}

#' @export
sim.tvpVAR.Gaussian.nonSV <- function(K = 3, p = 2, t_max = 1000,
                                   b0 = 0.5, a0 = 0.5, h = -3, sigma_ab = NULL,
                                   y0 = matrix(0, ncol = K, nrow = p),
                                   is_tv = NULL, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- B0_mat(b0, K, p)

  # Sample matrix corr A0
  A0 <- A0_mat(a0, K)
  tA0 <- t(A0)
  a0_vec <- tA0[upper.tri(tA0)]

  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # No tail

  if (is.null(sigma_ab)){
    V_ab <- (2e-3)
  }

  if (is.null(is_tv)){
    is_tv <- rep(1,K)
  }

  # Identify the tvp element
  k_alp <- K*(K-1)/2 # dimension of the impact matrix
  k_beta <- K^2*p + K # number of VAR coefficients
  k_beta_div_K <- K*p + 1
  n_tv <- sum(is_tv)# # of time-varying equations
  idx_b_tv <- (kronecker(is_tv,matrix(1,nrow = K*p+1,ncol = 1))==1)   # index for time-varying betas
  idx_a_tv <- matrix(FALSE, nrow = k_alp, ncol = 1)            # construct index for time-varying alphas
  count <- 0
  for (j in 2:K){
    if (is_tv[j] == 1){
      idx_a_tv[ (count+1):(count+j-1)] <- TRUE
    }
    count <- count + j-1
  }


  B_t <- reprow( as.numeric(t(B0)), t_max) + reprow(idx_b_tv, t_max) * apply(matrix(rnorm(t_max * k_beta) * V_ab, nrow = t_max), MARGIN = 2, FUN = cumsum)
  A_t <- reprow( as.numeric(t(a0_vec)), t_max) + reprow(idx_a_tv, t_max) * apply(matrix(rnorm(t_max * k_alp) * V_ab, nrow = t_max), MARGIN = 2, FUN = cumsum)
  H_t <- reprow( as.numeric(h), t_max)

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)
  eps <- matrix(rnorm(t_max*K), ncol = K)



  for (i in c(1:t_max)){
    inv_A0 <- solve(A0_mat(A_t[i,], K))
    Sigma_t <-  inv_A0 %*% diag(exp(0.5*H_t[i,]), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    xt <- rbind(1, matrixcalc::vec( t(ystar[(p+i-1):i,])))
    Bt_mat <- B0_mat(B_t[i,], K, p)
    y_mean[i,] <- inv_A0 %*% (Bt_mat %*% xt)
    ysim <-  y_mean[i,] + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       h = h, B_t = B_t, A_t = A_t, H_t = H_t,
       dist = "Gaussian", SV = FALSE)
}

#' @export
sim.tvpVAR.Student.nonSV <- function(K = 3, p = 2, t_max = 1000,
                                  b0 = 0.5, a0 = 0.5, h = -3, sigma_ab = NULL,
                                  y0 = matrix(0, ncol = K, nrow = p),
                                  nu = 6, is_tv = NULL, seednum = 0, burn_in = 0){
  t_max = t_max + burn_in
  set.seed(seednum)
  # Sample matrix coefficient B
  B0 <- B0_mat(b0, K, p)

  # Sample matrix corr A0
  A0 <- A0_mat(a0, K)
  tA0 <- t(A0)
  a0_vec <- tA0[upper.tri(tA0)]

  # Sample matrix variance h
  if (length(h) == 1){
    h <- rep(h,K)
  }
  # No skew
  # Tail of student
  w_t <- matrix(rinvgamma(t_max*K, shape = nu/2, rate = nu/2), ncol = K)
  w_sqrt_t <- sqrt(w_t)


  if (is.null(sigma_ab)){
    V_ab <- (2e-3)
  }

  if (is.null(is_tv)){
    is_tv <- rep(1,K)
  }

  # Identify the tvp element
  k_alp <- K*(K-1)/2 # dimension of the impact matrix
  k_beta <- K^2*p + K # number of VAR coefficients
  k_beta_div_K <- K*p + 1
  n_tv <- sum(is_tv)# # of time-varying equations
  idx_b_tv <- (kronecker(is_tv,matrix(1,nrow = K*p+1,ncol = 1))==1)   # index for time-varying betas
  idx_a_tv <- matrix(FALSE, nrow = k_alp, ncol = 1)            # construct index for time-varying alphas
  count <- 0
  for (j in 2:K){
    if (is_tv[j] == 1){
      idx_a_tv[ (count+1):(count+j-1)] <- TRUE
    }
    count <- count + j-1
  }


  B_t <- reprow( as.numeric(t(B0)), t_max) + reprow(idx_b_tv, t_max) * apply(matrix(rnorm(t_max * k_beta) * V_ab, nrow = t_max), MARGIN = 2, FUN = cumsum)
  A_t <- reprow( as.numeric(t(a0_vec)), t_max) + reprow(idx_a_tv, t_max) * apply(matrix(rnorm(t_max * k_alp) * V_ab, nrow = t_max), MARGIN = 2, FUN = cumsum)
  H_t <- reprow( as.numeric(h), t_max)

  ystar <- tail(y0, p)
  y_mean <- matrix(NA, nrow = t_max, ncol = K)
  y_var <- matrix(NA, nrow = t_max, ncol = 0.5*K*(K+1))
  volatility <- matrix(NA, nrow = t_max, ncol = K)

  eps <- matrix(rnorm(t_max*K), ncol = K)



  for (i in c(1:t_max)){
    inv_A0 <- solve(A0_mat(A_t[i,], K))
    Sigma_t <-  inv_A0 %*% diag(w_sqrt_t[i,] * exp(0.5*H_t[i,]), nrow = K)
    Sigma2_t <- Sigma_t %*% t(Sigma_t)
    y_var[i,] <- Sigma_t[lower.tri(Sigma_t, diag = T)]
    volatility[i,] <- diag(Sigma2_t)
    xt <- rbind(1, matrixcalc::vec( t(ystar[(p+i-1):i,])))
    Bt_mat <- B0_mat(B_t[i,], K, p)
    y_mean[i,] <- inv_A0 %*% (Bt_mat %*% xt)
    ysim <-  y_mean[i,] + Sigma_t %*% eps[i,]
    ystar <- rbind(ystar, t(ysim))
  }

  t_max = t_max - burn_in
  list(y = as.matrix(ystar[(p+burn_in+1):(p+burn_in+t_max),], nrow = t_max),
       y0 = y0, y_mean = y_mean, y_var = y_var, volatility = volatility,
       K = K, p = p, t_max = t_max,
       A0 = A0, B0 = matrix(B0, nrow = K),
       nu = nu, w = w_t[(burn_in+1):(burn_in+t_max)],
       h = h, B_t = B_t, A_t = A_t, H_t = H_t,
       dist = "Student", SV = FALSE)
}
