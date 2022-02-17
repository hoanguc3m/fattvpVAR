#' @export
ML_GaussTVPSVTest <- function(Chain, numCores = 4){
  #Chain <- G101_obj
  priors <- Chain$data$priors
  data = Chain$data
  Y0 <- data$y0
  shortY <- data$y
  t_max <- nrow(shortY)
  K <- ncol(shortY)
  Y <- matrixcalc::vec(t(shortY))
  p <- data$p
  is_tv <- data$inits$is_tv
  priors <- data$priors
  M <- 20000

  #is_tv = [1 1 1]'; # n-vector that denotes which equations have TVP; must have the same dim as the VAR
  # e.g., [1 0 1]: eq 1 and 3 have TVP, but eq 2 has const parameters

  k_alp <- K*(K-1)/2 # dimension of the impact matrix
  k_beta <- K^2*p + K # number of VAR coefficients
  k_beta_div_K <- K*p + 1
  n_tv <- sum(is_tv)# # of time-varying equations

  #n = size(store_Sigh,2); # It is our K

  ## compute and define a few things
  tmpY <- rbind(Y0, shortY)
  X2 <- matrix(0, nrow = t_max, ncol = K*p)
  for (ii in 1:p){
    X2[,((ii-1)*K+1):(ii*K)] <- tmpY[(p-ii+1):(t_max+p-ii),]
  }
  X2 <- cbind(rep(1, t_max), X2)

  idx_b_tv <- (kronecker(is_tv,matrix(1,nrow = K*p+1,ncol = 1))==1)   # index for time-varying betas
  idx_a_tv <- matrix(FALSE, nrow = k_alp, ncol = 1)            # construct index for time-varying alphas

  k_a_eq <- seq(0, K-1)
  id_a <- cbind(cumsum(k_a_eq) - k_a_eq + 1,cumsum(k_a_eq))
  count_seqa <- list(); count_seqb <- list()
  count_seqa[[1]] <- 0; count_seqb[[1]] <- seq(1, k_beta_div_K)
  for (ii in c(2:K)){
    count_seqa[[ii]] <- seq(id_a[ii,1], id_a[ii,2])
    count_seqb[[ii]] <- ((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)

    if (is_tv[ii] == 1){
      idx_a_tv[ count_seqa[[ii]]] <- TRUE
    }
  }

  # Sigma approx
  Sigma_h_list <- Gamma_approx(mcmc_sample = Chain$store_Sigh, ndraws = M)
  Sigma_h_gen <- Sigma_h_list$new_samples
  #sum_log_prop <- sum_log_prop + Sigma_h_gen$sum_log_prop


  Sigma_beta_gen <- matrix(0, ncol = ncol(Chain$store_Sigbeta), nrow = M) #
  if (sum(idx_b_tv) > 0 ){
    Sigma_beta_list <- Gamma_approx(mcmc_sample = Chain$store_Sigbeta[,idx_b_tv, drop = FALSE], ndraws = M)
    Sigma_beta_gen[,idx_b_tv] <- Sigma_beta_list$new_samples
  } else {
    Sigma_beta_list <- list(sum_log_prop = 0)
  }

  Sigma_alp_gen <- matrix(0, ncol = ncol(Chain$store_Sigalp), nrow = M) #
  if (sum(idx_a_tv) > 0 ){
    Sigma_alp_list <- Gamma_approx(mcmc_sample = Chain$store_Sigalp[,idx_a_tv, drop = FALSE], ndraws = M)
    Sigma_alp_gen[,idx_a_tv] <- Sigma_alp_list$new_samples
  } else {
    Sigma_alp_list <- list(sum_log_prop = 0)
  }

  # # beta0, alp0, h0 approx  and obtain IS draws
  # beta0_list <- Normal_approx(Chain$store_beta0, ndraws = M)
  # beta0_gen <- beta0_list$new_samples
  #
  # #alp0_gen <- Chain$store_alp0
  # alp0_list <- Normal_approx(Chain$store_alp0, ndraws = M)
  # alp0_gen <- alp0_list$new_samples

  # beta0, alp0, h0 approx  and obtain IS draws
  alpbeta0_list <- Normal_approx(cbind(Chain$store_beta0, Chain$store_alp0), ndraws = M)
  beta0_gen <- alpbeta0_list$new_samples[,1:k_beta]
  alp0_gen <- alpbeta0_list$new_samples[,(k_beta+1):(k_beta+k_alp)]

  #h0_gen <- Chain$store_h0
  h0_list <- Normal_approx(Chain$store_h0, ndraws = M)
  h0_gen <- h0_list$new_samples

  sum_log_prop <- Sigma_h_list$sum_log_prop + Sigma_beta_list$sum_log_prop + Sigma_alp_list$sum_log_prop +
    h0_list$sum_log_prop + alpbeta0_list$sum_log_prop # beta0_list$sum_log_prop + alp0_list$sum_log_prop


  ## prior borrow from inference
  abeta0 <- priors$b0
  Vbeta0 <- priors$V_b_prior
  aalp0 <- priors$a0
  Valp0 <- priors$V_a0_prior

  ah0 <- c(log(priors$sigma^2))
  Vh0 <- 4*matrix(1, nrow = K,ncol = 1)

  sum_log_prior <- apply( cbind(dgamma(Sigma_h_gen, shape = 0.5, rate = 0.5 / priors$hyper_h, log = T),
                                dgamma(Sigma_beta_gen[,idx_b_tv, drop = FALSE], shape = 0.5, rate = 0.5 / priors$hyper_ab, log = T),
                                dgamma(Sigma_alp_gen[,idx_a_tv, drop = FALSE], shape = 0.5, rate = 0.5 / priors$hyper_ab, log = T),
                                t(dnorm(t(h0_gen), mean = ah0, sd = sqrt(Vh0), log = T)),
                                t(dnorm(t(beta0_gen), mean = abeta0, sd = sqrt(Vbeta0), log = T)),
                                t(dnorm(t(alp0_gen), mean = aalp0, sd = sqrt(Valp0), log = T)) ),  MARGIN = 1, FUN = sum)

  #sum_log = rep(0, M);
  RhpcBLASctl::blas_set_num_threads(1)
  yt <- t(shortY)
  xt <- fatBVARS::makeRegressor(shortY, Y0, t_max, K, p)
  sum_log <- parallel::mclapply(1:M,
                                FUN = function(j) {
                                  Sigbeta = Sigma_beta_gen[j, ]
                                  Sigalp = Sigma_alp_gen[j, ]
                                  Sigh = Sigma_h_gen[j, ]
                                  beta0 = beta0_gen[j, ]
                                  alp0 = alp0_gen[j, ]
                                  h0 = h0_gen[j, ]


                                  A <- A0_mat(alp0_gen[j,], K)
                                  B <- solve(A) %*% matrix(beta0_gen[j, ], nrow = K, byrow = TRUE)

                                  sigma_h <- Sigma_h_gen[j, ]
                                  h0 <- h0_gen[j, ]
                                  ytilde <- as.numeric(A %*% (yt - B %*% xt)) # from dim K * t_max to (Kt_max) * 1



                                  llike <- fatBVARS::int_h_Gaussian(ytilde, h0 = h0, sigma_h = sigma_h, t_max = t_max, K=K, R = 100)
                                  llike

                                }, mc.cores = numCores)
  sum_log = unlist(sum_log)

  store_w = sum_log + sum_log_prior - sum_log_prop
  shortw =  matrix(store_w, nrow = M/20, ncol = 20)
  maxw = apply(shortw, MARGIN = 2 , FUN = max)
  bigml = log( apply(exp(shortw-reprow(maxw,M/20)), MARGIN = 2, FUN = mean )) + maxw
  ml = mean(bigml)
  mlstd = sd(bigml)/sqrt(20)
  return( list( LL = ml,
                std = mlstd,
                log_all = store_w,
                log_intll = sum_log,
                sum_log_prior = sum_log_prior,
                sum_log_prop = sum_log_prop ))
}
