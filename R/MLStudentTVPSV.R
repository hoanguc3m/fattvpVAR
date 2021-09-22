#' @export
ML_StudentTVPSV <- function(Chain, numCores = 4){
  #Chain <- T000_obj
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
  # X1 <- matrix(0, K*t_max,k_alp)
  # count <- 0
  # for (ii in 2:K){
  #   X1[seq(ii,K*t_max,by = K),(count+1):(count+ii-1)] <- - shortY[,1:ii-1]
  #   count <- count + ii-1
  # }

  idx_b_tv <- (kronecker(is_tv,matrix(1,nrow = K*p+1,ncol = 1))==1)   # index for time-varying betas
  idx_a_tv <- matrix(FALSE, nrow = k_alp, ncol = 1)            # construct index for time-varying alphas
  count <- 0
  for (ii in 2:K){
    if (is_tv[ii] == 1){
      idx_a_tv[ (count+1):(count+ii-1)] <- TRUE
    }
    count <- count + ii-1
  }

  # Sigma approx
  Sigma_h_gen <- Chain$store_Sigh # Always SV
  Sigma_h_list <- Normal_approx(Chain$store_Sigh^0.5, ndraws = M) # Change to normal
  Sigma_h_gen <- Sigma_h_list$new_samples^2 # Change to square

  Sigma_beta_gen <- matrix(0, ncol = ncol(Chain$store_Sigbeta), nrow = M) #
  if (sum(idx_b_tv) > 0 ){
    Sigma_beta_list <- Normal_approx(Chain$store_Sigbeta[,idx_b_tv, drop = FALSE]^0.5, ndraws = M) # Change to normal
    Sigma_beta_gen[,idx_b_tv] <- Sigma_beta_list$new_samples^2 # Change to square
  } else {
    Sigma_beta_list <- list(sum_log_prop = 0)
  }

  Sigma_alp_gen <- matrix(0, ncol = ncol(Chain$store_Sigalp), nrow = M) #
  if (sum(idx_a_tv) > 0 ){
    Sigma_alp_list <- Normal_approx(Chain$store_Sigalp[,idx_a_tv, drop = FALSE]^0.5, ndraws = M) # Change to normal
    Sigma_alp_gen[,idx_a_tv] <- Sigma_alp_list$new_samples^2 # Change to square
  } else {
    Sigma_alp_list <- list(sum_log_prop = 0)
  }

  # beta0, alp0, h0 approx  and obtain IS draws
  #beta0_gen <- Chain$store_beta0
  beta0_list <- Normal_approx(Chain$store_beta0, ndraws = M)
  beta0_gen <- beta0_list$new_samples

  #alp0_gen <- Chain$store_alp0
  alp0_list <- Normal_approx(Chain$store_alp0, ndraws = M)
  alp0_gen <- alp0_list$new_samples

  #h0_gen <- Chain$store_h0
  h0_list <- Normal_approx(Chain$store_h0, ndraws = M)
  h0_gen <- h0_list$new_samples

  Nu_gen_list <- Nu_Gamma_approx(Chain$store_nu, ndraws = M)
  Nu_gen <- Nu_gen_list$new_samples
  w_mean <- apply(Chain$store_w, MARGIN = c(2,3), mean)


  sum_log_prop <- Sigma_h_list$sum_log_prop + Sigma_beta_list$sum_log_prop + Sigma_alp_list$sum_log_prop +
    h0_list$sum_log_prop + beta0_list$sum_log_prop + alp0_list$sum_log_prop + Nu_gen_list$sum_log_prop


  ## prior borrow from inference
  abeta0 <- matrix(0, nrow = k_beta, ncol = 1);
  Vbeta0 <- 10*matrix(1, nrow = k_beta, ncol = 1);
  aalp0 <- matrix(0, nrow = k_alp, ncol = 1);
  Valp0 <- 10*matrix(1, nrow = k_alp, ncol = 1);
  nu_gam_a <- 2; nu_gam_b <- 0.1;

  EstMdl1 <- arima(shortY[,1] ,order = c(1,0,0))
  EstMdl2 <- arima(shortY[,2] ,order = c(1,0,0))
  EstMdl3 <- arima(shortY[,3] ,order = c(1,0,0))

  ah0 <- c(log(EstMdl1$sigma2), log(EstMdl1$sigma2), log(EstMdl1$sigma2))
  Vh0 <- 4*matrix(1, nrow = K,ncol = 1)


  sum_log_prior <- apply( cbind(dgamma(Sigma_h_gen, shape = 0.5, rate = 0.5 / priors$hyper_h, log = T),
                           dgamma(Sigma_beta_gen[,idx_b_tv], shape = 0.5, rate = 0.5 / priors$hyper_ab, log = T),
                           dgamma(Sigma_alp_gen[,idx_a_tv], shape = 0.5, rate = 0.5 / priors$hyper_ab, log = T),
                           dnorm(h0_gen, mean = ah0, sd = sqrt(Vh0), log = T),
                           dnorm(beta0_gen, mean = abeta0, sd = sqrt(Vbeta0), log = T),
                           dnorm(alp0_gen, mean = aalp0, sd = sqrt(Valp0), log = T),
                           dgamma(Nu_gen, shape = nu_gam_a, rate = nu_gam_b, log = T)),  MARGIN = 1, FUN = sum)

  #sum_log = rep(0, M);
  RhpcBLASctl::blas_set_num_threads(1)

  sum_log <- parallel::mclapply(1:M,
                     FUN = function(j) {
    Sigbeta = Sigma_beta_gen[j, ]
    Sigalp = Sigma_alp_gen[j, ]
    Sigh = Sigma_h_gen[j, ]
    beta0 = beta0_gen[j, ]
    alp0 = alp0_gen[j, ]
    h0 = h0_gen[j, ]
    nu = Nu_gen[j, ]

    count = 0
    llike = 0
    for (ii in c(1:K)){
      if (ii == 1){
        count_seq <- 0
        X <- X2
      } else {
        count_seq <- (count+1):(count+ii-1)
        X <- cbind(X2, -shortY[,1:(ii-1)])
      }


      if (is_tv[ii] == 1){
        bigXi = SURform(X)

        Sigthetai = c( Sigbeta[ ((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)], Sigalp[count_seq])
        thetai0 = c( beta0[((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)], alp0[count_seq] )
        llikei = intlike_tvpsv(Yi = shortY[,ii], Sigthetai = Sigthetai, Sig_hi = Sigh[ii], bigXi = bigXi, h0i = h0[ii], thetai0 = thetai0)
      }  else {
        bigXi = X
        thetai = c( beta0[((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)], alp0[count_seq] )
        llikei = intlike_Tvarsv(Yi = shortY[,ii], thetai = thetai, Sig_hi = Sigh[ii], bigXi = bigXi, h0i = h0[ii], nui = nu[ii], wi = w_mean[,ii]);
        # tmp1 <- sapply(1:200, FUN = function(o){ intlike_Tvarsv1(Yi = shortY[,ii], thetai = thetai, Sig_hi = Sigh[ii], bigXi = bigXi, h0i = h0[ii], nui = nu[ii], wi = w_mean[,ii]);})
        # tmp2 <- sapply(1:200, FUN = function(o){ intlike_Tvarsv2(Yi = shortY[,ii], thetai = thetai, Sig_hi = Sigh[ii], bigXi = bigXi, h0i = h0[ii], nui = nu[ii], wi = w_mean[,ii]);})
        # tmp3 <- sapply(1:200, FUN = function(o){ intlike_Tvarsv3(Yi = shortY[,ii], thetai = thetai, Sig_hi = Sigh[ii], bigXi = bigXi, h0i = h0[ii], nui = nu[ii], wi = w_mean[,ii]);})
        # c(mean(tmp1), mean(tmp2), mean(tmp3))
        # c(sd(tmp1), sd(tmp2), sd(tmp3))
        # #llikei = intlike_varsv(Yi = shortY[,ii], thetai = thetai, Sig_hi = Sigh[ii], bigXi = bigXi, h0i = h0[ii]);
      }

      llike = llike + llikei
      count = count + ii-1

    }

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
                std = mlstd))
}
