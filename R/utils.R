# > X
#       [,1] [,2]
# [1,]    1    4
# [2,]    2    5
# [3,]    3    6
#
# > Xout
# 3 x 6 sparse Matrix of class "dgCMatrix"
#
# [1,] 1 4 . . . .
# [2,] . . 2 5 . .
# [3,] . . . . 3 6
#' @export
SURform <- function(X){
  r <- nrow(X); c <- ncol(X);
  idi <- kronecker(c(1:r), rep(1,c));
  idj <- c(1: (r*c))
  Xout <- Matrix::sparseMatrix(i = idi,j = idj, x = as.numeric(t(X)))
  return(Xout)
}

#' @export
reprow<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#' @export
repcol<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#' @export
adaptamount <- function(iteration){
  return( min( 0.01, 1.0 / sqrt(iteration) ) );
}

#' @export
B0_mat <- function(b0, K, p){
  # matrix coefficient B
  B0 <- cbind(rep(0,K))
  if (length(b0) == 1) {
    for (i in c(1:p)){
      B0 <- cbind(B0, b0^i*diag(K))
    }
  } else {
    B0 <- matrix(b0, nrow = K, byrow = TRUE)
  }
  return(B0)
}


#' @export
A0_mat <- function(a0, K){
  # Sample matrix corr A0
  if (length(a0) == 1) {
    A0 <- matrix(a0, K, K)
    diag(A0) <- 1
    A0[upper.tri(A0)] <- 0
  } else {
    A0 <- matrix(0, nrow = K, ncol = K)
    A0[upper.tri(A0)] <- a0
    A0 <- t(A0)
    diag(A0) <- 1
  }
  return(A0)
}

#' @export
Sigma_sample <- function(Beta, Beta0, Prior_Beta, t_max){
  K <- ncol(Beta)
  if (K>1) {
    sse_2 <- apply( (Beta - rbind(Beta0,Beta[1:(t_max-1),]) )^2, MARGIN = 2, FUN = sum)
  } else {
    sse_2 <- sum( (Beta - c(Beta0,Beta[1:(t_max-1),]) )^2 )
  }

  # # Normal prior
  # # Equation 9 in https://doi.org/10.1016/j.csda.2013.01.002
  # sigma_post_a <- rep(t_max,K) # prior of sigma_h Gamma(1,0.0001)
  # sigma_post_b <- sse_2 # prior of sigma_h
  #
  # for (i in c(1:K)){
  #   sigma_new <- rinvgamma(1, shape = sigma_post_a[i] * 0.5, rate = sigma_post_b[i] * 0.5)
  #   alpha = 0.5 * (Sigma_Beta[i] - sigma_new) / Prior_Beta[i] + 0.5 * (log(sigma_new) - log(Sigma_Beta[i])) # B_sigma = 1
  #   temp = log(runif(1))
  #   if (alpha > temp){
  #     Sigma_Beta[i] <- sigma_new
  #   }
  #   #log_sigma_den[]
  # }

  sigma_h <- mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                          psi = 1/Prior_Beta )

  return(sigma_h)
}

# This function is borrowed from https://github.com/FK83/bvarsv
#' @export
getmix <- function(){
  # 7 components
  q <- c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)      # probabilities
  m <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) # means
  u2 <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)    #variances

  # # 10 components
  # q <- c(0.00609, 0.04775, 0.13057, 0.20674, 0.22715, 0.18842, 0.12047, 0.05591, 0.01575, 0.00115)      # probabilities
  # m <- c(3.13417, 2.55484, 1.94244, 1.23006, 0.35567, -0.76538, -2.26048, -4.34506, -7.47644, -13.44260)
  # u2 <- c(0.11265, 0.17788, 0.26768, 0.40611, 0.62699, 0.98583, 1.57469, 2.54498, 4.16591, 7.33342)    #variances
  # # m <- c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000) # means in the Omori paper
  return(list(q=q,m=m,u2=u2))
}

#' @export
sample_h_ele <- function(ytilde, sigma_h = 0.0001*diag(K), h0_mean = rep(0,K),
                         h = matrix(0, nrow = t_max, ncol = K), K, t_max, prior){
  tmp <- getmix()
  q <- tmp$q
  m_mean <- tmp$m
  u2 <- tmp$u2


  {
    h0 <- rep(0,K)
    cond_var_sigma <- rep(0,K)
    cond_mean_sigma <- rep(0,K)
    Zs <- matrix(1,t_max,1) %x% diag(1)
    for (i in c(1:K)){
      sigma_prmean <- h0_mean[i] # mean h_0
      sigma_prvar <- matrix(4)   # variance h_0
      aux <- sigmahelper4(t(ytilde[ i,, drop =FALSE]^2), q, m_mean, u2, h[ i,, drop =FALSE], Zs, matrix(sigma_h[i]), sigma_prmean, sigma_prvar)
      h[i,] <- aux$Sigtdraw
      h0[i] <- as.numeric(aux$h0)
    }
  }

  # sqrtvol <- aux$sigt
  # [TODO] fix this
  # sse_2 <- apply( (h[,2:t_max] - h[,1:(t_max-1)])^2, MARGIN = 1, FUN = sum)
  if (K>1) {
    sse_2 <- apply( (h[,1:t_max] - cbind(h0,h[,1:(t_max-1)]) )^2, MARGIN = 1, FUN = sum)
  } else {
    sse_2 <- sum( (h[,1:t_max] - c(h0,h[,1:(t_max-1)]) )^2)
  }


  sigma_h <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                          psi = 1/prior$hyper_h ) , nrow = K)


  aux <- list(sigma_h = sigma_h,
              h0 = h0,
              Sigtdraw = h,
              sigt = exp(0.5*h))
  return(aux)
}

#' @export
Normal_approx <- function(mcmc_sample, ndraws){
  mcmc_mean <- apply(mcmc_sample, 2, mean)
  mcmc_Sigma <- cov(mcmc_sample)
  nElements <- length(mcmc_mean)

  new_samples <- mvnfast::rmvn(ndraws, mu = mcmc_mean, sigma = mcmc_Sigma)
  colnames(new_samples) <- colnames(mcmc_sample)
  sum_log_prop <- mvnfast::dmvn(X = new_samples,
                                mu = as.numeric(mcmc_mean), sigma = mcmc_Sigma, log = T)
  return(list(new_samples = new_samples,
              sum_log_prop = sum_log_prop))
}

#' @export
LL_tnorm <- function(par, data){
- sum(log(dtruncnorm(x = data, a=0, b=Inf, mean = par[1], sd = par[2])))
}

#' @export
Normal_trunc_approx <- function(mcmc_sample, ndraws){

  nElements <- ncol(mcmc_sample)
  mcmc_mean <- rep(0, nElements)
  mcmc_Sigma <- rep(0, nElements)
  new_samples <- matrix(NA, ncol = nElements, nrow = ndraws)
  sum_log_prop <- rep(0, ndraws)

  for (i in c(1:nElements)){
    result <- optim(par = c(0, 1), fn = LL_tnorm, data = mcmc_sample[,i])
    new_samples[,i] <- rtruncnorm(ndraws, mean = result$par[1], sd = result$par[2], a = 0)
    sum_log_prop <- sum_log_prop +
      log(dtruncnorm(x = new_samples[,i], a=0, b=Inf, mean = result$par[1], sd = result$par[2])) - log(2) - log(new_samples[,i]) #Jacobian trans
  }

  colnames(new_samples) <- colnames(mcmc_sample)
  return(list(new_samples = new_samples,
              sum_log_prop = sum_log_prop))
}


#' @export
Nu_Gamma_approx <- function(mcmc_sample, ndraws){
  nElements <- ncol(mcmc_sample)
  new_samples <- matrix(NA, ncol = nElements, nrow = ndraws, dimnames = list(c(), colnames(mcmc_sample)))
  Density_prop <-  matrix(NA, ncol = nElements, nrow = ndraws)
  shape_param <- rep(0, nElements)
  rate_param <- rep(0, nElements)
  for (i in c(1:nElements)){
    fit.gamma <- fitdistrplus::fitdist(as.numeric(mcmc_sample[,i]), distr = "gamma", method = "mle")
    shape_param[i] <- fit.gamma$estimate[1]
    rate_param[i] <- fit.gamma$estimate[2]
    tmp <- rgamma(ndraws*2, shape = shape_param[i], rate_param[i])
    tmp <- tmp[tmp > 2]
    tmp <- tmp[tmp < 100]
    new_samples[,i] <- head(tmp,ndraws)
    Density_prop[,i] <- dgamma(new_samples[,i], shape = shape_param[i], rate = rate_param[i], log = T)
  }
  return(list(new_samples = new_samples,
              sum_log_prop = apply(Density_prop, 1, sum)))
}

#' @export
int_w_MultiOrthStudent2 <- function(y, xt, A, B, sigma, nu, gamma = rep(0,K), t_max, K, R = 100){
  u <- A %*% (t(y) - B %*%xt)
  llw <- rep(0,K)
  for (j in c(1:K)){
    llw[j] <- sum(ghyp::dghyp(u[j,], object = ghyp::ghyp(lambda = -0.5*nu[j], chi = nu[j], psi = 0, mu = 0, sigma = sigma[j],
                                                         gamma = gamma[j]), logvalue = TRUE)) # input here sigma, not sigma^2
  }
  return(sum(llw))
}

#' @export
get_post <- function(G_obj, element){
  if (element =="beta"){
    store_beta <- G_obj$store_beta
    store_beta0 <- G_obj$store_beta0
    store_Sigbeta <- G_obj$store_Sigbeta
    for (i in c(1:dim(G_obj$store_alp)[1])){
      store_beta[i,,] <- t(store_beta0[i,] + sqrt(store_Sigbeta[i,]) *  t(store_beta[i,,]))
    }
    return(store_beta)
  }

  if (element =="alpha"){
    store_alp <- G_obj$store_alp
    store_alp0 <- G_obj$store_alp0
    store_Sigalp <- G_obj$store_Sigalp
    for (i in c(1:dim(G_obj$store_alp)[1])){
      store_alp[i,,] <- t(store_alp0[i,] + sqrt(store_Sigalp[i,]) *  t(store_alp[i,,]))
    }
    return(store_alp)
  }

  if (element =="h"){
    return(G_obj$store_h)
  }
  if (element =="vol"){
    return(exp(G_obj$store_h))
  }
}

#' @export
get_last_post <- function(G_obj, element){
  if (element =="beta"){
    store_beta <- G_obj$store_beta[,dim(G_obj$store_beta)[2],]
    store_beta0 <- G_obj$store_beta0
    store_Sigbeta <- G_obj$store_Sigbeta
    store_beta <- store_beta0 + sqrt(store_Sigbeta) *  store_beta
    return(store_beta)
  }

  if (element =="alpha"){
    store_alp <- G_obj$store_alp[,dim(G_obj$store_alp)[2],]
    store_alp0 <- G_obj$store_alp0
    store_Sigalp <- G_obj$store_Sigalp
    store_alp <- store_alp0 + sqrt(store_Sigalp) *  store_alp

    return(store_alp)
  }
}


stability <- function(b, K, p){
  # Dimensions
  B <- matrix(b, nrow = K, byrow = TRUE)
  B <- B[,c(2:(K*p+1))]

  if (p > 1){
    # VAR matrices
    Bc <- matrix(0, K*p, K*p)
    Bc[1:K, ] <- B
    Bc[-(1:K), 1:(K*(p-1))] <- diag(K*(p-1))

  } else {
    Bc <- B
  }
  ee = max(abs(eigen(Bc)$values))
  return(ee<1)
}


check_stability <- function(store_beta, store_alpha, p, K){
  ndraws <- nrow(store_beta)
  stabi <- rep(FALSE, ndraws)

  for (i in c(1:ndraws)){
      inv_A0 <- solve(A0_mat(store_alpha[i,], K))
      B <- inv_A0 %*% matrix(store_beta[i,], nrow = K, byrow = TRUE) #b_mat[i,]
      Bc <- matrix(0, K*p, K*p)
      Bc[1:K, ] <- B[,c(2:(K*p+1))]
      Bc[-(1:K), 1:(K*(p-1))] <- diag(K*(p-1))
      stabi[i] <- max(abs(eigen(Bc)$values)) <1
  }
  #sum(stabi)
  return(stabi)
}


#' @export
Gamma_approx <- function(mcmc_sample, ndraws){

  # MASS::fitdistr(1/Sigma_mat[,1], "gamma")
  # MASS::fitdistr(1/Sigma_mat[,3]/1000, "gamma")
  # MASS::fitdistr(Sigma_mat[,3], "exponential")

  nElements <- ncol(mcmc_sample)
  new_samples <- matrix(NA, ncol = nElements, nrow = ndraws, dimnames = list(c(), colnames(mcmc_sample)))
  Density_prop <-  matrix(NA, ncol = nElements, nrow = ndraws)
  shape_param <- rep(0, nElements)
  rate_param <- rep(0, nElements)
  mcmc_mean <- apply(mcmc_sample, 2, mean)
  mcmc_sd <- apply(mcmc_sample, 2, sd)

  for (i in c(1:nElements)){
    if (mcmc_sd[i] > 0){
      fit.gamma <- tryCatch({
        fitdistrplus::fitdist(as.numeric(mcmc_sample[,i]), distr = "gamma", method = "mle")
      }, error = function(e) {
        fitdistrplus::fitdist(as.numeric(mcmc_sample[,i]), distr = "exp", method = "mle")
      })
      if (fit.gamma$distname == "gamma"){
        shape_param[i] <- fit.gamma$estimate[1]
        rate_param[i] <- fit.gamma$estimate[2]
        new_samples[,i] <- rgamma(ndraws, shape = shape_param[i], rate_param[i])
        Density_prop[,i] <- dgamma(new_samples[,i],
                                   shape = shape_param[i],
                                   rate = rate_param[i], log = T)
      }
      if (fit.gamma$distname == "exp"){
        rate_param[i] <- fit.gamma$estimate[1]
        new_samples[,i] <- rexp(ndraws, rate =  rate_param[i])
        Density_prop[,i] <- dexp(new_samples[,i], rate = rate_param[i], log = T)

      }

    } else {
      new_samples[,i] <- as.numeric(mcmc_sample[,i])
      Density_prop[,i] <- 0
    }

  }

  return(list(new_samples = new_samples,
              sum_log_prop = apply(Density_prop, 1, sum)))
}
