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
