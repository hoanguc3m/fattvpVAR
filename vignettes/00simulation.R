library(readxl)
library(Matrix)
library(fattvpVAR)
library(profvis)
library(invgamma)

datagenG <- sim.tvpVAR.SV(dist="Gaussian", h = -3, sigma_ab = 5e-3, sigma_h = 3e-2, b0 = 0.5, a0 = 0.1, is_tv = c(1,1,1) )
matplot(datagenG$y, type = "l")
datagenS <- sim.tvpVAR.SV(dist="Student", h = -3, sigma_ab = 5e-3, sigma_h = 3e-2, b0 = 0.5, a0 = 0.1 )
plot(datagenS$y[,1])

# datagenG <- sim.tvpVAR.nonSV(dist="Gaussian")
# plot(datagenG$y[,1])
# datagenS <- sim.tvpVAR.nonSV(dist="Student")
# plot(datagenS$y[,1])
datagenG <- sim.tvpVAR.SV(dist="Student", h = -3, sigma_ab = 5e-3, sigma_h = 3e-2, b0 = 0.5, a0 = 0.1, nu = 6 )


y <- datagenG$y
y0 <- datagenG$y0

inits <- list(samples = 10000,
              burnin = 1000,
              thin = 1,
              is_tv = c(1,1,1) )
priors <- list(  hyper_ab = 1,
                 hyper_h = 1)

p <- 2 # number of lags
K <- 3
k_alp <- K*(K-1)/2 # dimension of the impact matrix
k_beta <- K^2*p + K # number of VAR coefficients
k_beta_div_K <- K*p + 1
is_tv <- inits$is_tv
n_tv <- sum(is_tv)# # of time-varying equations

##########################################################################
beta <- datagenG$B_t
alp <- datagenG$A_t
h <- datagenG$H_t

Sigbeta <- matrix(0.01^2, nrow = k_beta, ncol = 1)
Sigalp <- matrix(0.01^2, nrow = k_alp, ncol = 1)
Sigh <- matrix( 0.01^2, nrow = K, ncol = 1)

beta0 <- matrix(t(datagenG$B0), ncol = 1)
alp0 <- matrix(t(datagenG$A0)[upper.tri(datagenG$A0)], ncol = 1)
h0 <- matrix(rep(-3,K), ncol = 1)

w <- datagenG$w
nu <- datagenG$nu

##########################################################################
tmp <- apply(store_beta, MARGIN = c(2,3), FUN = mean)
ii <- 6
plot(tmp[,ii], type = "l", col = "red", ylim = c(-1,1))
lines(datagenG$B_t[,ii])

lines(datagenG$H_t[,ii])
lines(h[,ii], col = "blue")

U <- matrix(0, nrow = t_max, ncol = K)
E <- matrix(0, nrow = t_max, ncol = K)
count <- 0


for (ii in 1:K){
  ki <- K*p+1+ii-1

  if (ii > 1) {
    X <- cbind(X2, -shortY[,1:(ii-1)])
  } else {
    X <- X2
  }

  if (is_tv[ii] == 1){

    bigXi <- SURform(X)
    # % S = sparse(i,j,s,m,n,nzmax) uses vectors i, j, and s to generate an
    # %   m-by-n sparse matrix such that S(i(k),j(k)) = s(k), with space
    # %   allocated for nzmax nonzeros.
    if (ii == 1){
      count_seq <- 0
    } else {
      count_seq <- (count+1):(count+ii-1)
    }

    thetai0 <- cbind( beta[,((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)], alp[,count_seq] )
    thetai <- vec(t(thetai0))
  }


  count <- count + ii-1
  U[,ii] <- as.matrix(shortY[,ii] - bigXi %*% thetai)
  E[,ii] <- U[,ii] / sqrt(w[,ii])
}

