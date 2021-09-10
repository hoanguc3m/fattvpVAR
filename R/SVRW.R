# This function implements the auxiliary mixture sampler to draw the
# log-volatility
#
# See:
# Chan, J.C.C. and Eisenstat, E. (2018). Comparing Hybrid Time-Varying
# Parameter VARs, Economics Letters, 171, 1-5.
#
#' @export
SVRW <- function(Ystar, h, sig, h0){
  t_max <- length(h)
  ## normal mixture
  pi <- c(0.0073, .10556, .00002, .04395, .34001, .24566, .2575)
  mi <- c(-10.12999, -3.97281, -8.56686, 2.77786, .61942, 1.79518, -1.08819) - 1.2704 ## means already adjusted!! ##
  sigi <- c(5.79596, 2.61369, 5.17950, .16735, .64009, .34023, 1.26261)
  sqrtsigi <- sqrt(sigi)

  ## sample S from a 7-point distrete distribution
  temprand <- matrix(runif(t_max), ncol = 1)
  q <- reprow(pi,t_max) * dnorm( repcol(Ystar,7), mean = repcol(h,7)+reprow(mi,t_max), sd = reprow(sqrtsigi,t_max))
  q <- q / repcol(rowSums(q),7)
  S <- matrix(7 - rowSums(repcol(temprand,7) <  t(apply(q, MARGIN = 1, cumsum)) )+1, ncol = 1)

  ## sample h
  # y^* = h + d + \epison, \epison \sim N(0,\Omega),
  # Hh = \alpha + \nu, \nu \ sim N(0,S),
  # where d_t = Ez_t, \Omega = diag(\omega_1,\ldots,\omega_n),
  # \omega_t = var z_t, S = diag(sig, \ldots, sig)
  # A = spdiags(B,d,m,n) which creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.

  Hh <-  Diagonal(t_max) - bandSparse(t_max,t_max, -1,list(rep(1,t_max-1)))
  invSh <- Diagonal(t_max,1/sig)
  dconst <- mi[S];

  invOmega <- 1./sigi[S] # diagonal
  alph <- Matrix::solve(Hh, sparseMatrix(i = 1, j = 1, x = h0, dims = c(t_max,1)))
  Kh <- Matrix::t(Hh) %*% invSh %*% Hh
  Ph <- Kh + Diagonal(t_max, invOmega)
  Ch <- Matrix::chol(Ph)
  hhat <- Matrix::solve(Ph, (Kh %*% alph + invOmega * (Ystar-dconst)))
  h <- hhat + Matrix::solve(Ch, rnorm(t_max))

  # invOmega <- Diagonal(t_max, 1./sigi[S])
  # alph <- solve(Hh, sparseMatrix(i = 1, j = 1, x = h0, dims = c(t_max,1)))
  # Kh <- t(Hh) %*% invSh %*% Hh
  # Ph <- Kh + invOmega
  # Ch <- chol(Ph)
  # hhat <- solve(Ph, (Kh %*% alph + invOmega %*% (Ystar-dconst)))
  # h <- hhat + solve(Ch, rnorm(t_max))
  return(as.matrix(h))
}
