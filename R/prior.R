#' @export
get_prior_minnesota <- function(y, p, intercept=TRUE, ...){
  # lambda1=0.2, lambda2=0.5, lambda3=1, lambda4=2
  # Chan paper on large hybrid VAR
  arguments <- eval(substitute(alist(...)))
  lambda1 <- ifelse(is.null(arguments$lambda1), 0.2, arguments$lambda1)
  lambda2 <- ifelse(is.null(arguments$lambda2), 0.5, arguments$lambda2)
  lambda3 <- ifelse(is.null(arguments$lambda3), 1, arguments$lambda3)
  lambda4 <- ifelse(is.null(arguments$lambda4), 10, arguments$lambda4)
  hyper_ab <- 1
  hyper_h <- 1

  t_max <- nrow(y)
  K   <- ncol(y)

  # Prior b0 mean
  constant = 0
  B0 <- NULL
  if(intercept == TRUE) {
    constant = 1
    B0 <- cbind(rep(0,K))
  }

  B0 <- cbind(B0, diag(K))
  if (p > 1){
    for (i in c(2:p)){
      B0 <- cbind(B0, matrix(0,K,K))
    }
  }
  b0 <- as.numeric(t(B0)) # B0 by row

  # Prior for covariance matrix

  sigmasq <- rep(0,K)
  for(ii in 1:K){
    EstMdl <- arima(y[,ii] ,order = c(p,0,0))
    sigmasq[ii] <- EstMdl$sigma2
  }


  # Covariance matrix for the prior

  M <- K * p + constant
  Vi <- rep(0, K * M)

  # Add Covariance coefficients for intercepts
  # Covariance for the intercept
  if (intercept==TRUE){
    for(ii in 1:K){
      Vi[ii]    <- lambda4 * sigmasq[ii]
    }
  }

  # without intercept
  for(jj in 1:p){ #loop over the jj-th lag
    for(kk in 1:K){ #kk-th variable
      for(ii in 1:K){ # loop over the ii-th equation
        indx <- (kk - 1) * K  + (jj - 1) * K * K + ii + constant * K
        if(ii==kk){
          Vi[indx] <- lambda1/jj^2
        } else{
          Vi[indx] <- (lambda1 * lambda2)/(jj ^2) * (sigmasq[ii] / sigmasq[kk])
        }
      }
    }
  }

  V_b_prior = as.numeric( t(matrix(Vi, nrow = K)))

  # Prior a0
  a0 <- rep(0, 0.5*K*(K-1))
  V_a0_prior <- lambda3 * crossprod(t(sigmasq), sigmasq^-1)
  V_a0_prior <- as.numeric(V_a0_prior[lower.tri(V_a0_prior)])

  pr <- list(type = "Minnesota",
             p = p,
             intercept = intercept,
             b0 = b0,
             V_b_prior = V_b_prior,
             a0 = a0,
             V_a0_prior = V_a0_prior,
             sigma = sqrt(sigmasq),
             hyper_ab = hyper_ab,
             hyper_h = hyper_h)
  pr$nu_gam_a = 2
  pr$nu_gam_b = 0.1

  return(pr)
}

######################################################################################

#' Prior distribution of BVAR model
#'
#' This function returns a prior specification of BVAR-SV-fatTail model.
#' @param y The input data as matrix T x K where T is the number of observations and K is the number of variables
#' @param p The number of lags in BVAR model.
#' @param priorStyle The prior style in BVAR model should be c("Minnesota", "OLS")
#' @param dist The variable specifies the BVAR error distribution. It should be one of
#' c("Gaussian","Student", "Skew.Student", "MT","MST","orthoStudent").
#' @param SV The indicator if this BVAR model has a Stochastic volatility part.
#' @return A list of prior specifications. \eqn{b \sim N(b0, V_b_prior)}, \eqn{a \sim N(0, 1000I)}, \eqn{sigmaSq \sim IG(0.5*sigma_T0, 0.5*sigma_S0)},
#' \eqn{nu \sim G(a,b)}, \eqn{gamma \sim N(0, I)}.
#' @export
#' @examples
#' \dontrun{
#' prior <- get_prior(y, p = 2, dist = c("Student"), SV = FALSE)
#' }
# get_prior <- function(y, p, priorStyle = c("Minnesota"),
#                       dist = c("Gaussian"),
#                       SV = FALSE, ...){
#   if (!(dist %in% c("Gaussian","Student") ))
#     stop("dist is not implemented.")
#
#   arguments <- eval(substitute(alist(...)))
#
#   K <- ncol(y)
#   M <- K + p*(K^2) # nr of beta parameters
#   numa <- 0.5*K*(K-1) # nr of VCV elements
#
#   # return b0, V_b_prior, sigma
#   if (priorStyle == "Minnesota"){
#     prior_sub <- get_prior_minnesota(y, p, intercept = TRUE, ...)
#   }
#
#   # B mean and var
#   b0 = prior_sub$b0
#   V_b_prior = prior_sub$V_b_prior
#
#   # A mean and var
#   a0 <- prior_sub$a0
#   a_Vprior <- prior_sub$V_a0_prior
#   V_a_prior <- diag(x = rep(a_Vprior,numa), nrow = numa, ncol = numa)
#
#   sigma <- prior_sub$sigma
#
#   #Gaussian
#   prior_collect <- list(b_prior=b0, V_b_prior = V_b_prior,
#                         a_prior = a0, V_a_prior = V_a_prior,
#                         sigma = sigma,
#                         sigma_T0 = 1, sigma_S0 = 1)
#   #Student
#   if (dist !="Gaussian"){
#     prior_collect$nu_gam_a = 2
#     prior_collect$nu_gam_b = 0.1
#   }
#
#   prior_collect$t_max <- nrow(y)
#   prior_collect$dist <- dist
#   prior_collect$SV <- SV
#   prior_collect$priorStyle <- priorStyle
#
#   if (SV) {
#     prior_collect$sigma_h <- 1;
#   }
#   return(prior_collect)
# }

