
#' Forecast future observations of tvpVAR model
#'
#' This function returns a forecast future observations of TVP-VAR-SV model.
#' @param Chain The TVP-VAR-SV object from command fitTVPGaussSV/fitTVPStudentSV
#' @param t_pred The time prediction horizon.
#' @return The list of forecast future observations of BVAR model
#' @export
#' @examples
#' \dontrun{
#' forecast1 <- get_forecast(Chain1)
#' }
#'
get_forecast <- function(Chain, y0 = NULL, t_pred = 24, Nfsample = NULL){
  K <- ncol(y0)
  p <- Chain$data$p

  dist <- Chain$data$dist

  ndraws <- nrow(Chain$store_h0)

  if (is.null(Nfsample) || (Nfsample == ndraws) ){
    Nfsample <- ndraws
  }

  store_beta <- get_last_post(Chain, element = "beta")
  store_alpha <- get_last_post(Chain, element = "alpha")
  stab_id <- check_stability(store_beta, store_alpha, p, K)


  frange <- sample( (1:ndraws)[stab_id], Nfsample, replace = T)


  if (is.null(y0)){
    y0 <- as.matrix(tail(Chain$data$y,p))
  }

  y_pred <- array(NA, dim = c(t_pred, K, Nfsample))
  ymean_pred <- array(NA, dim = c(t_pred, K, Nfsample)) # Predictive mean
  yvar_pred <- array(NA, dim = c(t_pred, 0.5*K*(K+1), Nfsample)) # Predictive Chol Sigma
  vol_pred <- array(NA, dim = c(t_pred, K, Nfsample)) # Predictive diag(Sigma)


  count_id <- 0
  seed_set <- sample.int(n = 50000, size = Nfsample, replace = F)

  store_h <- Chain$store_h
  store_h0 <- Chain$store_h0
  store_Sigh <- sqrt(Chain$store_Sigh) # get sd

  store_Sigalp <- sqrt(Chain$store_Sigalp) # get sd
  store_Sigbeta <- sqrt(Chain$store_Sigbeta) # get sd

  store_nu <- Chain$store_nu


  for (i in frange){
    count_id <- count_id+1

    h <- tail(store_h[i,,],1)
    Sigh <- store_Sigh[i,]

    alp <- store_alpha[i,]
    Sigalp <- store_Sigalp[i,]

    beta <- store_beta[i,]
    Sigbeta <- store_Sigbeta[i,]

    if (dist == "Gaussian") {
        b0 <- beta; a0 <- alp; h <- h; sigma_ab <- c(Sigbeta, Sigalp); sigma_h <- Sigh; t_max = t_pred
        pred_tmp <- sim.tvpVAR.Gaussian.SV(K = K, p = p, t_max = t_pred,
                                           b0 = beta, a0 = alp, h = h,
                                           sigma_ab = c(Sigbeta, Sigalp), sigma_h = Sigh,
                                           y0 = y0, is_tv = is_tv, seednum = seed_set[count_id], burn_in = 0)


      }

    if (dist == "Student") {
        nu <- store_nu[i,]
        b0 <- beta; a0 <- alp; h <- h; sigma_ab <- c(Sigbeta, Sigalp); sigma_h <- Sigh; t_max = t_pred
        pred_tmp <- sim.tvpVAR.Student.SV(K = K, p = p, t_max = t_pred,
                                          b0 = beta, a0 = alp, h = h,
                                          sigma_ab = c(Sigbeta, Sigalp), sigma_h = Sigh,
                                          y0 = y0, is_tv = is_tv, nu = nu,
                                          seednum = seed_set[count_id], burn_in = 0)
      }

    y_pred[,,count_id] <- pred_tmp$y
    ymean_pred[,,count_id] <- pred_tmp$y_mean
    yvar_pred[,,count_id] <- pred_tmp$y_var
    vol_pred[,,count_id] <- pred_tmp$volatility

  }
  return(list(y_pred = y_pred,
              ymean_pred = ymean_pred,
              yvar_pred = yvar_pred,
              vol_pred = vol_pred))
}

#' Forecast density of future observations of BVAR model
#'
#' This function returns a forecast density of future observations of BVAR-SV-fatTail model.
#' @param Chain The fatBVARSV object from command BVAR.
#' @param y_current The current values of y.
#' @param y_obs_future The future observable values of y.
#' @param t_pred The time prediction horizon.
#' @return The forecast density of future observations of BVAR model
#' @export
#' @examples
#' \dontrun{
#' forecast_dens1 <- forecast_density(Chain1, y_current, y_obs_future, time_current)
#' }
#'
#' @export
forecast_density <- function(Chain, y_current = NULL, y_obs_future, time_current){
  K <- ncol(Chain$data$y)
  p <- Chain$data$p

  if (! ncol(y_obs_future) == K) { stop("ncol(y_obs_future) != K") }
  if (is.null(y_current)) y_current <- Chain$y
  #if (is.null(time_current)) time_current = which(Chain$y == tail(y_current,p)[p,1] )

  t_pred = nrow(y_obs_future)
  Nfsample = 20000
  predictive_samples <- get_forecast(Chain = Chain, y0 = tail(y_current,p),
                                     t_pred = t_pred, Nfsample = Nfsample) # Nfsample

  #############################################################
  # Non cummulative
  #############################################################
  log_pred<- matrix(NA, nrow = t_pred, ncol = K)
  bilog_pred<- matrix(NA, nrow = t_pred, ncol = K*(K-1)*0.5)
  emp_CDF<- matrix(NA, nrow = t_pred, ncol = K)

  for (i in c(1:K)){
    for (j in c(1:t_pred)){
      y_obs <- as.numeric(y_obs_future[j,i])
      log_pred[j,i] <- log(mean(dnorm(y_obs, mean = predictive_samples$ymean_pred[j,i,], sd = sqrt(predictive_samples$vol_pred[j,i,]))))
      emp_CDF[j,i] <- ecdf(x = predictive_samples$y_pred[j,i,])(y_obs)
    }
  }
  # if (K > 1){
  #   Sigma2 <- array(0, dim = c(K, K, Nfsample))
  #   Sigma <- matrix(0, nrow = K, ncol = K)
  #   for (j in c(1:t_pred)){
  #     for (sample_id in c(1:Nfsample)){
  #       Sigma[lower.tri(Sigma, diag = T)] <- predictive_samples$yvar_pred[j,,sample_id]
  #       Sigma2[,,sample_id] <- Sigma %*% t(Sigma)
  #     }
  #
  #     for (i in c(1: (K-1))){
  #       for (k in c((i+1):K)){
  #         pred_bidens <- rep(NA, Nfsample)
  #         for (sample_id in c(1:Nfsample)){
  #           pred_bidens[sample_id] <- dmvn(X = c(y_obs_future[j,i], y_obs_future[j,k]),
  #                                          mu = predictive_samples$ymean_pred[j, c(i,k),sample_id],
  #                                          sigma = Sigma2[c(i,k),c(i,k), sample_id], log = FALSE, isChol = FALSE)
  #         }
  #         bilog_pred[j,(i-1)*K + k-i*(i+1)*0.5] <- log(mean(pred_bidens))
  #       }
  #     }
  #   }
  # }

  #############################################################
  # Cummulative
  #############################################################
  clog_pred<- matrix(NA, nrow = t_pred, ncol = K)
  cbilog_pred<- matrix(NA, nrow = t_pred, ncol = K*(K-1)*0.5)
  cemp_CDF<- matrix(NA, nrow = t_pred, ncol = K)
  cy_pred <- apply(predictive_samples$y_pred, MARGIN = c(2,3), FUN = cumsum)
  cy_obs_future <- apply(y_obs_future, MARGIN = 2, FUN = cumsum)
  # for (i in c(1:K)){
  #   for (j in c(1:t_pred)){
  #     y_obs <- as.numeric(cy_obs_future[j,i])
  #     predict_den <- stats::density(cy_pred[j,i,], adjust = 3)
  #     id <- which(predict_den$y == 0);
  #     if (length(id)>0) {
  #       for (fix in id){
  #         predict_den$y[fix] <- mean(predict_den$y[(fix-1):(fix+1)])
  #       }
  #     }
  #     appximate_density <- smooth.spline(x = predict_den$x, y = log(predict_den$y))
  #     pred_dens <- predict(appximate_density, y_obs)$y
  #     clog_pred[j,i] <- pred_dens
  #     cemp_CDF[j,i] <- ecdf(x = cy_pred[j,i,])(y_obs)
  #   }
  # }

  CRPS <- matrix(NA, nrow = t_pred, ncol = K)
  qwCRPS_2t <- matrix(NA, nrow = t_pred, ncol = K) # Quantile weighted CRPS - 2 tails
  qwCRPS_rt <- matrix(NA, nrow = t_pred, ncol = K) # Quantile weighted CRPS - right tail
  qwCRPS_lt <- matrix(NA, nrow = t_pred, ncol = K) # Quantile weighted CRPS - left tail
  alpha <- seq(0.05, 0.95, by = 0.05)
  for (i in c(1:K)){
    for (j in c(1:t_pred)){
      y_obs <- as.numeric(y_obs_future[j,i])
      samples1 <- predictive_samples$y_pred[j,i,]
      samples2 <- sample(samples1, size = length(samples1), replace = T)
      CRPS[j,i] <- - mean(abs(samples1 - y_obs)) + 0.5*mean(abs(samples1 - samples2))
      Q_tau <- quantile(predictive_samples$y_pred[j,i,], probs = alpha)
      QS <- (y_obs - Q_tau) * ( (y_obs < Q_tau) - alpha )
      qwCRPS_2t[j,i] <- 2/length(alpha) * sum(  (2 * alpha - 1)^2 * QS )
      qwCRPS_rt[j,i] <- 2/length(alpha) * sum(  alpha^2 * QS )
      qwCRPS_lt[j,i] <- 2/length(alpha) * sum(  (1 - alpha)^2 * QS )

    }
  }

  return(list(log_pred = log_pred,
              bilog_pred = bilog_pred,
              emp_CDF = emp_CDF,
              MSFE = (apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = mean) - y_obs_future)^2,
              MAFE = abs(apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = mean) - y_obs_future),
              # predictive_samples = predictive_samples,
              predictive_quantile05 = apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = quantile, probs = 0.05),
              predictive_quantile10 = apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = quantile, probs = 0.10),
              predictive_quantile50 = apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = quantile, probs = 0.50),
              predictive_quantile90 = apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = quantile, probs = 0.90),
              predictive_quantile95 = apply(predictive_samples$y_pred, MARGIN = c(1,2), FUN = quantile, probs = 0.95),
              time_current = time_current,
              clog_pred = clog_pred,
              cbilog_pred = cbilog_pred,
              cemp_CDF = cemp_CDF,
              cMSFE = (apply(cy_pred, MARGIN = c(1,2), FUN = mean) - cy_obs_future)^2,
              cMAFE = abs(apply(cy_pred, MARGIN = c(1,2), FUN = mean) - cy_obs_future),
              CRPS = CRPS,
              qwCRPS_2t = qwCRPS_2t,
              qwCRPS_rt = qwCRPS_rt,
              qwCRPS_lt = qwCRPS_lt
  ))
}

#' @export
recursive_seperate <- function(y_raw, is_tv, t_start = 100, t_pred = 24, K = ncol(y_raw), p = p, dist = "Gaussian", outname = NULL){
  t_max = nrow(y_raw)
  if (is.null(outname)) {
    outname = paste("Recursive_", dist, "_M", paste(is_tv, sep = '', collapse = "") ,
                    "_T", t_start+1, ".RData", sep = "")
  }

  time_current <- t_start # index of the longer y
  y_current <- as.matrix(y_raw[c((p+1):time_current), ], ncol = K)
  y0 <- as.matrix(y_raw[c(1:p), ], ncol = K)
  priors <- get_prior_minnesota(y = y_current, p = p, intercept=TRUE)
  inits <- list(samples = 50000, burnin = 20000, thin = 5)
  inits$is_tv <- is_tv

  if (dist == "Gaussian") {
    Chain <- fitTVPGaussSV(y = y_current, y0, p, priors, inits)
  }
  if (dist == "Student") {
    Chain <- fitTVPStudentSV(y = y_current, y0, p, priors, inits)
  }

  y_obs_future <- as.matrix(y_raw[c((time_current+1):(time_current+t_pred)), ],ncol = K)

  forecast_err <- forecast_density(Chain = Chain, y_current = y_current, y_obs_future = y_obs_future, time_current = time_current)
  out_recursive <- list(time_id = time_current+1,
                        forecast_err = forecast_err,
                        dist = dist,
                        y_current = y_current,
                        y_obs_future = y_obs_future)
  save(out_recursive, file = outname)
  return( out_recursive)

}
