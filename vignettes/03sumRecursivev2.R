################################################################################################
# Trivariate G000
################################################################################################
{
  library(gtools)
  t_pred <- 24
  folderDir <- "/home/hoanguc3m/MEGA/WP11/G000/"
  filenames <- mixedsort(sort(list.files(path = folderDir)))

  time_predict <- seq.Date(from = as.Date('2000-01-01'), to = as.Date('2019-09-01'), by = 'month') # Test
  length_recursive <- length(time_predict)

  G000 <- list(timeID_start = 561,
                     timeID_end = 797,
                     log_pred = array(NA, dim = c(t_pred,3,length_recursive)),
                     emp_CDF = array(NA, dim = c(t_pred,3,length_recursive)),
                     MSFE = array(NA, dim = c(t_pred,3,length_recursive)),
                     MAFE = array(NA, dim = c(t_pred,3,length_recursive)),
                     clog_pred = array(NA, dim = c(t_pred,3,length_recursive)),
                     cemp_CDF = array(NA, dim = c(t_pred,3,length_recursive)),
                     cMSFE = array(NA, dim = c(t_pred,3,length_recursive)),
                     cMAFE = array(NA, dim = c(t_pred,3,length_recursive)),
                     CRPS = array(NA, dim = c(t_pred,3,length_recursive)),
                     qwCRPS_2t = array(NA, dim = c(t_pred,3,length_recursive))
  )

  G011 = T011 = G001 = G111 = T000 = T001 = T111 = G000

  folderDir <- "/home/hoanguc3m/MEGA/WP11/G000/"
  filenames <- mixedsort(sort(list.files(path = folderDir)))

  for (i in c(1:length(filenames))){
    load(paste(folderDir,filenames[i], sep = ""))
    G000$log_pred[,,i] <- out_recursive$forecast_err$log_pred
    G000$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    G000$MSFE[,,i] <- out_recursive$forecast_err$MSFE
    G000$MAFE[,,i] <- out_recursive$forecast_err$MAFE
    G000$clog_pred[,,i] <- out_recursive$forecast_err$clog_pred
    G000$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    G000$cMSFE[,,i] <- out_recursive$forecast_err$cMSFE
    G000$cMAFE[,,i] <- out_recursive$forecast_err$cMAFE
    G000$CRPS[,,i] <- out_recursive$forecast_err$CRPS
    G000$qwCRPS_2t[,,i] <- out_recursive$forecast_err$qwCRPS_2t
  }

  G000$mlog_pred <- apply(G000$log_pred, MARGIN = c(1,2), FUN = mean)
  G000$mMSFE <- apply(G000$MSFE, MARGIN = c(1,2), FUN = mean)
  G000$mMAFE <- apply(G000$MAFE, MARGIN = c(1,2), FUN = mean)
  G000$mclog_pred <- apply(G000$clog_pred, MARGIN = c(1,2), FUN = mean)
  G000$mcMSFE <- apply(G000$cMSFE, MARGIN = c(1,2), FUN = mean)
  G000$mcMAFE <- apply(G000$cMAFE, MARGIN = c(1,2), FUN = mean)
  G000$mCRPS <- apply(G000$CRPS, MARGIN = c(1,2), FUN = mean)
  G000$mqwCRPS_2t <- apply(G000$qwCRPS_2t, MARGIN = c(1,2), FUN = mean)



  ##########################################
  folderDir <- "/home/hoanguc3m/MEGA/WP11/G001/"
  filenames <- mixedsort(sort(list.files(path = folderDir)))

  for (i in c(1:length(filenames))){
    load(paste(folderDir,filenames[i], sep = ""))
    G001$log_pred[,,i] <- out_recursive$forecast_err$log_pred
    G001$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    G001$MSFE[,,i] <- out_recursive$forecast_err$MSFE
    G001$MAFE[,,i] <- out_recursive$forecast_err$MAFE
    G001$clog_pred[,,i] <- out_recursive$forecast_err$clog_pred
    G001$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    G001$cMSFE[,,i] <- out_recursive$forecast_err$cMSFE
    G001$cMAFE[,,i] <- out_recursive$forecast_err$cMAFE
    G001$CRPS[,,i] <- out_recursive$forecast_err$CRPS
    G001$qwCRPS_2t[,,i] <- out_recursive$forecast_err$qwCRPS_2t
  }

  G001$mlog_pred <- apply(G001$log_pred, MARGIN = c(1,2), FUN = mean)
  G001$mMSFE <- apply(G001$MSFE, MARGIN = c(1,2), FUN = mean)
  G001$mMAFE <- apply(G001$MAFE, MARGIN = c(1,2), FUN = mean)
  G001$mclog_pred <- apply(G001$clog_pred, MARGIN = c(1,2), FUN = mean)
  G001$mcMSFE <- apply(G001$cMSFE, MARGIN = c(1,2), FUN = mean)
  G001$mcMAFE <- apply(G001$cMAFE, MARGIN = c(1,2), FUN = mean)
  G001$mCRPS <- apply(G001$CRPS, MARGIN = c(1,2), FUN = mean)
  G001$mqwCRPS_2t <- apply(G001$qwCRPS_2t, MARGIN = c(1,2), FUN = mean)

  ##########################################
  ##########################################
  folderDir <- "/home/hoanguc3m/MEGA/WP11/G011/"
  filenames <- mixedsort(sort(list.files(path = folderDir)))

  for (i in c(1:length(filenames))){
    load(paste(folderDir,filenames[i], sep = ""))
    G011$log_pred[,,i] <- out_recursive$forecast_err$log_pred
    G011$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    G011$MSFE[,,i] <- out_recursive$forecast_err$MSFE
    G011$MAFE[,,i] <- out_recursive$forecast_err$MAFE
    G011$clog_pred[,,i] <- out_recursive$forecast_err$clog_pred
    G011$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    G011$cMSFE[,,i] <- out_recursive$forecast_err$cMSFE
    G011$cMAFE[,,i] <- out_recursive$forecast_err$cMAFE
    G011$CRPS[,,i] <- out_recursive$forecast_err$CRPS
    G011$qwCRPS_2t[,,i] <- out_recursive$forecast_err$qwCRPS_2t
  }

  G011$mlog_pred <- apply(G011$log_pred, MARGIN = c(1,2), FUN = mean)
  G011$mMSFE <- apply(G011$MSFE, MARGIN = c(1,2), FUN = mean)
  G011$mMAFE <- apply(G011$MAFE, MARGIN = c(1,2), FUN = mean)
  G011$mclog_pred <- apply(G011$clog_pred, MARGIN = c(1,2), FUN = mean)
  G011$mcMSFE <- apply(G011$cMSFE, MARGIN = c(1,2), FUN = mean)
  G011$mcMAFE <- apply(G011$cMAFE, MARGIN = c(1,2), FUN = mean)
  G011$mCRPS <- apply(G011$CRPS, MARGIN = c(1,2), FUN = mean)
  G011$mqwCRPS_2t <- apply(G011$qwCRPS_2t, MARGIN = c(1,2), FUN = mean)

  ############################################

  folderDir <- "/home/hoanguc3m/MEGA/WP11/G111/"
  filenames <- mixedsort(sort(list.files(path = folderDir)))

  for (i in c(1:length(filenames))){
    load(paste(folderDir,filenames[i], sep = ""))
    G111$log_pred[,,i] <- out_recursive$forecast_err$log_pred
    G111$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    G111$MSFE[,,i] <- out_recursive$forecast_err$MSFE
    G111$MAFE[,,i] <- out_recursive$forecast_err$MAFE
    G111$clog_pred[,,i] <- out_recursive$forecast_err$clog_pred
    G111$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    G111$cMSFE[,,i] <- out_recursive$forecast_err$cMSFE
    G111$cMAFE[,,i] <- out_recursive$forecast_err$cMAFE
    G111$CRPS[,,i] <- out_recursive$forecast_err$CRPS
    G111$qwCRPS_2t[,,i] <- out_recursive$forecast_err$qwCRPS_2t
  }

  G111$mlog_pred <- apply(G111$log_pred, MARGIN = c(1,2), FUN = mean)
  G111$mMSFE <- apply(G111$MSFE, MARGIN = c(1,2), FUN = mean)
  G111$mMAFE <- apply(G111$MAFE, MARGIN = c(1,2), FUN = mean)
  G111$mclog_pred <- apply(G111$clog_pred, MARGIN = c(1,2), FUN = mean)
  G111$mcMSFE <- apply(G111$cMSFE, MARGIN = c(1,2), FUN = mean)
  G111$mcMAFE <- apply(G111$cMAFE, MARGIN = c(1,2), FUN = mean)
  G111$mCRPS <- apply(G111$CRPS, MARGIN = c(1,2), FUN = mean)
  G111$mqwCRPS_2t <- apply(G111$qwCRPS_2t, MARGIN = c(1,2), FUN = mean)

  ######################################
  folderDir <- "/home/hoanguc3m/MEGA/WP11/T000/"
  filenames <- mixedsort(sort(list.files(path = folderDir)))

  for (i in c(1:length(filenames))){
    load(paste(folderDir,filenames[i], sep = ""))
    T000$log_pred[,,i] <- out_recursive$forecast_err$log_pred
    T000$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    T000$MSFE[,,i] <- out_recursive$forecast_err$MSFE
    T000$MAFE[,,i] <- out_recursive$forecast_err$MAFE
    T000$clog_pred[,,i] <- out_recursive$forecast_err$clog_pred
    T000$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    T000$cMSFE[,,i] <- out_recursive$forecast_err$cMSFE
    T000$cMAFE[,,i] <- out_recursive$forecast_err$cMAFE
    T000$CRPS[,,i] <- out_recursive$forecast_err$CRPS
    T000$qwCRPS_2t[,,i] <- out_recursive$forecast_err$qwCRPS_2t
  }

  T000$mlog_pred <- apply(T000$log_pred, MARGIN = c(1,2), FUN = mean)
  T000$mMSFE <- apply(T000$MSFE, MARGIN = c(1,2), FUN = mean)
  T000$mMAFE <- apply(T000$MAFE, MARGIN = c(1,2), FUN = mean)
  T000$mclog_pred <- apply(T000$clog_pred, MARGIN = c(1,2), FUN = mean)
  T000$mcMSFE <- apply(T000$cMSFE, MARGIN = c(1,2), FUN = mean)
  T000$mcMAFE <- apply(T000$cMAFE, MARGIN = c(1,2), FUN = mean)
  T000$mCRPS <- apply(T000$CRPS, MARGIN = c(1,2), FUN = mean)
  T000$mqwCRPS_2t <- apply(T000$qwCRPS_2t, MARGIN = c(1,2), FUN = mean)

  #########################################
  folderDir <- "/home/hoanguc3m/MEGA/WP11/T001/"
  filenames <- mixedsort(sort(list.files(path = folderDir)))

  for (i in c(1:length(filenames))){
    load(paste(folderDir,filenames[i], sep = ""))
    T001$log_pred[,,i] <- out_recursive$forecast_err$log_pred
    T001$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    T001$MSFE[,,i] <- out_recursive$forecast_err$MSFE
    T001$MAFE[,,i] <- out_recursive$forecast_err$MAFE
    T001$clog_pred[,,i] <- out_recursive$forecast_err$clog_pred
    T001$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    T001$cMSFE[,,i] <- out_recursive$forecast_err$cMSFE
    T001$cMAFE[,,i] <- out_recursive$forecast_err$cMAFE
    T001$CRPS[,,i] <- out_recursive$forecast_err$CRPS
    T001$qwCRPS_2t[,,i] <- out_recursive$forecast_err$qwCRPS_2t
  }

  T001$mlog_pred <- apply(T001$log_pred, MARGIN = c(1,2), FUN = mean)
  T001$mMSFE <- apply(T001$MSFE, MARGIN = c(1,2), FUN = mean)
  T001$mMAFE <- apply(T001$MAFE, MARGIN = c(1,2), FUN = mean)
  T001$mclog_pred <- apply(T001$clog_pred, MARGIN = c(1,2), FUN = mean)
  T001$mcMSFE <- apply(T001$cMSFE, MARGIN = c(1,2), FUN = mean)
  T001$mcMAFE <- apply(T001$cMAFE, MARGIN = c(1,2), FUN = mean)
  T001$mCRPS <- apply(T001$CRPS, MARGIN = c(1,2), FUN = mean)
  T001$mqwCRPS_2t <- apply(T001$qwCRPS_2t, MARGIN = c(1,2), FUN = mean)

  ##########################################
  folderDir <- "/home/hoanguc3m/MEGA/WP11/T011/"
  filenames <- mixedsort(sort(list.files(path = folderDir)))

  for (i in c(1:length(filenames))){
    load(paste(folderDir,filenames[i], sep = ""))
    T011$log_pred[,,i] <- out_recursive$forecast_err$log_pred
    T011$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    T011$MSFE[,,i] <- out_recursive$forecast_err$MSFE
    T011$MAFE[,,i] <- out_recursive$forecast_err$MAFE
    T011$clog_pred[,,i] <- out_recursive$forecast_err$clog_pred
    T011$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    T011$cMSFE[,,i] <- out_recursive$forecast_err$cMSFE
    T011$cMAFE[,,i] <- out_recursive$forecast_err$cMAFE
    T011$CRPS[,,i] <- out_recursive$forecast_err$CRPS
    T011$qwCRPS_2t[,,i] <- out_recursive$forecast_err$qwCRPS_2t
  }

  T011$mlog_pred <- apply(T011$log_pred, MARGIN = c(1,2), FUN = mean)
  T011$mMSFE <- apply(T011$MSFE, MARGIN = c(1,2), FUN = mean)
  T011$mMAFE <- apply(T011$MAFE, MARGIN = c(1,2), FUN = mean)
  T011$mclog_pred <- apply(T011$clog_pred, MARGIN = c(1,2), FUN = mean)
  T011$mcMSFE <- apply(T011$cMSFE, MARGIN = c(1,2), FUN = mean)
  T011$mcMAFE <- apply(T011$cMAFE, MARGIN = c(1,2), FUN = mean)
  T011$mCRPS <- apply(T011$CRPS, MARGIN = c(1,2), FUN = mean)
  T011$mqwCRPS_2t <- apply(T011$qwCRPS_2t, MARGIN = c(1,2), FUN = mean)

  ###############################################
  folderDir <- "/home/hoanguc3m/MEGA/WP11/T111/"
  filenames <- mixedsort(sort(list.files(path = folderDir)))

  for (i in c(1:length(filenames))){
    load(paste(folderDir,filenames[i], sep = ""))
    T111$log_pred[,,i] <- out_recursive$forecast_err$log_pred
    T111$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    T111$MSFE[,,i] <- out_recursive$forecast_err$MSFE
    T111$MAFE[,,i] <- out_recursive$forecast_err$MAFE
    T111$clog_pred[,,i] <- out_recursive$forecast_err$clog_pred
    T111$emp_CDF[,,i] <- out_recursive$forecast_err$emp_CDF
    T111$cMSFE[,,i] <- out_recursive$forecast_err$cMSFE
    T111$cMAFE[,,i] <- out_recursive$forecast_err$cMAFE
    T111$CRPS[,,i] <- out_recursive$forecast_err$CRPS
    T111$qwCRPS_2t[,,i] <- out_recursive$forecast_err$qwCRPS_2t
  }

  T111$mlog_pred <- apply(T111$log_pred, MARGIN = c(1,2), FUN = mean)
  T111$mMSFE <- apply(T111$MSFE, MARGIN = c(1,2), FUN = mean)
  T111$mMAFE <- apply(T111$MAFE, MARGIN = c(1,2), FUN = mean)
  T111$mclog_pred <- apply(T111$clog_pred, MARGIN = c(1,2), FUN = mean)
  T111$mcMSFE <- apply(T111$cMSFE, MARGIN = c(1,2), FUN = mean)
  T111$mcMAFE <- apply(T111$cMAFE, MARGIN = c(1,2), FUN = mean)
  T111$mCRPS <- apply(T111$CRPS, MARGIN = c(1,2), FUN = mean)
  T111$mqwCRPS_2t <- apply(T111$qwCRPS_2t, MARGIN = c(1,2), FUN = mean)

}


################################################################################################
# Export tables
################################################################################################
{
library(sandwich)
Diebold.test <- function (e1, e2, alternative = c("two.sided", "less", "greater"), h = 1, power = 1) {
  alternative <- match.arg(alternative)
  d <- c(e1)^power - c(e2)^power
  d.cov <- acf(d, na.action = na.omit, lag.max = h - 1, type = "covariance",
               plot = FALSE)$acf[, , 1]
  d.var <- sum(c(d.cov[1], 2 * d.cov[-1]))/length(d)
  dv <- d.var
  if (dv > 0) {
    STATISTIC <- mean(d, na.rm = TRUE)/sqrt(dv)
  } else if (h == 1) {
    stop("Variance of DM statistic is zero")
  } else {
    warning("Variance is negative, using horizon h=1")
    return(dm.test(e1, e2, alternative, h = 1, power))
  }
  n <- length(d)
  k <- ((n + 1 - 2 * h + (h/n) * (h - 1))/n)^(1/2)
  STATISTIC <- STATISTIC * k
  names(STATISTIC) <- "DM"
  if (alternative == "two.sided") {
    PVAL <- 2 * pt(-abs(STATISTIC), df = n - 1)
  } else if (alternative == "less") {
    PVAL <- pt(STATISTIC, df = n - 1)
  } else if (alternative == "greater") {
    PVAL <- pt(STATISTIC, df = n - 1, lower.tail = FALSE)
  }
  PARAMETER <- c(h, power)
  names(PARAMETER) <- c("Forecast horizon", "Loss function power")
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 alternative = alternative, p.value = PVAL, method = "Diebold-Mariano Test",
                 data.name = c(deparse(substitute(e1)), deparse(substitute(e2)))),
            class = "htest")
}

ClackMcCracken.test <- function(e1, e2, alternative = c("two.sided", "less", "greater"), h = 1, power = 1){
  alternative <- match.arg(alternative)
  #     #d <- c(abs(e1))^power - c(abs(e2))^power
  # d <- c(e1)^power - c(e2)^power
  # d.cov <- acf(d, na.action = na.omit, lag.max = floor(1.5*h) , type = "covariance",
  #              plot = FALSE)$acf[, , 1] # NeweyWest max lag = 1.5 h # DM lag.max = h -1
  # m = length(d.cov)
  # d.var <- (d.cov[1] + 2 * sum( (1 - (1:(m-1)) /m) * d.cov[-1] ) ) / length(d)
  #     #d.var <- sum(c(d.cov[1], 2 * d.cov[-1]))/length(d)    # This is DM var

  d <- c(e1)^power - c(e2)^power

  fm <- lm(d ~ 1)
  d.var <- NeweyWest(fm, lag = floor(1.5*h), prewhite = FALSE, adjust = FALSE)
  #d.var <- NeweyWest(fm, lag = NULL, prewhite = FALSE, adjust = FALSE)

  dv <- d.var
  if (dv > 0) {
    STATISTIC <- mean(d, na.rm = TRUE)/sqrt(dv)
  } else if (h == 1) {
    stop("Variance of DM statistic is zero")
  } else {
    warning("Variance is negative, using horizon h=1")
    return(dm.test(e1, e2, alternative, h = 1, power))
  }
  n <- length(d)
  k <- ((n + 1 - 2 * h + (h/n) * (h - 1))/n)^(1/2)
  STATISTIC <- STATISTIC * k
  names(STATISTIC) <- "DM"
  if (alternative == "two.sided") {
    PVAL <- 2 * pt(-abs(STATISTIC), df = n - 1)
  }
  else if (alternative == "less") {
    PVAL <- pt(STATISTIC, df = n - 1)
  }
  else if (alternative == "greater") {
    PVAL <- pt(STATISTIC, df = n - 1, lower.tail = FALSE)
  }

  PARAMETER <- c(h, power)
  names(PARAMETER) <- c("Forecast horizon", "Loss function power")
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 alternative = alternative, p.value = PVAL, method = "DMW - Clark Test",
                 data.name = c(deparse(substitute(e1)), deparse(substitute(e2)))),
            class = "htest")
}

DM_Andrews.test <- function(e1, e2, alternative = c("two.sided", "less", "greater"), h = 1, power = 1){
  alternative <- match.arg(alternative)
  #d <- c(abs(e1))^power - c(abs(e2))^power
  d <- c(e1)^power - c(e2)^power

  fm <- lm(d ~ 1)
  #d.var <- kernHAC(fm, prewhite = 1, bw = bwAndrews, kernel = c("Quadratic Spectral"), approx = c("ARMA(1,1)"))
  d.var <- kernHAC(fm, prewhite = 1, bw = bwAndrews, kernel = c("Quadratic Spectral"), approx = c("AR(1)"))

  dv <- d.var
  if (dv > 0) {
    STATISTIC <- mean(d, na.rm = TRUE)/sqrt(dv)
  } else if (h == 1) {
    stop("Variance of DM statistic is zero")
  } else {
    warning("Variance is negative, using horizon h=1")
    return(dm.test(e1, e2, alternative, h = 1, power))
  }
  n <- length(d)
  k <- ((n + 1 - 2 * h + (h/n) * (h - 1))/n)^(1/2)
  STATISTIC <- STATISTIC * k
  names(STATISTIC) <- "DM"
  if (alternative == "two.sided") {
    PVAL <- 2 * pt(-abs(STATISTIC), df = n - 1)
  }
  else if (alternative == "less") {
    PVAL <- pt(STATISTIC, df = n - 1)
  }
  else if (alternative == "greater") {
    PVAL <- pt(STATISTIC, df = n - 1, lower.tail = FALSE)
  }

  PARAMETER <- c(h, power)
  names(PARAMETER) <- c("Forecast horizon", "Loss function power")
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 alternative = alternative, p.value = PVAL, method = "DM Andrews Test",
                 data.name = c(deparse(substitute(e1)), deparse(substitute(e2)))),
            class = "htest")
}

Print_DmTableMSE <- function(G000MSFE,G001MSFE,G011MSFE, G111MSFE,
                             T000MSFE,T001MSFE,T011MSFE, T111MSFE,
                             Table3_MSFE_Var1){
  Table3_MSFE_Var1_add <- Table3_MSFE_Var1
  id <- c(1,2,3,6,12,24)
  for (i in c(1:length(id))){
    tmpMat <- cbind(G000MSFE[id[i],],G001MSFE[id[i],],G011MSFE[id[i],],G111MSFE[id[i],],
                    T000MSFE[id[i],],T001MSFE[id[i],],T011MSFE[id[i],],T111MSFE[id[i],])
    for (j in c(2:8)){
      #tmpDmtest <- Diebold.test(tmpMat[,1] , tmpMat[,j], h=id[i], power = 1, alternative = "greater")
      #tmpDmtest <- DM_Andrews.test(tmpMat[,1] , tmpMat[,j], h=id[i], power = 1, alternative = "greater")
      tmpDmtest <- ClackMcCracken.test(tmpMat[,1] , tmpMat[,j], h=id[i], power = 1, alternative = "greater")
      greater <- (mean(tmpMat[,1]) - mean(tmpMat[,j])) > 0
      if (tmpDmtest$p.value < 0.10 && greater) Table3_MSFE_Var1_add[j,i] <- paste(Table3_MSFE_Var1[j,i], "*", sep = "")
      if (tmpDmtest$p.value < 0.05 && greater) Table3_MSFE_Var1_add[j,i] <- paste(Table3_MSFE_Var1[j,i], "**", sep = "")
      if (tmpDmtest$p.value < 0.01 && greater) Table3_MSFE_Var1_add[j,i] <- paste(Table3_MSFE_Var1[j,i], "***", sep = "")
    }

  }
  print(xtable::xtable(Table3_MSFE_Var1_add, digits = 3), include.rownames=T, type = "latex", sanitize.text.function = function(x){x})

}

Print_DmTableDens <- function(G000MSFE,G001MSFE,G011MSFE, G111MSFE,
                              T000MSFE,T001MSFE,T011MSFE, T111MSFE,
                              Table3_MSFE_Var1){
  Table3_MSFE_Var1_add <- Table3_MSFE_Var1
  id <- c(1,2,3,6,12,24)
  for (i in c(1:length(id))){
    tmpMat <- cbind(G000MSFE[id[i],],G001MSFE[id[i],],G011MSFE[id[i],],G111MSFE[id[i],],
                    T000MSFE[id[i],],T001MSFE[id[i],],T011MSFE[id[i],],T111MSFE[id[i],])
    for (j in c(2:8)){
      #tmpDmtest <- Diebold.test(tmpMat[,1] , tmpMat[,j], h=id[i], power = 1, alternative = "less")
      #tmpDmtest <- DM_Andrews.test(tmpMat[,1] , tmpMat[,j], h=id[i], power = 1, alternative = "less")
      tmpDmtest <- ClackMcCracken.test(tmpMat[,1] , tmpMat[,j], h=id[i], power = 1, alternative = "less")
      less <- (mean(tmpMat[,1]) - mean(tmpMat[,j])) < 0
      if (tmpDmtest$p.value < 0.10 && less) Table3_MSFE_Var1_add[j,i] <- paste(Table3_MSFE_Var1[j,i], "*", sep = "")
      if (tmpDmtest$p.value < 0.05 && less) Table3_MSFE_Var1_add[j,i] <- paste(Table3_MSFE_Var1[j,i], "**", sep = "")
      if (tmpDmtest$p.value < 0.01 && less) Table3_MSFE_Var1_add[j,i] <- paste(Table3_MSFE_Var1[j,i], "***", sep = "")
    }

  }
  print(xtable::xtable(Table3_MSFE_Var1_add, digits = 3), include.rownames=T, type = "latex", sanitize.text.function = function(x){x})

}
}
# MSE
{
  Table3_MSFE_Var1 <- matrix(NA, ncol = 6, nrow = 8)
  Table3_MSFE_Var1[1,] <- G000$mMSFE[c(1,2,3,6,12,24),1]
  Table3_MSFE_Var1[2,] <- G001$mMSFE[c(1,2,3,6,12,24),1] / G000$mMSFE[c(1,2,3,6,12,24),1]
  Table3_MSFE_Var1[3,] <- G011$mMSFE[c(1,2,3,6,12,24),1] / G000$mMSFE[c(1,2,3,6,12,24),1]
  Table3_MSFE_Var1[4,] <- G111$mMSFE[c(1,2,3,6,12,24),1] / G000$mMSFE[c(1,2,3,6,12,24),1]
  Table3_MSFE_Var1[5,] <- T000$mMSFE[c(1,2,3,6,12,24),1] / G000$mMSFE[c(1,2,3,6,12,24),1]
  Table3_MSFE_Var1[6,] <- T001$mMSFE[c(1,2,3,6,12,24),1] / G000$mMSFE[c(1,2,3,6,12,24),1]
  Table3_MSFE_Var1[7,] <- T011$mMSFE[c(1,2,3,6,12,24),1] / G000$mMSFE[c(1,2,3,6,12,24),1]
  Table3_MSFE_Var1[8,] <- T111$mMSFE[c(1,2,3,6,12,24),1] / G000$mMSFE[c(1,2,3,6,12,24),1]

  Table3_MSFE_Var1 <- matrix(sprintf("%.3f", Table3_MSFE_Var1), ncol = 6, nrow = 8)
  rownames(Table3_MSFE_Var1) <- c("G000", "G001", "G011", "G111", "T000", "T001", "T011", "T111")
  colnames(Table3_MSFE_Var1) <- c("1", "2", "3", "6", "12", "24")


  Print_DmTableMSE(G000$MSFE[,1,],
                   G001$MSFE[,1,],
                   G011$MSFE[,1,],
                   G111$MSFE[,1,],
                   T000$MSFE[,1,],
                   T001$MSFE[,1,],
                   T011$MSFE[,1,],
                   T111$MSFE[,1,], Table3_MSFE_Var1)

  Table3_MSFE_Var2 <- matrix(NA, ncol = 6, nrow = 8)
  Table3_MSFE_Var2[1,] <- G000$mMSFE[c(1,2,3,6,12,24),2]
  Table3_MSFE_Var2[2,] <- G001$mMSFE[c(1,2,3,6,12,24),2] / G000$mMSFE[c(1,2,3,6,12,24),2]
  Table3_MSFE_Var2[3,] <- G011$mMSFE[c(1,2,3,6,12,24),2] / G000$mMSFE[c(1,2,3,6,12,24),2]
  Table3_MSFE_Var2[4,] <- G111$mMSFE[c(1,2,3,6,12,24),2] / G000$mMSFE[c(1,2,3,6,12,24),2]
  Table3_MSFE_Var2[5,] <- T000$mMSFE[c(1,2,3,6,12,24),2] / G000$mMSFE[c(1,2,3,6,12,24),2]
  Table3_MSFE_Var2[6,] <- T001$mMSFE[c(1,2,3,6,12,24),2] / G000$mMSFE[c(1,2,3,6,12,24),2]
  Table3_MSFE_Var2[7,] <- T011$mMSFE[c(1,2,3,6,12,24),2] / G000$mMSFE[c(1,2,3,6,12,24),2]
  Table3_MSFE_Var2[8,] <- T111$mMSFE[c(1,2,3,6,12,24),2] / G000$mMSFE[c(1,2,3,6,12,24),2]

  Table3_MSFE_Var2 <- matrix(sprintf("%.3f", Table3_MSFE_Var2), ncol = 6, nrow = 8)
  rownames(Table3_MSFE_Var2) <- c("G000", "G001", "G011", "G111", "T000", "T001", "T011", "T111")
  colnames(Table3_MSFE_Var2) <- c("1", "2", "3", "6", "12", "24")


  Print_DmTableMSE(G000$MSFE[,2,],
                   G001$MSFE[,2,],
                   G011$MSFE[,2,],
                   G111$MSFE[,2,],
                   T000$MSFE[,2,],
                   T001$MSFE[,2,],
                   T011$MSFE[,2,],
                   T111$MSFE[,2,], Table3_MSFE_Var2)

  Table3_MSFE_Var3 <- matrix(NA, ncol = 6, nrow = 8)
  Table3_MSFE_Var3[1,] <- G000$mMSFE[c(1,2,3,6,12,24),3]
  Table3_MSFE_Var3[2,] <- G001$mMSFE[c(1,2,3,6,12,24),3] / G000$mMSFE[c(1,2,3,6,12,24),3]
  Table3_MSFE_Var3[3,] <- G011$mMSFE[c(1,2,3,6,12,24),3] / G000$mMSFE[c(1,2,3,6,12,24),3]
  Table3_MSFE_Var3[4,] <- G111$mMSFE[c(1,2,3,6,12,24),3] / G000$mMSFE[c(1,2,3,6,12,24),3]
  Table3_MSFE_Var3[5,] <- T000$mMSFE[c(1,2,3,6,12,24),3] / G000$mMSFE[c(1,2,3,6,12,24),3]
  Table3_MSFE_Var3[6,] <- T001$mMSFE[c(1,2,3,6,12,24),3] / G000$mMSFE[c(1,2,3,6,12,24),3]
  Table3_MSFE_Var3[7,] <- T011$mMSFE[c(1,2,3,6,12,24),3] / G000$mMSFE[c(1,2,3,6,12,24),3]
  Table3_MSFE_Var3[8,] <- T111$mMSFE[c(1,2,3,6,12,24),3] / G000$mMSFE[c(1,2,3,6,12,24),3]
  Table3_MSFE_Var3 <- matrix(sprintf("%.3f", Table3_MSFE_Var3), ncol = 6, nrow = 8)
  rownames(Table3_MSFE_Var3) <- c("G000", "G001", "G011", "G111", "T000", "T001", "T011", "T111")
  colnames(Table3_MSFE_Var3) <- c("1", "2", "3", "6", "12", "24")


  Print_DmTableMSE(G000$MSFE[,3,],
                   G001$MSFE[,3,],
                   G011$MSFE[,3,],
                   G111$MSFE[,3,],
                   T000$MSFE[,3,],
                   T001$MSFE[,3,],
                   T011$MSFE[,3,],
                   T111$MSFE[,3,], Table3_MSFE_Var3)
}

# LP
{
  Table4_LP_Var1 <-  matrix(NA, ncol = 6, nrow = 8)

  Table4_LP_Var1[1,] <- G000$mlog_pred[c(1,2,3,6,12,24),1]
  Table4_LP_Var1[2,] <- G001$mlog_pred[c(1,2,3,6,12,24),1] - G000$mlog_pred[c(1,2,3,6,12,24),1]
  Table4_LP_Var1[3,] <- G011$mlog_pred[c(1,2,3,6,12,24),1] - G000$mlog_pred[c(1,2,3,6,12,24),1]
  Table4_LP_Var1[4,] <- G111$mlog_pred[c(1,2,3,6,12,24),1] - G000$mlog_pred[c(1,2,3,6,12,24),1]
  Table4_LP_Var1[5,] <- T000$mlog_pred[c(1,2,3,6,12,24),1] - G000$mlog_pred[c(1,2,3,6,12,24),1]
  Table4_LP_Var1[6,] <- T001$mlog_pred[c(1,2,3,6,12,24),1] - G000$mlog_pred[c(1,2,3,6,12,24),1]
  Table4_LP_Var1[7,] <- T011$mlog_pred[c(1,2,3,6,12,24),1] - G000$mlog_pred[c(1,2,3,6,12,24),1]
  Table4_LP_Var1[8,] <- T111$mlog_pred[c(1,2,3,6,12,24),1] - G000$mlog_pred[c(1,2,3,6,12,24),1]

  Table4_LP_Var1 <- matrix(sprintf("%.3f", Table4_LP_Var1), ncol = 6, nrow = 8)
  row.names(Table4_LP_Var1) <- c("G000", "G001", "G011", "G111", "T000", "T001", "T011", "T111")

  Print_DmTableDens(G000$log_pred[,1,],
                    G001$log_pred[,1,],
                    G011$log_pred[,1,],
                    G111$log_pred[,1,],
                    T000$log_pred[,1,],
                    T001$log_pred[,1,],
                    T011$log_pred[,1,],
                    T111$log_pred[,1,],
                    Table4_LP_Var1)

  Table4_LP_Var2 <-  matrix(NA, ncol = 6, nrow = 8)

  Table4_LP_Var2[1,] <- G000$mlog_pred[c(1,2,3,6,12,24),2]
  Table4_LP_Var2[2,] <- G001$mlog_pred[c(1,2,3,6,12,24),2] - G000$mlog_pred[c(1,2,3,6,12,24),2]
  Table4_LP_Var2[3,] <- G011$mlog_pred[c(1,2,3,6,12,24),2] - G000$mlog_pred[c(1,2,3,6,12,24),2]
  Table4_LP_Var2[4,] <- G111$mlog_pred[c(1,2,3,6,12,24),2] - G000$mlog_pred[c(1,2,3,6,12,24),2]
  Table4_LP_Var2[5,] <- T000$mlog_pred[c(1,2,3,6,12,24),2] - G000$mlog_pred[c(1,2,3,6,12,24),2]
  Table4_LP_Var2[6,] <- T001$mlog_pred[c(1,2,3,6,12,24),2] - G000$mlog_pred[c(1,2,3,6,12,24),2]
  Table4_LP_Var2[7,] <- T011$mlog_pred[c(1,2,3,6,12,24),2] - G000$mlog_pred[c(1,2,3,6,12,24),2]
  Table4_LP_Var2[8,] <- T111$mlog_pred[c(1,2,3,6,12,24),2] - G000$mlog_pred[c(1,2,3,6,12,24),2]

  Table4_LP_Var2 <- matrix(sprintf("%.3f", Table4_LP_Var2), ncol = 6, nrow = 8)
  row.names(Table4_LP_Var2) <- c("G000", "G001", "G011", "G111", "T000", "T001", "T011", "T111")

  Print_DmTableDens(G000$log_pred[,2,],
                    G001$log_pred[,2,],
                    G011$log_pred[,2,],
                    G111$log_pred[,2,],
                    T000$log_pred[,2,],
                    T001$log_pred[,2,],
                    T011$log_pred[,2,],
                    T111$log_pred[,2,],
                    Table4_LP_Var2)

  Table4_LP_Var3 <-  matrix(NA, ncol = 6, nrow = 8)

  Table4_LP_Var3[1,] <- G000$mlog_pred[c(1,2,3,6,12,24),3]
  Table4_LP_Var3[2,] <- G001$mlog_pred[c(1,2,3,6,12,24),3] - G000$mlog_pred[c(1,2,3,6,12,24),3]
  Table4_LP_Var3[3,] <- G011$mlog_pred[c(1,2,3,6,12,24),3] - G000$mlog_pred[c(1,2,3,6,12,24),3]
  Table4_LP_Var3[4,] <- G111$mlog_pred[c(1,2,3,6,12,24),3] - G000$mlog_pred[c(1,2,3,6,12,24),3]
  Table4_LP_Var3[5,] <- T000$mlog_pred[c(1,2,3,6,12,24),3] - G000$mlog_pred[c(1,2,3,6,12,24),3]
  Table4_LP_Var3[6,] <- T001$mlog_pred[c(1,2,3,6,12,24),3] - G000$mlog_pred[c(1,2,3,6,12,24),3]
  Table4_LP_Var3[7,] <- T011$mlog_pred[c(1,2,3,6,12,24),3] - G000$mlog_pred[c(1,2,3,6,12,24),3]
  Table4_LP_Var3[8,] <- T111$mlog_pred[c(1,2,3,6,12,24),3] - G000$mlog_pred[c(1,2,3,6,12,24),3]

  Table4_LP_Var3 <- matrix(sprintf("%.3f", Table4_LP_Var3), ncol = 6, nrow = 8)
  row.names(Table4_LP_Var3) <- c("G000", "G001", "G011", "G111", "T000", "T001", "T011", "T111")

  Print_DmTableDens(G000$log_pred[,3,],
                    G001$log_pred[,3,],
                    G011$log_pred[,3,],
                    G111$log_pred[,3,],
                    T000$log_pred[,3,],
                    T001$log_pred[,3,],
                    T011$log_pred[,3,],
                    T111$log_pred[,3,],
                    Table4_LP_Var3)
}

# CRPS
{
  Table5_qwCRPS_2t_Var1 <-  matrix(NA, ncol = 6, nrow = 8)

  Table5_qwCRPS_2t_Var1[1,] <- G000$mqwCRPS_2t[c(1,2,3,6,12,24),1]
  Table5_qwCRPS_2t_Var1[2,] <- G001$mqwCRPS_2t[c(1,2,3,6,12,24),1] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),1]
  Table5_qwCRPS_2t_Var1[3,] <- G111$mqwCRPS_2t[c(1,2,3,6,12,24),1] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),1]
  Table5_qwCRPS_2t_Var1[4,] <- T000$mqwCRPS_2t[c(1,2,3,6,12,24),1] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),1]
  Table5_qwCRPS_2t_Var1[5,] <- T001$mqwCRPS_2t[c(1,2,3,6,12,24),1] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),1]
  Table5_qwCRPS_2t_Var1[6,] <- T111$mqwCRPS_2t[c(1,2,3,6,12,24),1] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),1]

  Table5_qwCRPS_2t_Var1 <- matrix(sprintf("%.3f", Table5_qwCRPS_2t_Var1), ncol = 6, nrow = 8)
  row.names(Table5_qwCRPS_2t_Var1) <- c("G000", "G001", "G011", "G111", "T000", "T001", "T011", "T111")

  Print_DmTableDens(G000$qwCRPS_2t[,1,],
                    G001$qwCRPS_2t[,1,],
                    G111$qwCRPS_2t[,1,],
                    T000$qwCRPS_2t[,1,],
                    T001$qwCRPS_2t[,1,],
                    T111$qwCRPS_2t[,1,],
                    Table5_qwCRPS_2t_Var1)

  Table5_qwCRPS_2t_Var2 <-  matrix(NA, ncol = 6, nrow = 8)

  Table5_qwCRPS_2t_Var2[1,] <- G000$mqwCRPS_2t[c(1,2,3,6,12,24),2]
  Table5_qwCRPS_2t_Var2[2,] <- G001$mqwCRPS_2t[c(1,2,3,6,12,24),2] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),2]
  Table5_qwCRPS_2t_Var2[3,] <- G111$mqwCRPS_2t[c(1,2,3,6,12,24),2] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),2]
  Table5_qwCRPS_2t_Var2[4,] <- T000$mqwCRPS_2t[c(1,2,3,6,12,24),2] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),2]
  Table5_qwCRPS_2t_Var2[5,] <- T001$mqwCRPS_2t[c(1,2,3,6,12,24),2] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),2]
  Table5_qwCRPS_2t_Var2[6,] <- T111$mqwCRPS_2t[c(1,2,3,6,12,24),2] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),2]

  Table5_qwCRPS_2t_Var2 <- matrix(sprintf("%.3f", Table5_qwCRPS_2t_Var2), ncol = 6, nrow = 8)
  row.names(Table5_qwCRPS_2t_Var2) <- c("G000", "G001", "G011", "G111", "T000", "T001", "T011", "T111")

  Print_DmTableDens(G000$qwCRPS_2t[,2,],
                    G001$qwCRPS_2t[,2,],
                    G111$qwCRPS_2t[,2,],
                    T000$qwCRPS_2t[,2,],
                    T001$qwCRPS_2t[,2,],
                    T111$qwCRPS_2t[,2,],
                    Table5_qwCRPS_2t_Var2)

  Table5_qwCRPS_2t_Var3 <-  matrix(NA, ncol = 6, nrow = 8)

  Table5_qwCRPS_2t_Var3[1,] <- G000$mqwCRPS_2t[c(1,2,3,6,12,24),3]
  Table5_qwCRPS_2t_Var3[2,] <- G001$mqwCRPS_2t[c(1,2,3,6,12,24),3] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),3]
  Table5_qwCRPS_2t_Var3[3,] <- G111$mqwCRPS_2t[c(1,2,3,6,12,24),3] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),3]
  Table5_qwCRPS_2t_Var3[4,] <- T000$mqwCRPS_2t[c(1,2,3,6,12,24),3] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),3]
  Table5_qwCRPS_2t_Var3[5,] <- T001$mqwCRPS_2t[c(1,2,3,6,12,24),3] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),3]
  Table5_qwCRPS_2t_Var3[6,] <- T111$mqwCRPS_2t[c(1,2,3,6,12,24),3] - G000$mqwCRPS_2t[c(1,2,3,6,12,24),3]

  Table5_qwCRPS_2t_Var3 <- matrix(sprintf("%.3f", Table5_qwCRPS_2t_Var3), ncol = 6, nrow = 8)
  row.names(Table5_qwCRPS_2t_Var3) <- c("G000", "G001", "G011", "G111", "T000", "T001", "T011", "T111")

  Print_DmTableDens(G000$qwCRPS_2t[,3,],
                    G001$qwCRPS_2t[,3,],
                    G111$qwCRPS_2t[,3,],
                    T000$qwCRPS_2t[,3,],
                    T001$qwCRPS_2t[,3,],
                    T111$qwCRPS_2t[,3,],
                    Table5_qwCRPS_2t_Var3)
}

setwd("/home/hoanguc3m/Dropbox/WP11/")
{
  library(gridExtra)
  library(ggthemr)
  ggthemr('light')

  t_step <- 3
  df.G000 <- cbind(G000$emp_CDF[t_step,1,],
                    G000$emp_CDF[t_step,2,],
                    G000$emp_CDF[t_step,3,])
  colnames(df.G000) <-  c("TBILL 3M", "S10-3M", "BAA spread")
  df.G000 <- as.data.frame(df.G000)

  p1 <- ggplot(df.G000, aes(x=`TBILL 3M`)) + geom_histogram(binwidth=0.11) + xlab("TBILL 3M") + ylab("G000")  + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p2 <- ggplot(df.G000, aes(x=`S10-3M`)) + geom_histogram(binwidth=0.11) + xlab("S10-3M") + ylab("") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p3 <- ggplot(df.G000, aes(x=`BAA spread`)) + geom_histogram(binwidth=0.11) + xlab("BAA spread") + ylab("") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")


  df.G001 <- t(G001$emp_CDF[t_step,,])
  colnames(df.G001) <- c("TBILL 3M", "S10-3M", "BAA spread")
  df.G001 <- as.data.frame(df.G001)

  p4 <- ggplot(df.G001, aes(x=`TBILL 3M`)) + geom_histogram(binwidth=0.11) + ylab("G001") + xlab("TBILL 3M") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p5 <- ggplot(df.G001, aes(x=`S10-3M`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("S10-3M") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p6 <- ggplot(df.G001, aes(x=`BAA spread`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("BAA spread") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")


  df.G111 <- t(G111$emp_CDF[t_step,,])
  colnames(df.G111) <- c("TBILL 3M", "S10-3M", "BAA spread")
  df.G111 <- as.data.frame(df.G111)

  p7 <- ggplot(df.G111, aes(x=`TBILL 3M`)) + geom_histogram(binwidth=0.11) + ylab("G111") + xlab("TBILL 3M") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p8 <- ggplot(df.G111, aes(x=`S10-3M`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("S10-3M") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p9 <- ggplot(df.G111, aes(x=`BAA spread`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("BAA spread") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")


  df.T000 <- cbind(T000$emp_CDF[t_step,1,],
                    T000$emp_CDF[t_step,2,],
                    T000$emp_CDF[t_step,3,])
  colnames(df.T000) <-  c("TBILL 3M", "S10-3M", "BAA spread")
  df.T000 <- as.data.frame(df.T000)

  p10 <- ggplot(df.T000, aes(x=`TBILL 3M`)) + geom_histogram(binwidth=0.11) + ylab("T000") + xlab("TBILL 3M") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p11 <- ggplot(df.T000, aes(x=`S10-3M`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("S10-3M") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p12 <- ggplot(df.T000, aes(x=`BAA spread`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("BAA spread") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")


  df.T001 <- t(T001$emp_CDF[t_step,,])
  colnames(df.T001) <- c("TBILL 3M", "S10-3M", "BAA spread")
  df.T001 <- as.data.frame(df.T001)

  p13 <- ggplot(df.T001, aes(x=`TBILL 3M`)) + geom_histogram(binwidth=0.11) + ylab("T001") + xlab("TBILL 3M") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p14 <- ggplot(df.T001, aes(x=`S10-3M`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("S10-3M") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p15 <- ggplot(df.T001, aes(x=`BAA spread`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("BAA spread") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")


  df.T111 <- t(T111$emp_CDF[t_step,,])
  colnames(df.T111) <- c("TBILL 3M", "S10-3M", "BAA spread")
  df.T111 <- as.data.frame(df.T111)

  p16 <- ggplot(df.T111, aes(x=`TBILL 3M`)) + geom_histogram(binwidth=0.11) + ylab("T111") + xlab("TBILL 3M") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p17 <- ggplot(df.T111, aes(x=`S10-3M`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("S10-3M") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")
  p18 <- ggplot(df.T111, aes(x=`BAA spread`)) + geom_histogram(binwidth=0.11) + ylab("") + xlab("BAA spread") + geom_hline(yintercept=237/10, linetype="dashed", color = "red")


  pdf(file='img/PIT.pdf', width = 12, height = 16)

  par(mar=c(2,1,3,5))

  grid.arrange(p1, p2, p3,
               p4, p5, p6,
               p7, p8, p9,
               p10, p11, p12,
               p13, p14, p15,
               p16, p17, p18, nrow = 6, ncol = 3)

  dev.off()

}

# LP plot
################################################
{
  recessions.df = read.table(textConnection(
    "Peak, Trough
  2001-03-01, 2001-11-01
  2007-12-01, 2009-06-01"), sep=',',
    colClasses=c('Date', 'Date'), header=TRUE)


  t_pred <- 3
  time_recursive <- seq(as.Date("2000/01/01"), as.Date("2019/09/01"), "month") + 365 * t_pred/12
  time_full <- seq(as.Date("1953-04-01"), as.Date("2020-09-01"), "month")

  load("/home/hoanguc3m/MEGA/WP11/G000/Recursive_Gaussian_M000_T562.RData")
  out_recursive$y_obs_future
  tobs <- which(y[,3] == out_recursive$y_obs_future[t_pred,3]) # 54
  tobs <- which(y[,1] == out_recursive$y_obs_future[t_pred,1]) # 54
  tobs <- 559 + t_pred - 1

  dfLP <- data.frame(date = time_recursive,
                     LP_TBILL = G001$log_pred[t_pred,1,] - G000$log_pred[t_pred,1,],
                     LP_S10 = G001$log_pred[t_pred,2,] - G000$log_pred[t_pred,2,],
                     LP_BAA = G001$log_pred[t_pred,3,] - G000$log_pred[t_pred,3,],
                     TBILL = y[tobs:(tobs+length(time_recursive)-1),1],
                     S10 = y[tobs:(tobs+length(time_recursive)-1),2],
                     BAA = y[tobs:(tobs+length(time_recursive)-1),3]
  )
  library(dplyr)
  library(ggplot2)
  dfLP <- mutate( dfLP,
                  upT=ifelse(LP_TBILL>0,LP_TBILL,0),
                  downT=ifelse(LP_TBILL<0,LP_TBILL,0),
                  up10=ifelse(LP_S10>0,LP_S10,0),
                  down10=ifelse(LP_S10<0,LP_S10,0),
                  upS=ifelse(LP_BAA>0,LP_BAA,0),
                  downS=ifelse(LP_BAA<0,LP_BAA,0)
  )



  p1 <- ggplot(data=dfLP,aes(x=date,y=LP_TBILL))+
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) +
    geom_line(color="black")+
    geom_line(linetype=2,aes(y=TBILL/100)) +
    geom_ribbon(aes(ymin=LP_TBILL,ymax=downT),fill="#d73027",alpha=0.5)+
    geom_ribbon(aes(ymin=LP_TBILL,ymax=upT),fill="#4575b4",alpha=0.5)  +
    scale_x_date(date_breaks="5 years",date_labels="%Y")+
    theme_minimal(base_size=10)+
    theme(legend.position="top",
          plot.caption=element_text(hjust=0),
          plot.subtitle=element_text(face="italic"),
          plot.title=element_text(size=16,face="bold"))+
    labs(x="",y="TBILL 3M")
  p2 <- ggplot(data=dfLP,aes(x=date,y=LP_S10))+
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) +
    geom_line(color="black")+
    geom_line(linetype=2,aes(y=S10/100)) +
    geom_ribbon(aes(ymin=LP_S10,ymax=down10),fill="#d73027",alpha=0.5)+
    geom_ribbon(aes(ymin=LP_S10,ymax=up10),fill="#4575b4",alpha=0.5)  +
    scale_x_date(date_breaks="5 years",date_labels="%Y")+
    theme_minimal(base_size=10)+
    theme(legend.position="top",
          plot.caption=element_text(hjust=0),
          plot.subtitle=element_text(face="italic"),
          plot.title=element_text(size=16,face="bold"))+
    labs(x="",y="S10-3M")


  p3 <- ggplot(data=dfLP,aes(x=date,y=LP_BAA))+
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) +
    geom_line(color="black")+
    geom_line(linetype=2,aes(y=BAA/100)) +
    geom_ribbon(aes(ymin=LP_BAA,ymax=downS),fill="#d73027",alpha=0.5)+
    geom_ribbon(aes(ymin=LP_BAA,ymax=upS),fill="#4575b4",alpha=0.5)  +
    scale_x_date(date_breaks="5 years",date_labels="%Y")+
    theme_minimal(base_size=10)+
    theme(legend.position="top",
          plot.caption=element_text(hjust=0),
          plot.subtitle=element_text(face="italic"),
          plot.title=element_text(size=16,face="bold"))+
    labs(x="",y="BAA spread")

  setwd("/home/hoanguc3m/Dropbox/WP11/")
  library(gridExtra)
  library(ggthemr)
  ggthemr('light')

  pdf(file='img/diffLPG001-G000.pdf', width = 8, height = 9)

  par(mar=c(2,1,3,5))
  grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
  dev.off()
}




############################################################
# Cumulative log Bayes factors
############################################################

{
  library(readxl)
  library(Matrix)
  library(fattvpVAR)
  library(profvis)
  library(invgamma)
  setwd("/home/hoanguc3m/Dropbox/WP11/")

  dataraw <- read_excel("/home/hoanguc3m/MEGA/HybridVAR/EconLetter/temp/Data210927.xlsx",
                        col_types = c("text", "numeric", "numeric", "numeric"))
  priors <- list(  hyper_ab = 1,
                   hyper_h = 1)

  p <- 3 # number of lags
  atT <- 821


  y <- data.matrix(dataraw[(p+1):atT,c(2:4)])
  y0 <- data.matrix(dataraw[1:p,c(2:4)])

  recessions.df = read.table(textConnection(
    "Peak, Trough
    2001-03-01, 2001-11-01
    2007-12-01, 2009-06-01"), sep=',',
    colClasses=c('Date', 'Date'), header=TRUE)

  setwd("/home/hoanguc3m/Dropbox/WP11/")

  t_pred <- 3
  load("/home/hoanguc3m/Dropbox/WP11/G000/Recursive_Gaussian_M000_T562.RData")
  out_recursive$y_obs_future
  tobs <- which(y[,3] == out_recursive$y_obs_future[t_pred,3]) # 54
  tobs <- which(y[,1] == out_recursive$y_obs_future[t_pred,1]) # 54
  tobs <- 559 + t_pred - 1

  dfLP <- data.frame(date = time_recursive,
                     LP_TBILL = cumsum(G001$log_pred[t_pred,1,] - G000$log_pred[t_pred,1,]),
                     LP_S10 = cumsum(G001$log_pred[t_pred,2,] - G000$log_pred[t_pred,2,]),
                     LP_BAA = cumsum(G001$log_pred[t_pred,3,] - G000$log_pred[t_pred,3,]),
                     TBILL = y[tobs:(tobs+length(time_recursive)-1),1],
                     S10 = y[tobs:(tobs+length(time_recursive)-1),2],
                     BAA = y[tobs:(tobs+length(time_recursive)-1),3]
  )


  library(dplyr)
  library(ggplot2)
  dfLP <- mutate( dfLP,
                  upT=ifelse(LP_TBILL>0,LP_TBILL,0),
                  downT=ifelse(LP_TBILL<0,LP_TBILL,0),
                  up10=ifelse(LP_S10>0,LP_S10,0),
                  down10=ifelse(LP_S10<0,LP_S10,0),
                  upS=ifelse(LP_BAA>0,LP_BAA,0),
                  downS=ifelse(LP_BAA<0,LP_BAA,0)
  )



  p1 <- ggplot(data=dfLP,aes(x=date,y=LP_TBILL))+
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) +
    geom_line(color="black")+
    geom_line(linetype=2,aes(y=TBILL)) +
    geom_ribbon(aes(ymin=LP_TBILL,ymax=downT),fill="#d73027",alpha=0.5)+
    geom_ribbon(aes(ymin=LP_TBILL,ymax=upT),fill="#4575b4",alpha=0.5)  +
    scale_x_date(date_breaks="5 years",date_labels="%Y")+
    theme_minimal(base_size=10)+
    theme(legend.position="top",
          plot.caption=element_text(hjust=0),
          plot.subtitle=element_text(face="italic"),
          plot.title=element_text(size=16,face="bold"))+
    labs(x="",y="TBILL 3M")
  p2 <- ggplot(data=dfLP,aes(x=date,y=LP_S10))+
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) +
    geom_line(color="black")+
    geom_line(linetype=2,aes(y=S10)) +
    geom_ribbon(aes(ymin=LP_S10,ymax=down10),fill="#d73027",alpha=0.5)+
    geom_ribbon(aes(ymin=LP_S10,ymax=up10),fill="#4575b4",alpha=0.5)  +
    scale_x_date(date_breaks="5 years",date_labels="%Y")+
    theme_minimal(base_size=10)+
    theme(legend.position="top",
          plot.caption=element_text(hjust=0),
          plot.subtitle=element_text(face="italic"),
          plot.title=element_text(size=16,face="bold"))+
    labs(x="",y="S10-3M")


  p3 <- ggplot(data=dfLP,aes(x=date,y=LP_BAA))+
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) +
    geom_line(color="black")+
    geom_line(linetype=2,aes(y=BAA)) +
    geom_ribbon(aes(ymin=LP_BAA,ymax=downS),fill="#d73027",alpha=0.5)+
    geom_ribbon(aes(ymin=LP_BAA,ymax=upS),fill="#4575b4",alpha=0.5)  +
    scale_x_date(date_breaks="5 years",date_labels="%Y")+
    theme_minimal(base_size=10)+
    theme(legend.position="top",
          plot.caption=element_text(hjust=0),
          plot.subtitle=element_text(face="italic"),
          plot.title=element_text(size=16,face="bold"))+
    labs(x="",y="BAA spread")

  library(gridExtra)
  library(ggthemr)
  ggthemr('light')

  pdf(file='img/cumLogBF3-G1-G0.pdf', width = 8, height = 9)
  #pdf(file='img/diff_log.pdf', width = 11, height = 6)

  par(mar=c(2,1,3,5))
  grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
  dev.off()
}

{
  library(readxl)
  library(Matrix)
  library(fattvpVAR)
  library(profvis)
  library(invgamma)
  setwd("/home/hoanguc3m/Dropbox/WP11/")
  dataraw <- read_excel("/home/hoanguc3m/MEGA/HybridVAR/EconLetter/temp/Data210927.xlsx",
                        col_types = c("text", "numeric", "numeric", "numeric"))
  priors <- list(  hyper_ab = 1,
                   hyper_h = 1)

  p <- 3 # number of lags
  atT <- 821


  y <- data.matrix(dataraw[(p+1):atT,c(2:4)])
  y0 <- data.matrix(dataraw[1:p,c(2:4)])

  recessions.df = read.table(textConnection(
    "Peak, Trough
    2001-03-01, 2001-11-01
    2007-12-01, 2009-06-01"), sep=',',
    colClasses=c('Date', 'Date'), header=TRUE)

  setwd("/home/hoanguc3m/Dropbox/WP11/")

  t_pred <- 3
  load("/home/hoanguc3m/Dropbox/WP11/G000/Recursive_Gaussian_M000_T562.RData")
  out_recursive$y_obs_future
  tobs <- which(y[,3] == out_recursive$y_obs_future[t_pred,3]) # 54
  tobs <- which(y[,1] == out_recursive$y_obs_future[t_pred,1]) # 54
  tobs <- 559 + t_pred - 1

  dfLP <- data.frame(date = time_recursive,
                     LP_TBILL = cumsum(T001$log_pred[t_pred,1,] - G001$log_pred[t_pred,1,]),
                     LP_S10 = cumsum(T001$log_pred[t_pred,2,] - G001$log_pred[t_pred,2,]),
                     LP_BAA = cumsum(T001$log_pred[t_pred,3,] - G001$log_pred[t_pred,3,]),
                     TBILL = y[tobs:(tobs+length(time_recursive)-1),1],
                     S10 = y[tobs:(tobs+length(time_recursive)-1),2],
                     BAA = y[tobs:(tobs+length(time_recursive)-1),3]
  )


  library(dplyr)
  library(ggplot2)
  dfLP <- mutate( dfLP,
                  upT=ifelse(LP_TBILL>0,LP_TBILL,0),
                  downT=ifelse(LP_TBILL<0,LP_TBILL,0),
                  up10=ifelse(LP_S10>0,LP_S10,0),
                  down10=ifelse(LP_S10<0,LP_S10,0),
                  upS=ifelse(LP_BAA>0,LP_BAA,0),
                  downS=ifelse(LP_BAA<0,LP_BAA,0)
  )



  p1 <- ggplot(data=dfLP,aes(x=date,y=LP_TBILL))+
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) +
    geom_line(color="black")+
    geom_line(linetype=2,aes(y=TBILL)) +
    geom_ribbon(aes(ymin=LP_TBILL,ymax=downT),fill="#d73027",alpha=0.5)+
    geom_ribbon(aes(ymin=LP_TBILL,ymax=upT),fill="#4575b4",alpha=0.5)  +
    scale_x_date(date_breaks="5 years",date_labels="%Y")+
    theme_minimal(base_size=10)+
    theme(legend.position="top",
          plot.caption=element_text(hjust=0),
          plot.subtitle=element_text(face="italic"),
          plot.title=element_text(size=16,face="bold"))+
    labs(x="",y="TBILL 3M")
  p2 <- ggplot(data=dfLP,aes(x=date,y=LP_S10))+
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) +
    geom_line(color="black")+
    geom_line(linetype=2,aes(y=S10)) +
    geom_ribbon(aes(ymin=LP_S10,ymax=down10),fill="#d73027",alpha=0.5)+
    geom_ribbon(aes(ymin=LP_S10,ymax=up10),fill="#4575b4",alpha=0.5)  +
    scale_x_date(date_breaks="5 years",date_labels="%Y")+
    theme_minimal(base_size=10)+
    theme(legend.position="top",
          plot.caption=element_text(hjust=0),
          plot.subtitle=element_text(face="italic"),
          plot.title=element_text(size=16,face="bold"))+
    labs(x="",y="S10-3M")


  p3 <- ggplot(data=dfLP,aes(x=date,y=LP_BAA))+
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) +
    geom_line(color="black")+
    geom_line(linetype=2,aes(y=BAA)) +
    geom_ribbon(aes(ymin=LP_BAA,ymax=downS),fill="#d73027",alpha=0.5)+
    geom_ribbon(aes(ymin=LP_BAA,ymax=upS),fill="#4575b4",alpha=0.5)  +
    scale_x_date(date_breaks="5 years",date_labels="%Y")+
    theme_minimal(base_size=10)+
    theme(legend.position="top",
          plot.caption=element_text(hjust=0),
          plot.subtitle=element_text(face="italic"),
          plot.title=element_text(size=16,face="bold"))+
    labs(x="",y="BAA spread")

  library(gridExtra)
  library(ggthemr)
  ggthemr('light')

  pdf(file='img/cumLogBF3-S1-G1.pdf', width = 8, height = 9)
  #pdf(file='img/diff_log.pdf', width = 11, height = 6)

  par(mar=c(2,1,3,5))
  grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
  dev.off()
}

