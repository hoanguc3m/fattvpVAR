#' @export
MLChib_TVPGSV <- function(Chain, samples = 20000, burnin = 10000, thin = 1, numCores = 4){
  # samples = 6000; burnin = 3000; thin = 1;
  # numCores <- 16;

  P1 <- P2 <- P3 <- P4 <- P5 <- 0

  priors <- Chain$data$priors
  inits <- Chain$data$inits
  dist <- Chain$data$dist

  B_mat <- Chain$store_beta0
  A_mat <- Chain$store_alp0
  SigmaB_mat <- Chain$store_Sigbeta
  SigmaA_mat <- Chain$store_Sigalp
  Sigmah_mat <- Chain$store_Sigh

  #Nu_mat <- Chain$store_nu
  H0_mat <- Chain$store_h0

  B_med <- apply(B_mat, MARGIN = 2, FUN = mean)
  A_med <- apply(A_mat, MARGIN = 2, FUN = mean)
  SigmaB_med <- apply(sqrt(SigmaB_mat), MARGIN = 2, FUN = mean)^2
  SigmaA_med <- apply(sqrt(SigmaA_mat), MARGIN = 2, FUN = mean)^2
  Sigmah_med <- apply(sqrt(Sigmah_mat), MARGIN = 2, FUN = mean)^2
  #Nu_med <- apply(Nu_mat, MARGIN = 2, FUN = mean)
  H0_med <- apply(H0_mat, MARGIN = 2, FUN = mean)

  inits$samples <- samples # Reset
  inits$burnin <- burnin # Reset
  inits$thin <- thin # Reset

  # Commom parameters
  K <- ncol(Chain$data$y)
  p <- Chain$data$p
  inits$a0 <- A_med
  inits$b0 <- B_med
  inits$SigmaB <- SigmaB_med
  inits$SigmaA <- SigmaA_med
  inits$sigma_h <- Sigmah_med
  inits$h0 <- H0_med

  str(inits)

  Start = Sys.time()

  ChibLLP_chain <- ChibLLP_TVP_GSV(Chain, ndraws = 1000, numCores = numCores)
  RhpcBLASctl::blas_set_num_threads(RhpcBLASctl::get_num_cores())

  y <- Chain$data$y
  y0 <- Chain$data$y0

  cat("Pi( B*, A*, SigmaB*, SigmaA* | y, SigmaH*) \n ")
  P1 <- Chib_fitTVP_GSV(y, K = K, p = p, y0 = y0, priors = priors, inits = inits,
                         samples = samples, burnin = burnin, thin = thin,
                         fix_B = FALSE, fix_Sigma = TRUE,
                         cal_B = TRUE,  cal_Sigma = FALSE)
  cat("Pi( SigmaH* | y) \n ")
  P2 <- Chib_fitTVP_GSV(y, K = K, p = p, y0 = y0, priors = priors, inits = inits,
                         samples = samples, burnin = burnin, thin = thin,
                         fix_B = FALSE, fix_Sigma = FALSE,
                         cal_B = FALSE, cal_Sigma = TRUE)

  CML_Chain <- ChibLLP_chain$LL - sum(logmeanexp(P1)) - logmeanexp(P2)


  elapsedTime = Sys.time() - Start
  print(elapsedTime)
  out <- list(CML = CML_Chain,
              LLP = ChibLLP_chain$LL,
              lprior = ChibLLP_chain$lprior,
              lc = ChibLLP_chain$lc,
              LLP_chain = ChibLLP_chain,
              aP1 = logmeanexp(P1), aP2 = logmeanexp(P2), aP3 = logmeanexp(P3),
              aP4 = sum(logmeanexp(P4)), aP5 = sum(logmeanexp(P5)),
              P1 = P1, P2 = P2, P3 = P3, P4 = P4, P5 = P5,
              esttime = elapsedTime)

  return(out)
}
#' @export
Chib_fitTVP_GSV <- function(y, K, p, y0 = NULL, priors = NULL, inits = NULL,
                             samples = 1000, burnin = 0, thin = 1,
                             fix_B = FALSE, fix_Sigma = FALSE,
                             cal_B = FALSE, cal_Sigma = FALSE){

  Y0 <- y0      # p x K matrix
  shortY <- y   # T x K matrix
  t_max <- nrow(shortY)
  K <- ncol(shortY)
  Y <- matrixcalc::vec(t(shortY)) # Y stack by obs and time


  # samples <- inits$samples
  # burnin <- inits$burnin
  # thin <- inits$thin
  if (cal_B) {
    samples <- 120000
    burnin <- 20000
    thin <- 10
  }
  samplesDiv <- floor(samples / thin)
  is_tv <- inits$is_tv

  #is_tv = [1 1 1]'; # n-vector that denotes which equations have TVP; must have the same dim as the VAR
  # e.g., [1 0 1]: eq 1 and 3 have TVP, but eq 2 has const parameters

  k_alp <- K*(K-1)/2 # dimension of the impact matrix
  k_beta <- K^2*p + K # number of VAR coefficients
  k_beta_div_K <- K*p + 1
  n_tv <- sum(is_tv)# # of time-varying equations

  ## prior
  abeta0 <- priors$b0
  Vbeta0 <- priors$V_b_prior
  aalp0 <- priors$a0
  Valp0 <- priors$V_a0_prior

  ah0 <- c(log(priors$sigma^2))
  Vh0 <- 4*matrix(1, nrow = K,ncol = 1)

  # Variance of the tvp
  hyper_ab <- priors$hyper_ab
  hyper_h <- priors$hyper_h
  Sbeta0 <- hyper_ab*matrix(1, nrow = k_beta, ncol = 1)
  Sbeta0[seq(1,k_beta, by =  K*p+1)] <- hyper_ab # intercepts more diffused prior 10*hyper_ab
  Sbeta0 <- Sbeta0*kronecker(is_tv, matrix(1, nrow = K*p+1, ncol = 1) )# set to 0 for time-invariant equations

  Salp0 <- hyper_ab*matrix(1, nrow = k_alp, ncol = 1)

  k_a_eq <- seq(0, K-1)
  id_a <- cbind(cumsum(k_a_eq) - k_a_eq + 1,cumsum(k_a_eq))
  count_seqa <- list(); count_seqb <- list()
  count_seqa[[1]] <- 0; count_seqb[[1]] <- seq(1, k_beta_div_K)
  for (ii in c(2:K)){
    count_seqa[[ii]] <- seq(id_a[ii,1], id_a[ii,2])
    count_seqb[[ii]] <- ((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)

    if (is_tv[ii] == 0){
      Salp0[count_seqa[[ii]]] <- 0 # set to 0 for time-invariant equations
    }
  }

  #Sh0 <- hyper_h*matrix(1, nrow = K, ncol = 1)


  ## compute and define a few things
  tmpY <- rbind(Y0, shortY)
  X2 <- matrix(0, nrow = t_max, ncol = K*p)
  for (ii in 1:p){
    X2[,((ii-1)*K+1):(ii*K)] <- tmpY[(p-ii+1):(t_max+p-ii),]
  }
  X2 <- cbind(rep(1, t_max), X2)

  idx_b_tv <- (kronecker(is_tv,matrix(1,nrow = K*p+1,ncol = 1))==1)   # index for time-varying betas
  idx_a_tv <- matrix(FALSE, nrow = k_alp, ncol = 1)            # construct index for time-varying alphas
  for (ii in 2:K){
    if (is_tv[ii] == 1){
      idx_a_tv[ count_seqa[[ii]]] <- TRUE
    }
  }


  ## initialize the Markov chain

  Sigbeta <- inits$SigmaB
  Sigalp <- inits$SigmaA
  Sigh <- diag(inits$sigma_h, K)

  beta0 <- inits$b0
  alp0 <- inits$a0
  h0 <- inits$h0

  beta <- matrix(0, nrow = t_max, ncol = k_beta)
  alp <- matrix(0, nrow = t_max, ncol = k_alp)
  h <- kronecker(t(h0), rep(1,t_max))

  ## MCMC starts here
  set.seed(NULL)
  sprintf('Starting MCMC for (hybrid) Gaussian TVP-VAR-SV with ID: %s ', paste(as.character(is_tv), collapse = ' ') )
  start_time <- Sys.time()

  lpost <- rep(0, (samples - burnin)%/% thin)
  if (cal_B) lpost <- matrix(0,nrow = (samples - burnin)%/% thin, ncol = K)

  for (j in c(1:samples)){
    if(!fix_B) {
      # sample alp0 and beta0 - equation by equation
      U1 <- shortY # U = Y - Z * tilde(theta)
      U2 <- matrix(0, nrow = t_max, ncol = K) # U = Y - X * theta
      U <- matrix(0, nrow = t_max, ncol = K) # U = Y - X * theta - Z * tilde(theta)

      for (ii in 1:K){
        ki <- K*p+1+ii-1 # Number of theta in equation ii
        B_star <- inits$b0[ count_seqb[[ii]] ]
        A_star <- inits$a0[ count_seqa[[ii]] ]
        SigmaB_star <- inits$SigmaB[ count_seqb[[ii]] ]
        SigmaA_star <- inits$SigmaA[ count_seqa[[ii]] ]

        if (ii > 1) {
          X <- cbind(X2, -shortY[,1:(ii-1)])
        } else {
          X <- X2
        }
        if (is_tv[ii] == 1){
          thetai <- cbind(beta[, count_seqb[[ii]] ], alp[,count_seqa[[ii]] ])

          invvol <- as.vector(exp(-h[,ii]/2))
          y.tilde <- as.vector( U1[,ii] ) * invvol
          x.tilde <- cbind(X, X*thetai) * invvol

          V_b_prior_inv <- diag(1/c(Vbeta0[ count_seqb[[ii]] ], Valp0[count_seqa[[ii]] ],
                                    Sbeta0[ count_seqb[[ii]] ], Salp0[count_seqa[[ii]] ] ))
          theta.prior.precmean <- c(abeta0[count_seqb[[ii]] ], aalp0[count_seqa[[ii]] ], rep(0,ki) ) /
            c(Vbeta0[ count_seqb[[ii]] ], Valp0[count_seqa[[ii]] ],
              Sbeta0[ count_seqb[[ii]] ], Salp0[count_seqa[[ii]] ] )

          theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
          while (TRUE) {
            thetai0 <- backsolve( theta.prec.chol,
                                  backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                             upper.tri = T, transpose = T )
                                  + rnorm(ki) )
            ab_sample <- thetai0[1:ki]
            Sab_sample <- thetai0[(ki+1):(2*ki)]
            if (all(Sab_sample > 0) ) {
              break
            }
          }
          beta0[count_seqb[[ii]] ] <- ab_sample[1:k_beta_div_K]
          Sigbeta[count_seqb[[ii]] ] <- Sab_sample[1:k_beta_div_K]^2
          if ( ii > 1){
            alp0[count_seqa[[ii]] ] <- ab_sample[(k_beta_div_K+1):ki]
            Sigalp[count_seqa[[ii]] ] <- Sab_sample[(k_beta_div_K+1):ki]^2
          }
          U2[,ii] <- shortY[,ii] - X %*% ab_sample

        } else {

          invvol <- as.vector(exp(-h[,ii]/2))
          y.tilde <- as.vector( U1[,ii] ) * invvol
          x.tilde <- X * invvol

          V_b_prior_inv <- diag(1/c(Vbeta0[ count_seqb[[ii]] ], Valp0[count_seqa[[ii]] ]))
          theta.prior.precmean <- c(abeta0[count_seqb[[ii]] ], aalp0[count_seqa[[ii]] ] ) /
            c(Vbeta0[ count_seqb[[ii]] ], Valp0[count_seqa[[ii]] ])

          theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
          thetai0 <- backsolve( theta.prec.chol,
                                backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                           upper.tri = T, transpose = T )
                                + rnorm(ki) )
          beta0[count_seqb[[ii]] ] <- thetai0[1:k_beta_div_K];
          if ( ii > 1){
            alp0[count_seqa[[ii]] ] <- thetai0[(k_beta_div_K+1):ki]
          }
          U2[,ii] <- shortY[,ii] - X %*% thetai0


        }


      }
    }

    if (cal_B & (j > burnin) & (j %% thin == 0)){
      for (ii in 1:K){
        B_star <- inits$b0[ count_seqb[[ii]] ]
        A_star <- inits$a0[ count_seqa[[ii]] ]
        SigmaB_star <- inits$SigmaB[ count_seqb[[ii]] ]
        SigmaA_star <- inits$SigmaA[ count_seqa[[ii]] ]

        ki <- K*p+1+ii-1 # Number of theta in equation ii

        if (ii > 1) {
          X <- cbind(X2, -shortY[,1:(ii-1)])
        } else {
          X <- X2
        }
        if (is_tv[ii] == 1){
          # if (ii < 3) {
          ABS_star <- c(B_star, A_star,
                        sqrt(SigmaB_star), sqrt(SigmaA_star) )

          thetai <- cbind(beta[, count_seqb[[ii]] ], alp[,count_seqa[[ii]] ])

          invvol <- as.vector(exp(-h[,ii]/2))
          y.tilde <- as.vector( U1[,ii] ) * invvol
          x.tilde <- cbind(X, X*thetai) * invvol

          V_b_prior_inv <- diag(1/c(Vbeta0[ count_seqb[[ii]] ], Valp0[count_seqa[[ii]] ],
                                    Sbeta0[ count_seqb[[ii]] ], Salp0[count_seqa[[ii]] ] ))
          theta.prior.precmean <- c(abeta0[count_seqb[[ii]] ], aalp0[count_seqa[[ii]] ], rep(0,ki) ) /
            c(Vbeta0[ count_seqb[[ii]] ], Valp0[count_seqa[[ii]] ],
              Sbeta0[ count_seqb[[ii]] ], Salp0[count_seqa[[ii]] ] )

          theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
          b_star <- backsolve( theta.prec.chol,
                                backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                           upper.tri = T, transpose = T ))

          lpost[(j - burnin) %/% thin,ii] <- - length(ABS_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
                              0.5 * t(ABS_star - b_star) %*% (V_b_prior_inv + crossprod(x.tilde)) %*% (ABS_star - b_star) +
                              - sum(log(ABS_star[(ki+1):(2*ki)])) # Jacobian transformation

        } else {
          ABS_star <- c(B_star, A_star)

          invvol <- as.vector(exp(-h[,ii]/2))
          y.tilde <- as.vector( U1[,ii] ) * invvol
          x.tilde <- X * invvol

          V_b_prior_inv <- diag(1/c(Vbeta0[ count_seqb[[ii]] ], Valp0[count_seqa[[ii]] ]))
          theta.prior.precmean <- c(abeta0[count_seqb[[ii]] ], aalp0[count_seqa[[ii]] ] ) /
            c(Vbeta0[ count_seqb[[ii]] ], Valp0[count_seqa[[ii]] ])

          theta.prec.chol <- chol( V_b_prior_inv + crossprod(x.tilde) )
          b_star <- backsolve( theta.prec.chol,
                                backsolve( theta.prec.chol, theta.prior.precmean + crossprod( x.tilde, y.tilde ),
                                           upper.tri = T, transpose = T ))
          lpost[(j - burnin) %/% thin, ii] <- - length(ABS_star) * 0.5 * log(2*pi) + sum(log(diag(theta.prec.chol))) -
                                0.5 * t(ABS_star - b_star) %*% (V_b_prior_inv + crossprod(x.tilde)) %*% (ABS_star - b_star)

        }


      }

    }
    # matrix(beta0, nrow = K, byrow = T); alp0

    # sample latent alp and beta - equation by equation
    U <- U2
    for (ii in 1:K){
      if (is_tv[ii] == 1){
        ki <- K*p+1+ii-1 # Number of theta in equation ii
        if (ii > 1) {
          X <- cbind(X2, -shortY[,1:(ii-1)])
        } else {
          X <- X2
        }

        Sigthetai <- c( Sigbeta[count_seqb[[ii]] ], Sigalp[count_seqa[[ii]] ])
        bigXi <- SURform(X * reprow(sqrt(Sigthetai), t_max) ) # "dgCMatrix"
        Tthetai <- Matrix::Diagonal(t_max*ki) - Matrix::sparseMatrix( i = (ki+1):(t_max*ki),
                                                                      j = 1:((t_max-1)*ki),
                                                                      x = rep(1, (t_max-1)*ki),
                                                                      dims = c(t_max*ki,t_max*ki))

        thetai0 <- rep(0,ki)

        x.tilde = Matrix::t(Matrix::t(bigXi) %*% Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-0.5*h[,ii])))
        TiST_thetai <- Matrix::t(Tthetai) %*% Tthetai

        y.tilde = Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-0.5*h[,ii])) %*% U2[,ii]
        #theta.prior.precmean = Matrix(data = 0, nrow = t_max*length(thetai0), ncol = 1)
        theta.prec.chol <- Matrix::chol( TiST_thetai + Matrix::crossprod(x.tilde) )
        thetai <- Matrix::solve( theta.prec.chol,
                                 Matrix::solve( Matrix::t(theta.prec.chol), Matrix::crossprod( x.tilde, y.tilde ))
                                 + rnorm(t_max*ki) )


        Thetai <- t(matrix(thetai, nrow = ki, ncol = t_max))

        beta[,count_seqb[[ii]] ] <- Thetai[,1:k_beta_div_K]
        if ( ii > 1){
          alp[,count_seqa[[ii]] ] <- Thetai[,(k_beta_div_K+1):ki]
        }
        U[,ii] <- as.matrix(U[,ii] - bigXi %*% thetai)
      }
    }

    # sample latent h
    aux <- sample_h_mod(ytilde = t(U), sigma_h = Sigh, h0_mean = h0,
                        h = t(h), K = K, t_max = t_max, prior = priors)
    h <- t(aux$Sigtdraw)
    h0 <- as.numeric(aux$h0)


    if(!fix_Sigma) {
      if (K>1) {
        sse_2 <- apply( (h[1:t_max,] - rbind(h0,h[1:(t_max-1),]) )^2, MARGIN = 2, FUN = sum)
      } else {
        sse_2 <- sum( (h[1:t_max,] - c(h0,h[1:(t_max-1),]) )^2)
      }

      Sigh <- diag(mapply( GIGrvg::rgig, n = 1, lambda = - (t_max - 1)*0.5, chi = sse_2,
                              psi = 1/priors$hyper_h ), nrow = K)

    }

    if (cal_Sigma & (j > burnin) & (j %% thin == 0)){

      if (K>1) {
        sse_2 <- apply( (h[1:t_max,] - rbind(h0,h[1:(t_max-1),]) )^2, MARGIN = 2, FUN = sum)
      } else {
        sse_2 <- sum( (h[1:t_max,] - c(h0,h[1:(t_max-1),]) )^2)
      }


      lpost[(j - burnin) %/% thin] <- sum(mapply(GIGrvg::dgig, x = inits$sigma_h, lambda = - (t_max - 1)*0.5,
                                                 chi = sse_2, psi = 1/priors$hyper_h, log = T))
    }

    if ( j %% 1000 == 0){
      cat(' Iteration ... ', j, '\n')
    }

  }
  end_time <- Sys.time()
  print( end_time - start_time)


  return(lpost)
}

#' @export
ChibLLP_TVP_GSV <- function(Chain, ndraws = 1000, numCores = NULL){
  RhpcBLASctl::blas_set_num_threads(1)
  priors <- Chain$data$priors
  data = Chain$data
  Y0 <- data$y0
  shortY <- data$y
  t_max <- nrow(shortY)
  K <- ncol(shortY)
  Y <- matrixcalc::vec(t(shortY))
  p <- data$p
  is_tv <- data$inits$is_tv


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

  B_med <- apply(Chain$store_beta0, MARGIN = 2, FUN = mean)
  A_med <- apply(Chain$store_alp0, MARGIN = 2, FUN = mean)
  SigmaB_med <- apply(sqrt(Chain$store_Sigbeta), MARGIN = 2, FUN = mean)^2
  SigmaA_med <- apply(sqrt(Chain$store_Sigalp), MARGIN = 2, FUN = mean)^2
  Sigmah_med <- apply(sqrt(Chain$store_Sigh), MARGIN = 2, FUN = mean)^2
  H0_med <- apply(Chain$store_h0, MARGIN = 2, FUN = mean)

  ## prior borrow from inference
  abeta0 <- priors$b0
  Vbeta0 <- priors$V_b_prior
  aalp0 <- priors$a0
  Valp0 <- priors$V_a0_prior

  ah0 <- c(log(priors$sigma^2))
  Vh0 <- 4*matrix(1, nrow = K,ncol = 1)

  sum_log_prior <- sum(dgamma(Sigmah_med, shape = 0.5, rate = 0.5 / priors$hyper_h, log = T),
                                dgamma(SigmaB_med[idx_b_tv], shape = 0.5, rate = 0.5 / priors$hyper_ab, log = T),
                                dgamma(SigmaA_med[idx_a_tv], shape = 0.5, rate = 0.5 / priors$hyper_ab, log = T),
                                #dnorm(H0_med, mean = ah0, sd = sqrt(Vh0), log = T),
                                dnorm(B_med, mean = abeta0, sd = sqrt(Vbeta0), log = T),
                                dnorm(A_med, mean = aalp0, sd = sqrt(Valp0), log = T) )

  dpriorABS <- rep(0, K)
  for (ii in c(1:K) ){
    dpriorABS[ii] <- sum(dgamma(SigmaB_med[count_seqb[[ii]]] [idx_b_tv[count_seqb[[ii]]]], shape = 0.5, rate = 0.5 / priors$hyper_ab, log = T)) +
      sum(dgamma(SigmaA_med[count_seqa[[ii]]] [idx_a_tv[count_seqa[[ii]]]], shape = 0.5, rate = 0.5 / priors$hyper_ab, log = T)) +
      sum(dnorm(B_med[count_seqb[[ii]] ], mean = abeta0, sd = sqrt(Vbeta0), log = T)) +
      sum(dnorm(A_med[count_seqa[[ii]] ], mean = aalp0, sd = sqrt(Valp0), log = T))
  }


  #sum_log = rep(0, M);
  RhpcBLASctl::blas_set_num_threads(1)

  sum_log <- parallel::mclapply(1:ndraws,
                                FUN = function(j) {
                                  Sigbeta = SigmaB_med
                                  Sigalp = SigmaA_med
                                  Sigh = Sigmah_med
                                  beta0 = B_med
                                  alp0 = A_med
                                  h0 = H0_med

                                  llike = 0
                                  for (ii in c(1:K)){
                                    ki <- K*p+1+ii-1 # Number of theta in equation ii

                                    if (ii > 1) {
                                      X <- cbind(X2, -shortY[,1:(ii-1)])
                                    } else {
                                      X <- X2
                                    }

                                    if (is_tv[ii] == 1){
                                      Sigthetai <- c( Sigbeta[count_seqb[[ii]] ], Sigalp[count_seqa[[ii]] ])
                                      bigXi <- SURform(X * reprow(sqrt(Sigthetai), t_max) ) # "dgCMatrix"
                                      thetaXi0 <- c(beta0[count_seqb[[ii]] ], alp0[count_seqa[[ii]] ])
                                      Yi <- shortY[,ii] - X %*% thetaXi0

                                      thetai0 <- rep(0, ki)
                                      Sigthetai <- rep(1, ki)
                                      llikei = Chib_intlike_tvpsv(Yi = Yi, Sigthetai = Sigthetai, Sig_hi = Sigh[ii], bigXi = bigXi, h0i = h0[ii], thetai0 = thetai0)
                                    }  else {
                                      bigXi = X
                                      thetai0 <- c(beta0[count_seqb[[ii]] ], alp0[ count_seqa[[ii]] ])
                                      llikei = Chib_intlike_varsv(Yi = shortY[,ii], thetai0 = thetai0, Sig_hi = Sigh[ii], bigXi = bigXi, h0i = h0[ii])
                                    }

                                    llike = llike + llikei

                                  }

                                  llike

                                }, mc.cores = numCores)
  sum_log = unlist(sum_log)

  store_w = sum_log
  shortw =  matrix(store_w, nrow = ndraws/20, ncol = 20)
  maxw = apply(shortw, MARGIN = 2 , FUN = max)
  bigml = log( apply(exp(shortw-reprow(maxw,ndraws/20)), MARGIN = 2, FUN = mean )) + maxw
  ml = mean(bigml)
  mlstd = sd(bigml)/sqrt(20)
  return( list( lprior = as.numeric(sum_log_prior),
                LL = as.numeric(ml + sum_log_prior),
                std = mlstd,
                lc = ml,
                sum_log = sum_log,
                dpriorABS = dpriorABS))
}


