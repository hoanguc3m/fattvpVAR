#' @export
fitTVPStudentSV <- function(y, y0, p, priors, inits){
  TARGACCEPT = 0.3
  batchlength = 10

  Y0 <- y0      # p x K matrix
  shortY <- y   # T x K matrix
  t_max <- nrow(shortY)
  K <- ncol(shortY)
  Y <- matrixcalc::vec(t(shortY)) # Y stack by obs and time


  samples <- inits$samples
  burnin <- inits$burnin
  thin <- inits$thin
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
  Sh0 <- matrix(1, nrow = K, ncol = 1)
  nu_gam_a <- 2; nu_gam_b <- 0.1;


  # Variance of the tvp
  hyper_ab <- 1
  hyper_h <- 1

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

  ## initialize for storage
  store_alp <- array(0, dim = c(samplesDiv,t_max,k_alp))
  store_beta <- array(0, dim = c(samplesDiv,t_max,k_beta))
  store_h <- array(0, dim = c(samplesDiv,t_max,K))
  store_Sigbeta <- array(0, dim = c(samplesDiv,k_beta))
  store_Sigalp <- array(0, dim = c(samplesDiv,k_alp))
  store_Sigh <- array(0, dim = c(samplesDiv,K))
  store_beta0 <- array(0, dim = c(samplesDiv,k_beta))
  store_alp0 <- array(0, dim = c(samplesDiv,k_alp))
  store_h0 <- array(0, dim = c(samplesDiv,K))
  store_w <- array(0, dim = c(samplesDiv,t_max,K))
  store_nu <- array(0, dim = c(samplesDiv,K))

  logsigma_nu <- rep(0,K)
  acount_nu <- rep(0,K)

  ## initialize the Markov chain

  Sigbeta <- matrix(0, nrow = k_beta, ncol = 1)
  Sigbeta[idx_b_tv] <- .01 * matrix(runif( k_beta_div_K*n_tv), ncol = 1)
  Sigalp <- matrix(0, nrow = k_alp, ncol = 1)
  Sigalp[idx_a_tv] <- .01 * matrix(runif(sum(idx_a_tv)), ncol = 1)
  Sigh <- .1*matrix( runif(K), ncol = 1)

  beta0 <- abeta0
  alp0 <- aalp0
  h0 <- matrix(log( diag(var(shortY))),ncol = 1)

  beta <- matrix(0, nrow = t_max, ncol = k_beta)
  alp <- matrix(0, nrow = t_max, ncol = k_alp)
  h <- kronecker(t(h0), rep(1,t_max))

  nu <- rep(10,K)
  w <- matrix(1, nrow = t_max, ncol = K)
  sqrt_w <- sqrt(w)


  ## MCMC starts here
  set.seed(NULL)
  sprintf('Starting MCMC for (hybrid) Student TVP-VAR-SV with ID: %s ', paste(as.character(is_tv), collapse = ' ') )
  start_time <- Sys.time()

  for (isim in c(1:(samples + burnin))){
    # sample alp0 and beta0 - equation by equation
    U1 <- shortY # U = Y - Z * tilde(theta)
    U2 <- matrix(0, nrow = t_max, ncol = K) # U = Y - X * theta
    U <- matrix(0, nrow = t_max, ncol = K) # U = Y - X * theta - Z * tilde(theta)
    E <- matrix(0, nrow = t_max, ncol = K)
    for (ii in 1:K){
      ki <- K*p+1+ii-1 # Number of theta in equation ii

      if (ii > 1) {
        X <- cbind(X2, -shortY[,1:(ii-1)])
      } else {
        X <- X2
      }
      if (is_tv[ii] == 1){
        thetai <- cbind(beta[, count_seqb[[ii]] ], alp[,count_seqa[[ii]] ])

        invvol <- as.vector(exp(-h[,ii]/2)) / sqrt_w[,ii]
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

        invvol <- as.vector(exp(-h[,ii]/2)) / sqrt_w[,ii]
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
    # matrix(beta0, nrow = K, byrow = T); alp0

    # sample alp and beta - equation by equation
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

        # % S = sparse(i,j,s,m,n,nzmax) uses vectors i, j, and s to generate an
        # %   m-by-n sparse matrix such that S(i(k),j(k)) = s(k), with space
        # %   allocated for nzmax nonzeros.

        Tthetai <- Matrix::Diagonal(t_max*ki) - Matrix::sparseMatrix( i = (ki+1):(t_max*ki),
                                                                      j = 1:((t_max-1)*ki),
                                                                      x = rep(1, (t_max-1)*ki),
                                                                      dims = c(t_max*ki,t_max*ki))

        thetai0 <- rep(0,ki)

        x.tilde = Matrix::t(Matrix::t(bigXi) %*% Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-0.5*h[,ii]) / sqrt_w[,ii]))
        TiST_thetai <- Matrix::t(Tthetai) %*% Tthetai

        y.tilde = Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-0.5*h[,ii]) / sqrt_w[,ii]) %*% U2[,ii]
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
    E = U / sqrt_w
    aux <- sample_h_ele(ytilde = t(E), sigma_h = Sigh, h0_mean = ah0, h = t(h), K = K, t_max = t_max, prior = priors)
    h <- t(aux$Sigtdraw)
    h0 <- as.numeric(aux$h0)
    Sigh <- diag(aux$sigma_h)


    # Sample w
    w_post_a <- reprow((nu+1)*0.5, t_max)
    w_post_b <- reprow(nu*0.5, t_max) + 0.5*U^2/exp(h)
    w <- matrix(mapply(rinvgamma, n = 1, shape = w_post_a, rate = w_post_b), ncol = K)
    sqrt_w = sqrt(w)
    # Sample nu
    nu_temp = nu + exp(logsigma_nu)*rnorm(K)
    for (k in c(1:K)){
      if (nu_temp[k] > 4 && nu_temp[k] < 100){
        num_mh = dgamma(nu_temp[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[,k], shape = nu_temp[k]*0.5, rate = nu_temp[k]*0.5, log = T))
        denum_mh = dgamma(nu[k], shape = nu_gam_a, rate = nu_gam_b, log = T) +
          sum(dinvgamma(w[,k], shape = nu[k]*0.5, rate = nu[k]*0.5, log = T))
        alpha = num_mh - denum_mh;
        temp = log(runif(1));
        if (alpha > temp){
          nu[k] = nu_temp[k]
          acount_nu[k] = acount_nu[k] + 1
        }

      }
    }
    if(isim %% batchlength == 0 ){
      for (jj in c(1:K)) {
        if (acount_nu[jj] > batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] + adaptamount(isim %/% batchlength);
        }
        if (acount_nu[jj] < batchlength * TARGACCEPT){
          logsigma_nu[jj] = logsigma_nu[jj] - adaptamount(isim %/% batchlength);
        }
        acount_nu[jj] = 0
      }
    }

    if ((isim > burnin) & (isim %% thin == 0) ) {
      isave <- (isim - burnin) %/% thin
      store_h[isave,,] <- h
      store_alp[isave,,] <- alp
      store_beta[isave,,] <- beta
      store_Sigbeta[isave,] <- Sigbeta
      store_Sigalp[isave,] <- Sigalp
      store_Sigh[isave,] <- Sigh
      store_beta0[isave,] <- beta0
      store_alp0[isave,] <- alp0
      store_h0[isave,] <- h0
      store_nu[isave,] <- nu
      store_w[isave,,] <- w
    }

    if ( isim %% 1000 == 0){
      cat(' Iteration ... ', isim, ' ', round(nu,1), '\n')
    }

  }
  end_time <- Sys.time()
  print( end_time - start_time)



  output <- list(store_h = store_h, store_alp = store_alp, store_beta = store_beta,
                 store_Sigbeta = store_Sigbeta, store_Sigalp = store_Sigalp, store_Sigh = store_Sigh,
                 store_beta0 = store_beta0, store_alp0 = store_alp0, store_h0 = store_h0,
                 store_nu = store_nu, store_w = store_w,
                 data = list(y = y, y0 = y0, p = p, priors = priors, inits = inits, dist = "Student"),
                 esttime = end_time - start_time,
                 class = "TVPStudentSV")
  return(output)
}
