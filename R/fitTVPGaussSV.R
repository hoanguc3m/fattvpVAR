#' @export
fitTVPGaussSV <- function(y, y0, p, priors, inits){

  Y0 <- y0
  shortY <- y
  t_max <- nrow(shortY)
  K <- ncol(shortY)
  Y <- matrixcalc::vec(t(shortY))


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
  abeta0 <- matrix(0, nrow = k_beta, ncol = 1);
  Vbeta0 <- 10*matrix(1, nrow = k_beta, ncol = 1);
  aalp0 <- matrix(0, nrow = k_alp, ncol = 1);
  Valp0 <- 10*matrix(1, nrow = k_alp, ncol = 1);

  EstMdl1 <- arima(y[,1] ,order = c(p,0,0))
  EstMdl2 <- arima(y[,2] ,order = c(p,0,0))
  EstMdl3 <- arima(y[,3] ,order = c(p,0,0))

  ah0 <- c(log(EstMdl1$sigma2), log(EstMdl2$sigma2), log(EstMdl3$sigma2))
  Vh0 <- 4*matrix(1, nrow = K,ncol = 1)
  # EstMdl1 <- var(y[,1])
  # EstMdl2 <- var(y[,2])
  # EstMdl3 <- var(y[,3])
  #
  # ah0 <- c(log(EstMdl1), log(EstMdl2), log(EstMdl3))
  # Vh0 <- 4*matrix(1, nrow = K,ncol = 1)

  hyper_ab <- priors$hyper_ab
  hyper_h <- priors$hyper_h

  # nuh0 <- 5*matrix(1, nrow = K, ncol = 1);
  # Sh0 <- hyper_h^2*matrix(1, nrow = K, ncol = 1)*(nuh0-1)
  Sh0 <- hyper_h*matrix(1, nrow = K, ncol = 1)

  # nubeta0 <- 5*matrix(1, nrow = k_beta, ncol = 1);
  # Sbeta0 <- hyper_ab^2*matrix(1, nrow = k_beta, ncol = 1)*(nubeta0-1)
  Sbeta0 <- hyper_ab*matrix(1, nrow = k_beta, ncol = 1)
  Sbeta0[seq(1,k_beta, by =  K*p+1)] <- 10*hyper_ab # intercepts more diffused prior 10*hyper_ab
  Sbeta0 <- Sbeta0*kronecker(is_tv, matrix(1, nrow = K*p+1, ncol = 1) )# set to 0 for time-invariant equations

  Salp0 <- hyper_ab*matrix(1, nrow = k_alp, ncol = 1)
  count <- 0
  for (ii in c(2:K)){
    if (is_tv[ii] == 0){
      Salp0[(count+1):(count+ii-1)] <- 0 # set to 0 for time-invariant equations
    }
    count <- count + ii - 1
  }


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
  for (j in 2:K){
    if (is_tv[j] == 1){
      idx_a_tv[ (count+1):(count+j-1)] <- TRUE
    }
    count <- count + j-1
  }

  # cpri <- -.5*(K+k_beta+k_alp)*log(2*pi) -.5*sum(log(Vbeta0)) -.5*sum(log(Vh0)) +
  #         sum(nubeta0[idx_b_tv]*log(Sbeta0[idx_b_tv])) - sum(lgamma(nubeta0[idx_b_tv])) +
  #         sum(nualp0[idx_a_tv]*log(Salp0[idx_a_tv])) - sum(lgamma(nualp0[idx_a_tv])) +
  #         sum(nuh0*log(Sh0)) - sum(lgamma(nuh0))
  # priorcalc <- function(sb,sa,sh,b0,a0,c0){
  #   return( cpri -.5*t(b0-abeta0) %*% ((b0-abeta0)/Vbeta0) -
  #           .5*t(a0-aalp0) %*% ((a0-aalp0)/Valp0) -.5*t((c0-ah0)/Vh0) %*% (c0-ah0) -
  #           t(nubeta0[idx_b_tv]+1) %*% log(sb[idx_b_tv]) -sum(Sbeta0[idx_b_tv]/sb[idx_b_tv]) -
  #           t(nualp0[idx_a_tv]+1) %*% log(sa[idx_a_tv]) -sum(Salp0[idx_a_tv]/sa[idx_a_tv]) -
  #           t(nuh0+1) %*% log(sh) - sum(Sh0/sh)
  #   )
  # }

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

  ## MCMC starts here
  set.seed(NULL)
  sprintf('Starting MCMC for (hybrid) Gaussian TVP-VAR-SV with ID: %s ', paste(as.character(is_tv), collapse = ' ') )
  start_time <- Sys.time()

  for (isim in c(1:(samples + burnin))){
    # sample alp and beta - equation by equation
    count <- 0
    U <- matrix(0, nrow = t_max, ncol = K)
    for (ii in 1:K){
      ki <- K*p+1+ii-1 # Number of theta in equation ii

      if (ii > 1) {
        X <- cbind(X2, -shortY[,1:(ii-1)])
      } else {
        X <- X2
      }

      if (is_tv[ii] == 1){

        bigXi <- SURform(X) # "dgCMatrix"
        # % S = sparse(i,j,s,m,n,nzmax) uses vectors i, j, and s to generate an
        # %   m-by-n sparse matrix such that S(i(k),j(k)) = s(k), with space
        # %   allocated for nzmax nonzeros.

        Tthetai <- Matrix::Diagonal(t_max*ki) - Matrix::sparseMatrix( i = (ki+1):(t_max*ki),
                                                                      j = 1:((t_max-1)*ki),
                                                                      x = rep(1, (t_max-1)*ki),
                                                                      dims = c(t_max*ki,t_max*ki))
        if (ii == 1){
          count_seq <- 0
        } else {
          count_seq <- (count+1):(count+ii-1)
        }

        Sigthetai <- c( Sigbeta[ ((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)], Sigalp[count_seq])
        thetai0 <- c( beta0[((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)], alp0[count_seq] )

        # XiSig <- Matrix::t(bigXi) %*% Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-h[,ii]))
        # TiST_thetai <- Matrix::t(Tthetai) %*% Matrix::sparseMatrix(i = 1:(t_max*ki), j = 1:(t_max*ki), x = rep(1./Sigthetai,t_max)) %*% Tthetai
        # Kthetai <- TiST_thetai + XiSig %*% bigXi
        #
        # thetai_hat <- Matrix::solve(Kthetai, TiST_thetai %*% Matrix(data = thetai0, nrow = t_max*length(thetai0), ncol = 1) + XiSig %*% shortY[,ii]  )
        # thetai <- thetai_hat + Matrix::solve(Matrix::chol(Kthetai), Matrix(rnorm(t_max*ki), ncol = 1) )

        x.tilde = Matrix::t(Matrix::t(bigXi) %*% Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-0.5*h[,ii])))
        TiST_thetai <- Matrix::t(Tthetai) %*% Matrix::sparseMatrix(i = 1:(t_max*ki), j = 1:(t_max*ki), x = rep(1./Sigthetai,t_max)) %*% Tthetai

        y.tilde = Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-0.5*h[,ii])) %*% shortY[,ii]
        theta.prior.precmean = TiST_thetai %*% Matrix(data = thetai0, nrow = t_max*length(thetai0), ncol = 1)
        theta.prec.chol <- Matrix::chol( TiST_thetai + Matrix::crossprod(x.tilde) )
        thetai <- Matrix::solve( theta.prec.chol,
                                 Matrix::solve( Matrix::t(theta.prec.chol), theta.prior.precmean + Matrix::crossprod( x.tilde, y.tilde ))
                             + rnorm(t_max*ki) )


        Thetai <- t(matrix(thetai, nrow = ki, ncol = t_max))

        beta[,((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)] <- Thetai[,1:k_beta_div_K]
        if ( (k_beta_div_K) < ncol(Thetai)){
          alp[,count_seq] <- Thetai[,(k_beta_div_K+1):ncol(Thetai)]
        }

      } else {
        if (ii == 1){
          count_seq <- 0
        } else {
          count_seq <- (count+1):(count+ii-1)
        }
        bigXi <- X
        # XiSig <- Matrix::t(bigXi) %*% Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-h[,ii]))
        # Vthetai <- c(Vbeta0[ ((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K),], Valp0[count_seq])
        #
        # thetai0 <- c(abeta0[((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)], aalp0[count_seq] )
        # Kthetai <- 1/Vthetai + XiSig %*% bigXi
        # Kthetai <- 0.5 * ( Kthetai + Matrix::t(Kthetai))
        # thetai_hat <- Matrix::solve(Kthetai, thetai0/Vthetai + XiSig %*% shortY[,ii])
        # thetai <- as.vector(thetai_hat + Matrix::solve(Matrix::chol(Kthetai) , Matrix(rnorm(ki), ncol = 1) ))

        x.tilde = Matrix::t(Matrix::t(bigXi) %*% Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-0.5*h[,ii])))
        Vthetai <- c(Vbeta0[ ((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K),], Valp0[count_seq])
        thetai0 <- c(abeta0[((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)], aalp0[count_seq] )

        y.tilde = Matrix::sparseMatrix(i = 1:t_max, j = 1:t_max, x = exp(-0.5*h[,ii])) %*% shortY[,ii]
        theta.prior.precmean = thetai0/Vthetai
        theta.prec.chol <- Matrix::chol( diag(1/Vthetai) + Matrix::crossprod(x.tilde) )
        thetai <- backsolve( theta.prec.chol,
                               backsolve( theta.prec.chol, theta.prior.precmean + Matrix::crossprod( x.tilde, y.tilde ),
                                          upper.tri = T, transpose = T )
                               + rnorm(ki) )

        betai <- thetai[1:k_beta_div_K];
        beta0[((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)] <- betai
        beta[,((ii-1)*k_beta_div_K+1):(ii*k_beta_div_K)] <- reprow(betai,t_max)
        if ( length(thetai) >= k_beta_div_K+1){
          alpi <- thetai[(k_beta_div_K+1):length(thetai)]
          alp0[count_seq] <- alpi
          alp[,count_seq] <- reprow(alpi,t_max)
        }


      }
      count <- count + ii-1
      U[,ii] <- as.matrix(shortY[,ii] - bigXi %*% thetai)
    }


    # sample h
    for (ii in c(1:K)){
      Ystar <- log(U[,ii]^2 + .0001)
      h[,ii] <- SVRW(Ystar = Ystar, h = h[,ii], sig = Sigh[ii], h0 = h0[ii])
    }

    # sample beta0
    if (sum(idx_b_tv) > 0){
      Kbeta0_tv <- Matrix::sparseMatrix(i = 1:(k_beta_div_K*n_tv),j = 1:(k_beta_div_K*n_tv),
                                        x = 1/Sigbeta[idx_b_tv] + 1/Vbeta0[idx_b_tv])

      beta0_tv_hat <- Matrix::solve(Kbeta0_tv, (abeta0[idx_b_tv]/Vbeta0[idx_b_tv] + beta[1,idx_b_tv]/Sigbeta[idx_b_tv]) )
      beta0[idx_b_tv] <- as.vector(beta0_tv_hat + Matrix::solve( Matrix::chol(Kbeta0_tv), rnorm(k_beta_div_K*n_tv) ))
    }

    # sample alp0
    if (sum(idx_a_tv) > 0){
      Kalp0_tv <- Matrix::sparseMatrix(i = 1:sum(idx_a_tv), j = 1:sum(idx_a_tv),
                               x = 1/Sigalp[idx_a_tv] + 1/Valp0[idx_a_tv])
      alp0_tv_hat <- Matrix::solve(Kalp0_tv, (aalp0[idx_a_tv] / Valp0[idx_a_tv] + alp[1,idx_a_tv]/Sigalp[idx_a_tv]))
      alp0[idx_a_tv] <- as.vector(alp0_tv_hat + Matrix::solve( Matrix::chol(Kalp0_tv) , rnorm(sum(idx_a_tv))))
    }
    # sample h0
    Kh0 <- Matrix::sparseMatrix(i = 1:K, j = 1:K, x = as.numeric(1/Sigh + 1/Vh0))
    h0_hat <- Matrix::solve(Kh0, (ah0/Vh0 + h[1,]/Sigh))
    h0 <- as.matrix(h0_hat + Matrix::solve(Matrix::chol(Kh0), rnorm(K)))

    # sample Sigbeta - InvGamma conjugate prior
    # E_beta <- beta[,idx_b_tv] - rbind(beta0[idx_b_tv], beta[1:(t_max-1),idx_b_tv])
    # Sigbeta[idx_b_tv] <- 1/ mapply(FUN = rgamma, n = 1, shape = t_max/2, rate = Sbeta0[idx_b_tv] + colSums(E_beta^2)/2)
    if (sum(idx_b_tv) > 0){
      Sigbeta[idx_b_tv] <- Sigma_sample(Beta = beta[,idx_b_tv, drop = FALSE],
                                      Beta0 = as.numeric(beta0[idx_b_tv]),
                                      Prior_Beta = Sbeta0[idx_b_tv],
                                      t_max = t_max)
    }
    # sample Sigalp - InvGamma conjugate prior
    # E_alp <- alp[,idx_a_tv] - rbind(alp0[idx_a_tv], alp[1:(t_max-1),idx_a_tv])
    # Sigalp[idx_a_tv] <- 1/ mapply(FUN = rgamma, n = 1, shape = nualp0[idx_a_tv]+t_max/2, rate = Salp0[idx_a_tv] + colSums(E_alp^2)/2)
    if (sum(idx_a_tv) > 0){
      Sigalp[idx_a_tv] <- Sigma_sample(Beta = alp[,idx_a_tv, drop = FALSE],
                                     Beta0 = as.numeric(alp0[idx_a_tv]),
                                     Prior_Beta = Salp0[idx_a_tv],
                                     t_max = t_max)

    }
    # sample Sigh
    # E_h <- h - rbind(t(h0), h[1:(t_max-1),])
    # Sigh <- 1/ mapply(FUN = rgamma, n = 1, shape = nuh0+t_max/2, rate = Sh0 + colSums(E_h^2)/2)
    Sigh <- Sigma_sample(Beta = h,
                         Beta0 = as.numeric(h0),
                         Prior_Beta = Sh0,
                         t_max = t_max)

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
    }

    if ( isim %% 1000 == 0){
      cat(' Iteration ... ', isim, '\n')
    }

  }
  end_time <- Sys.time()
  print( end_time - start_time)



  output <- list(store_h = store_h, store_alp = store_alp, store_beta = store_beta,
                 store_Sigbeta = store_Sigbeta, store_Sigalp = store_Sigalp, store_Sigh = store_Sigh,
                 store_beta0 = store_beta0, store_alp0 = store_alp0, store_h0 = store_h0,
                 data = list(y = y, y0 = y0, p = p, priors = priors, inits = inits, dist = "Gaussian"),
                 esttime = end_time - start_time,
                 class = "TVPGaussSV")
  return(output)
}

