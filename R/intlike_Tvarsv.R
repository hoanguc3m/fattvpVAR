#' @export
intlike_Tvarsv <- function(Yi,thetai0,Sig_hi,bigXi,h0i,nui){
  max_loop = 100
  K = length(Sig_hi)
  t_max = length(Yi)/K;
  # obtain the proposal density
  Hh = sparseMatrix(i = 1:(t_max*K),
                    j = 1:(t_max*K),
                    x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))

  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = rep(1/Sig_hi, t_max))
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh

  alph = Matrix::solve(Hh, Matrix::sparseMatrix(i = 1:K, j = rep(1,K), x = h0i, dims = c(t_max*K,1)))
  e = Yi - bigXi %*% thetai0  # Uni student
  s2 = e^2
  ht = log(s2 + 0.001)# alph + .01*rnorm(t_max)

  nu_vec <- rep(nui, t_max)
  errh_out = 1;
  while (errh_out> 0.001){
    # E-step
    s2_mod <- s2 * exp(-ht)  # Univariate student (y^2 exp(-h))
    Eilam = (nu_vec+1)/(nu_vec + s2_mod)
    s2Eilam = s2*Eilam
    # M-step
    htt = ht
    errh_in = 1
    while (errh_in> 0.001){
      eht = exp(htt)
      sieht = s2Eilam/eht
      fh = -.5 + .5*sieht
      Gh = .5*sieht
      Kh = HinvSH_h + Matrix::sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = as.numeric(Gh))
      newht = Matrix::solve(Kh, fh+Gh*htt+ HinvSH_h %*% alph)
      errh_in = max(abs(newht-htt))
      htt = newht
    }

    errh_out = max(abs(ht-htt))
    ht = htt
  }

  # compute negative Hessian
  Gh = (nu_vec+1)/(2*nu_vec) * (s2*exp(ht)) / ((exp(ht) + s2/nu_vec )^2)
  Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
  CKh = Matrix::t(Matrix::chol(Kh))

  # evaluate the importance weights
  c_pri = -t_max*K/2*log(2*pi) -.5*t_max*sum(log(Sig_hi))
  c_IS = -t_max*K/2*log(2*pi) + sum(log( Matrix::diag(CKh)))


  R = 50
  store_llike = rep(0, R)
  for (i in c(1:R)){
        hc = as.numeric(ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K)))

        llike = sum(dt(e*exp(-0.5*hc), df = nui, log = T)) - 0.5 * sum(hc)
        #sum(dnorm(e*exp(-0.5*hc), log = T)) - 0.5 * sum(hc)
        store_llike[i] = as.numeric(llike + c_pri - 0.5*Matrix::t(hc-alph) %*% HinvSH_h %*% (hc-alph) -
                                      (c_IS - 0.5*Matrix::t(hc-ht) %*% Kh %*% (hc-ht) ))
  }
  # increase simulation size if the variance of the log-likelihood > 1
  var_llike = var(store_llike)/R
  maxllike = max(store_llike)
  intlike = log(mean(exp(store_llike-maxllike))) + maxllike
  return(intlike)
}

#' @export
Chib_intlike_Tvarsv <- function(Yi,thetai0,Sig_hi,bigXi,h0i,nui){
  max_loop = 100
  K = length(Sig_hi)
  t_max = length(Yi)/K;
  prior_var_h0 <- 4
  # obtain the proposal density
  Hh = sparseMatrix(i = 1:(t_max*K),
                    j = 1:(t_max*K),
                    x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))

  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = c(1./ (Sig_hi + prior_var_h0), rep(1./Sig_hi, t_max-1))) # Prior for h1 \sim N(2 log sigma, sigmah^2 + 4 )
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh

  alph = Matrix::solve(Hh, Matrix::sparseMatrix(i = 1:K, j = rep(1,K), x = h0i, dims = c(t_max*K,1)))
  e = Yi - bigXi %*% thetai0  # Uni student
  s2 = e^2
  ht = log(s2 + 0.001)# alph + .01*rnorm(t_max)

  nu_vec <- rep(nui, t_max)
  errh_out = 1;
  while (errh_out> 0.001){
    # E-step
    s2_mod <- s2 * exp(-ht)  # Univariate student (y^2 exp(-h))
    Eilam = (nu_vec+1)/(nu_vec + s2_mod)
    s2Eilam = s2*Eilam
    # M-step
    htt = ht
    errh_in = 1
    while (errh_in> 0.001){
      eht = exp(htt)
      sieht = s2Eilam/eht
      fh = -.5 + .5*sieht
      Gh = .5*sieht
      Kh = HinvSH_h + Matrix::sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = as.numeric(Gh))
      newht = Matrix::solve(Kh, fh+Gh*htt+ HinvSH_h %*% alph)
      errh_in = max(abs(newht-htt))
      htt = newht
    }

    errh_out = max(abs(ht-htt))
    ht = htt
  }

  # compute negative Hessian
  Gh = (nu_vec+1)/(2*nu_vec) * (s2*exp(ht)) / ((exp(ht) + s2/nu_vec )^2)
  Kh = HinvSH_h + sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = Gh)
  CKh = Matrix::t(Matrix::chol(Kh))

  # evaluate the importance weights
  c_pri = -t_max*K/2*log(2*pi) -.5*(t_max-1)*sum(log(Sig_hi)) -.5*sum(log(Sig_hi + prior_var_h0))
  c_IS = -t_max*K/2*log(2*pi) + sum(log( Matrix::diag(CKh)))


  R = 50
  store_llike = rep(0, R)
  for (i in c(1:R)){
    hc = as.numeric(ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K)))

    llike = sum(dt(e*exp(-0.5*hc), df = nui, log = T)) - 0.5 * sum(hc)
    #sum(dnorm(e*exp(-0.5*hc), log = T)) - 0.5 * sum(hc)
    store_llike[i] = as.numeric(llike + c_pri - 0.5*Matrix::t(hc-alph) %*% HinvSH_h %*% (hc-alph) -
                                  (c_IS - 0.5*Matrix::t(hc-ht) %*% Kh %*% (hc-ht) ))
  }
  # increase simulation size if the variance of the log-likelihood > 1
  var_llike = var(store_llike)/R
  maxllike = max(store_llike)
  intlike = log(mean(exp(store_llike-maxllike))) + maxllike
  return(intlike)
}
