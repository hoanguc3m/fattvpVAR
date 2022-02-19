#' @export
deny_h <- function(Yi,hc,Sigthetai,bigXi,thetai0){

  K = 1
  t_max = length(hc)

  ktheta = length(Sigthetai);
  Htheta = sparseMatrix(i = 1:(t_max*ktheta),
                        j = 1:(t_max*ktheta),
                        x = rep(1,t_max*ktheta)) -
    sparseMatrix( i = (ktheta+1):(t_max*ktheta),
                  j = 1:((t_max-1)*ktheta),
                  x = rep(1,(t_max-1)*ktheta),
                  dims =  c(t_max*ktheta, t_max*ktheta))

  invSig = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = exp(-hc))


  invStheta = sparseMatrix(i = 1:(t_max*ktheta),
                      j = 1:(t_max*ktheta),
                      x = rep(1/Sigthetai,t_max),
                      dims = c(t_max*ktheta,t_max*ktheta))

  XinvSig = Matrix::t(bigXi) %*% invSig
  HinvSH = Matrix::t(Htheta) %*% invStheta %*% Htheta;
  alptheta = Matrix::solve(Htheta, matrix(c(thetai0, rep(0, (t_max-1)*ktheta)), ncol = 1))
  Ktheta = HinvSH + XinvSig %*% bigXi;
  dtheta = XinvSig %*% Yi + HinvSH %*% alptheta;

  llike = -t_max*K*0.5*log(2*pi) - t_max/2*sum(log(Sigthetai)) - .5*sum(hc) -
   #sum(log(Matrix::diag(Matrix::chol(Ktheta)))) - .5*( Matrix::t(Yi) %*% invSig %*% Yi + Matrix::t(alptheta) %*% HinvSH %*% alptheta -
    sum(log(Matrix::diag(Matrix::chol(Ktheta)))) - .5*( sum(Yi^2 * exp(-hc)) + Matrix::t(alptheta) %*% HinvSH %*% alptheta -
  Matrix::t(dtheta) %*% Matrix::solve(Ktheta, dtheta))
  return(as.numeric(llike))
}

#' @export
intlike_tvpsv <- function(Yi,Sigthetai,Sig_hi,bigXi,h0i,thetai0){
  K = length(Sig_hi)
  t_max = length(Yi)/K;
  ktheta = length(Sigthetai)

  # obtain the mode of the marginal density of h (unconditional on theta)
  Htheta = sparseMatrix(i = 1:(t_max*ktheta),
                        j = 1:(t_max*ktheta),
                        x = rep(1,t_max*ktheta)) -
    sparseMatrix( i = (ktheta+1):(t_max*ktheta),
                  j = 1:((t_max-1)*ktheta),
                  x = rep(1,(t_max-1)*ktheta),
                  dims =  c(t_max*ktheta, t_max*ktheta))

  invStheta = sparseMatrix(i = 1:(t_max*ktheta), j = 1:(t_max*ktheta), x = rep(1/Sigthetai, t_max))
  alptheta = Matrix::solve(Htheta, matrix(c(thetai0, rep(0, (t_max-1)*ktheta))));

  Hh =  sparseMatrix(i = 1:(t_max*K),
                     j = 1:(t_max*K),
                     x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = rep(1./Sig_hi, t_max))
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0i, dims = c(t_max*K,1)))

  e_h = 1; ht = as.numeric(log(Yi^2+0.001))
  countout = 0;
  while ( e_h> .01 & countout < 100){
    # E-step
    invSig = sparseMatrix(i = 1:(t_max*K),j = 1:(t_max*K), x = exp(-ht))
    XinvSig = Matrix::t(bigXi) %*% invSig
    HinvSH_theta = Matrix::t(Htheta) %*% invStheta %*% Htheta;
    Ktheta = HinvSH_theta + XinvSig %*% bigXi;
    dtheta = XinvSig %*% Yi + HinvSH_theta %*% alptheta;
    # thetahat = Matrix::solve(Ktheta,dtheta)
    # CKtheta = Matrix::chol(Ktheta)
    # zhat = apply((bigXi %*% Matrix::solve(CKtheta) )^2, MARGIN = 1, FUN = sum) + (Yi-bigXi%*%thetahat)^2;
    # more robust and slightly faster - Sune Karlsson
    qq = bigXi %*% Matrix::solve(Ktheta, cbind(Matrix::t(bigXi), dtheta) );
    zhat = as.numeric( Matrix::diag(qq) + (Yi-qq[,t_max+1])^2)


    # M-step
    e_hj = 1; htt = ht; countin = 0;
    while ( e_hj> .01 & countin < 1000){

      einvhttzhat = exp(-htt)*zhat;
      gQ = -HinvSH_h %*% (htt-alph) -.5*(1-einvhttzhat);
      HQ = -HinvSH_h -.5 * Matrix::sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = einvhttzhat)
      newhtt = htt - Matrix::solve(HQ,gQ);
      e_hj = mean(abs(as.numeric(newhtt-htt)));
      htt = newhtt;
      countin = countin + 1;
      }
    if (countin < 1000){
        e_h = mean(abs(as.numeric(ht-htt)));
        ht = htt;
      }
    countout = countout + 1;
    }

    if (countout == 100){
      ht = rep(h0i,t_max);
      invSig = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = exp(-ht))
      XinvSig = Matrix::t(bigXi) %*% invSig;
      HinvSH_theta = Matrix::t(Htheta) %*% invStheta %*% Htheta;
      Ktheta = HinvSH_theta + XinvSig %*% bigXi;
      dtheta = XinvSig %*% Yi + HinvSH_theta %*% alptheta;
      # thetahat = Matrix::solve(Ktheta,dtheta)
      #CKtheta = Matrix::chol(Ktheta)
      #zhat = apply((CKtheta %*% solve(bigXi)^2),MARGIN = 2, FUN = sum) + (Yi-bigXi %*% thetahat)^2;
      qq = bigXi %*% Matrix::solve(Ktheta, cbind(Matrix::t(bigXi), dtheta) );
      zhat = Matrix::diag(qq) + (Yi-qq[,t_max+1])^2;

      einvhttzhat = exp(-ht)*zhat;
      HQ = -HinvSH_h -.5*sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = einvhttzhat)
    }

    Z = Matrix::t(XinvSig) %*% Matrix::solve(Ktheta,Matrix::t(bigXi))
    HH = -.5*Matrix::t(Z)*(Diagonal(K*t_max)-Z);
    Kh = -(HQ+HH);
    Cg = Matrix::t(Matrix::chol(Kh))

    # evaluate the importance weights
    c_pri = -t_max*K/2*log(2*pi) -.5*t_max*sum(log(Sig_hi))
    c_IS = -t_max*K/2*log(2*pi) + sum(log(Matrix::diag(Cg)))

    R = 20
    store_llike = rep(0, R)
    for (i in c(1:R)){
      hc = ht + Matrix::solve(Matrix::t(Cg), rnorm(t_max*K))
      #shorthc = matrix(hc, ncol = K, nrow = t_max);
            store_llike[i] = deny_h(Yi,hc,Sigthetai,bigXi,thetai0) + c_pri - 0.5*Matrix::t(hc-alph) %*% HinvSH_h %*% (hc-alph) -
              (c_IS - 0.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht) );
    }
    # increase simulation size if the variance of the log-likelihood > 1
    var_llike = var(store_llike)/R
    var_llike

    maxllike = max(store_llike);
    intlike = log(mean(exp(store_llike-maxllike))) + maxllike;
    return(intlike)
}

#' @export
Chib_intlike_tvpsv <- function(Yi,Sigthetai,Sig_hi,bigXi,h0i,thetai0){
  K = length(Sig_hi)
  t_max = length(Yi)/K;
  ktheta = length(Sigthetai)
  prior_var_h0 <- 4

  # obtain the mode of the marginal density of h (unconditional on theta)
  Htheta = sparseMatrix(i = 1:(t_max*ktheta),
                        j = 1:(t_max*ktheta),
                        x = rep(1,t_max*ktheta)) -
    sparseMatrix( i = (ktheta+1):(t_max*ktheta),
                  j = 1:((t_max-1)*ktheta),
                  x = rep(1,(t_max-1)*ktheta),
                  dims =  c(t_max*ktheta, t_max*ktheta))

  invStheta = sparseMatrix(i = 1:(t_max*ktheta), j = 1:(t_max*ktheta), x = rep(1/Sigthetai, t_max))
  alptheta = Matrix::solve(Htheta, matrix(c(thetai0, rep(0, (t_max-1)*ktheta))));

  Hh =  sparseMatrix(i = 1:(t_max*K),
                     j = 1:(t_max*K),
                     x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = c(1./ (Sig_hi + prior_var_h0), rep(1./Sig_hi, t_max-1))) # Prior for h1 \sim N(2 log sigma, sigmah^2 + 4 )

  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0i, dims = c(t_max*K,1)))

  e_h = 1; ht = as.numeric(log(Yi^2+0.001))
  countout = 0;
  while ( e_h> .01 & countout < 100){
    # E-step
    invSig = sparseMatrix(i = 1:(t_max*K),j = 1:(t_max*K), x = exp(-ht))
    XinvSig = Matrix::t(bigXi) %*% invSig
    HinvSH_theta = Matrix::t(Htheta) %*% invStheta %*% Htheta;
    Ktheta = HinvSH_theta + XinvSig %*% bigXi;
    dtheta = XinvSig %*% Yi + HinvSH_theta %*% alptheta;
    # thetahat = Matrix::solve(Ktheta,dtheta)
    # CKtheta = Matrix::chol(Ktheta)
    # zhat = apply((bigXi %*% Matrix::solve(CKtheta) )^2, MARGIN = 1, FUN = sum) + (Yi-bigXi%*%thetahat)^2;
    # more robust and slightly faster - Sune Karlsson
    qq = bigXi %*% Matrix::solve(Ktheta, cbind(Matrix::t(bigXi), dtheta) );
    zhat = as.numeric( Matrix::diag(qq) + (Yi-qq[,t_max+1])^2)


    # M-step
    e_hj = 1; htt = ht; countin = 0;
    while ( e_hj> .01 & countin < 1000){

      einvhttzhat = exp(-htt)*zhat;
      gQ = -HinvSH_h %*% (htt-alph) -.5*(1-einvhttzhat);
      HQ = -HinvSH_h -.5 * Matrix::sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = einvhttzhat)
      newhtt = htt - Matrix::solve(HQ,gQ);
      e_hj = mean(abs(as.numeric(newhtt-htt)));
      htt = newhtt;
      countin = countin + 1;
    }
    if (countin < 1000){
      e_h = mean(abs(as.numeric(ht-htt)));
      ht = htt;
    }
    countout = countout + 1;
  }

  if (countout == 100){
    ht = rep(h0i,t_max);
    invSig = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = exp(-ht))
    XinvSig = Matrix::t(bigXi) %*% invSig;
    HinvSH_theta = Matrix::t(Htheta) %*% invStheta %*% Htheta;
    Ktheta = HinvSH_theta + XinvSig %*% bigXi;
    dtheta = XinvSig %*% Yi + HinvSH_theta %*% alptheta;
    # thetahat = Matrix::solve(Ktheta,dtheta)
    #CKtheta = Matrix::chol(Ktheta)
    #zhat = apply((CKtheta %*% solve(bigXi)^2),MARGIN = 2, FUN = sum) + (Yi-bigXi %*% thetahat)^2;
    qq = bigXi %*% Matrix::solve(Ktheta, cbind(Matrix::t(bigXi), dtheta) );
    zhat = Matrix::diag(qq) + (Yi-qq[,t_max+1])^2;

    einvhttzhat = exp(-ht)*zhat;
    HQ = -HinvSH_h -.5*sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = einvhttzhat)
  }

  Z = Matrix::t(XinvSig) %*% Matrix::solve(Ktheta,Matrix::t(bigXi))
  HH = -.5*Matrix::t(Z)*(Diagonal(K*t_max)-Z);
  Kh = -(HQ+HH);
  Cg = Matrix::t(Matrix::chol(Kh))

  # evaluate the importance weights
  c_pri = -t_max*K/2*log(2*pi) -.5*(t_max-1)*sum(log(Sig_hi)) -.5*sum(log(Sig_hi + prior_var_h0))
  c_IS = -t_max*K/2*log(2*pi) + sum(log(Matrix::diag(Cg)))

  R = 20
  store_llike = rep(0, R)
  for (i in c(1:R)){
    hc = ht + Matrix::solve(Matrix::t(Cg), rnorm(t_max*K))
    #shorthc = matrix(hc, ncol = K, nrow = t_max);
    store_llike[i] = deny_h(Yi,hc,Sigthetai,bigXi,thetai0) + c_pri - 0.5*Matrix::t(hc-alph) %*% HinvSH_h %*% (hc-alph) -
      (c_IS - 0.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht) );
  }
  # increase simulation size if the variance of the log-likelihood > 1
  var_llike = var(store_llike)/R
  var_llike

  maxllike = max(store_llike);
  intlike = log(mean(exp(store_llike-maxllike))) + maxllike;
  return(intlike)
}
