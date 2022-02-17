#' @export
denTy_h <- function(Yi,hc,Sigthetai,bigXi,thetai0,ws){

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

  invSig = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = exp(-hc) / ws )


  invStheta = sparseMatrix(i = 1:(t_max*ktheta),
                           j = 1:(t_max*ktheta),
                           x = rep(1/Sigthetai,t_max),
                           dims = c(t_max*ktheta,t_max*ktheta))

  XinvSig = Matrix::t(bigXi) %*% invSig
  HinvSH = Matrix::t(Htheta) %*% invStheta %*% Htheta;
  alptheta = Matrix::solve(Htheta, matrix(c(thetai0, rep(0, (t_max-1)*ktheta)), ncol = 1))
  Ktheta = HinvSH + XinvSig %*% bigXi;
  dtheta = XinvSig %*% Yi + HinvSH %*% alptheta;

  llike = -t_max*K*0.5*log(2*pi) - t_max/2*sum(log(Sigthetai)) - .5*sum(hc) - 0.5 * sum(log(ws)) -
    #sum(log(Matrix::diag(Matrix::chol(Ktheta)))) - .5*( Matrix::t(Yi) %*% invSig %*% Yi + Matrix::t(alptheta) %*% HinvSH %*% alptheta -
    sum(log(Matrix::diag(Matrix::chol(Ktheta)))) - .5*( sum(Yi^2 * exp(-hc) / ws) + Matrix::t(alptheta) %*% HinvSH %*% alptheta -
                                                          Matrix::t(dtheta) %*% Matrix::solve(Ktheta, dtheta))
  return(as.numeric(llike))
}

#' @export
intlike_Ttvpsv <- function(Yi, bigXi, Sigthetai, Sig_hi, h0i, thetai0, nui){

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
  e_w = 1; wt = rep(1,t_max); logwt = rep(0, t_max);
  countout = 0;
  while ( (e_h> .01 || e_w > 0.01) & countout < 100){

    # E-step
    invSig = sparseMatrix(i = 1:(t_max*K),j = 1:(t_max*K), x = exp(-ht) / wt )
    XinvSig = Matrix::t(bigXi) %*% invSig
    HinvSH_theta = Matrix::t(Htheta) %*% invStheta %*% Htheta;
    Ktheta = HinvSH_theta + XinvSig %*% bigXi;
    dtheta = XinvSig %*% Yi + HinvSH_theta %*% alptheta;
    # thetahat = Matrix::solve(Ktheta,dtheta)
    # CKtheta = Matrix::chol(Ktheta)
    # zhat = apply((bigXi %*% Matrix::solve(CKtheta) )^2, MARGIN = 1, FUN = sum) + (Yi-bigXi%*%thetahat)^2;
    # more robust and slightly faster - Sune Karlsson
    qq = bigXi %*% Matrix::solve(Ktheta, cbind(Matrix::t(bigXi), dtheta) );
    zhat = Matrix::diag(qq) + (Yi-qq[,t_max+1])^2;

    # Q_hwhw = - 0.5 * t((ht-alph)) %*% HinvSH_h %*% (ht-alph) - (nui *0.5 + 1) * sum(logwt) - sum(exp(-logwt) * nui * 0.5) -
    #   0.5 * ( sum(ht) + sum(logwt) ) - 0.5 * sum(exp(-ht) * exp(-logwt) * zhat)



    # M-step
    e_hj = 1; htt = ht; countin = 0;
    e_wj = 1; wtt = wt; logwtt = logwt;
    while ( (e_hj> .01 || e_wj > 0.1) & countin < 1000){

      einvhttzhat = exp(-htt)/wtt*zhat;
      gQ = -HinvSH_h %*% (htt-alph) -.5*(1-einvhttzhat);
      #gW = matrix(-(nui*0.5+1) + 0.5*nui/wtt - 0.5*(1-einvhttzhat))
      #gW = matrix(-(nui*0.5-0.5) + 0.5*nui/wtt - 0.5*(1-einvhttzhat))
      gW = matrix(-(nui*0.5-0.5) + 0.5*(nui+0.0)/wtt - 0.5*(1-einvhttzhat))



      gHW = rbind(gQ, gW)

      HHW = kronecker(matrix(c(1,0,0,0), ncol = 2, nrow = 2), -HinvSH_h) + sparseMatrix(i = (t_max+1):(t_max*2),j = (t_max+1):(t_max*2), x = - 0.5*nui/wtt) +
        kronecker(matrix(1, ncol = 2, nrow = 2), Diagonal(n = t_max, x = - 0.5*einvhttzhat))

      # HQ = -HinvSH_h -.5 * sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = einvhttzhat)
      # newhtt = htt - Matrix::solve(HQ,gQ);
      gdHW = Matrix::solve(HHW,gHW);
      newhtt = htt - gdHW[1:t_max];

      e_hj = mean(abs(as.numeric(newhtt-htt)));

      # gwQ = - (nui*0.5+1) + 0.5*nui/wtt - 0.5*(1-einvhttzhat)
      # HwQ = - 0.5*nui/wtt - 0.5*einvhttzhat
      # newlogwtt = logwtt - gwQ/HwQ;
      newlogwtt = logwtt - gdHW[(t_max+1):(t_max*2)];
      newwtt = exp(newlogwtt)
      e_wj = mean(abs(as.numeric(newlogwtt-logwtt)));

      htt = newhtt;
      logwtt = newlogwtt; wtt = newwtt;

      countin = countin + 1;
    }
    if (countin < 1000){
      e_h = mean(abs(as.numeric(ht-htt)));
      ht = htt;

      e_w = mean(abs(as.numeric(logwtt-logwt)));
      logwt = logwtt; wt = wtt;
    }
    countout = countout + 1;
  }


  invSig = sparseMatrix(i = 1:(t_max*K),j = 1:(t_max*K), x = exp(-ht) / wt )
  XinvSig = Matrix::t(bigXi) %*% invSig
  Z = Matrix::t(XinvSig) %*% Matrix::solve(Ktheta,Matrix::t(bigXi))
  HH = kronecker(matrix(1, ncol = 2, nrow = 2), -.5*Matrix::t(Z)*(Diagonal(K*t_max)-Z));
  Kh = -(HHW+HH);
  Cg = Matrix::t(Matrix::chol(Kh))

  invCg = Matrix::solve(Matrix::t(Cg))
  hw_mean = c(ht,logwt)

  # evaluate the importance weights
  c_pri = -t_max*K/2*log(2*pi) -.5*t_max*sum(log(Sig_hi))
  c_IS = -t_max*K/2*log(2*pi)*2 + sum(log(Matrix::diag(Cg)))

  R = 50
  store_llike = rep(0, R)

  for (i in c(1:R)){
    hw_sim = hw_mean + invCg %*% rnorm(t_max*K*2)
    hc = hw_sim[1:t_max]
    logws = hw_sim[(t_max+1):(2*t_max)]
    ws = exp(logws)

    store_llike[i] = denTy_h(Yi = Yi, hc = hc, Sigthetai = Sigthetai,
                              bigXi = bigXi, thetai0 = thetai0, ws = ws) +
      (c_pri - 0.5*Matrix::t(hc-alph) %*% HinvSH_h %*% (hc-alph)) +
      sum(dinvgamma(ws, shape = nui*0.5, rate = nui*0.5, log = TRUE)) -
      (c_IS - 0.5*Matrix::t(hw_sim-hw_mean)%*%Kh%*%(hw_sim-hw_mean) - sum(logws)) # Jacobian
  }
  # increase simulation size if the variance of the log-likelihood > 1
  var_llike = var(store_llike)/R

  maxllike = max(store_llike);
  intlike = log(mean(exp(store_llike-maxllike))) + maxllike;intlike

  return(intlike)
}
