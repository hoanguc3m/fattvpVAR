#' @export
deny_h <- function(Y,h,Sigtheta,bigX,theta0){

  K = ncol(h)
  t_max = nrow(h)

  k = length(Sigtheta);
  Htheta = sparseMatrix(i = 1:(t_max*K),
                        j = 1:(t_max*K),
                        x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))

  invSig = sparse(1:T*n,1:T*n,as.numeric(exp(-h)))

  invS = sparseMatrix(i = 1:(t_max*K),
                      j = 1:(t_max*K),
                      rep(1/Sigtheta,t_max),
                      dims = c(t_max*K,t_max*K))
  XinvSig = Matrix::t(bigX) %*% invSig
  HinvSH = Matrix::t(Htheta) %*% invS %*% Htheta;
  alptheta = Matrix::solve(Htheta, matrix(c(theta0, rep(0, (t_max-1)*K)), ncol = 1))
  Ktheta = HinvSH + XinvSig %*% bigX;
  dtheta = XinvSig %*% Y + HinvSH %*% alptheta;

  llike = -t_max*K/2*log(2*pi) - t_max/2*sum(log(Sigtheta)) - .5*sum(sum(h)) -
  - sum(log(diag(chol(Ktheta)))) - .5*(Matrix::t(Y) %*% invSig %*% Y + Matrix::t(alptheta) %*% HinvSH %*% alptheta -
  Matrix::t(dtheta) %*% Matrix::solve(Ktheta, dtheta))
  return(llike)
}

#' @export
intlike_tvpsv <- function(Y,Sigtheta,Sigh,bigX,h0,theta0,R){
  K = size(Sigh,1);
  t_max = size(Y,1)/K;
  k = size(Sigtheta,1);
  m = K*(K-1)/2;
  # obtain the mode of the marginal density of h (unconditional on theta)
  Htheta = sparseMatrix(i = 1:(t_max*K),
                        j = 1:(t_max*K),
                        x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))

  invS = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = rep(1./sigma_h, t_max))

  alptheta = Matrix::solve(Htheta, matrix(theta0, rep(1, (t_max-1)*K)));
  Hh =  sparseMatrix(i = 1:(t_max*K),
                     j = 1:(t_max*K),
                     x = rep(1,t_max*K)) -
    sparseMatrix( i = (K+1):(t_max*K),
                  j = 1:((t_max-1)*K),
                  x = rep(1,(t_max-1)*K),
                  dims =  c(t_max*K, t_max*K))
  SH = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = rep(1./sigma_h, t_max))
  HinvSH_h = Matrix::t(Hh) %*% SH %*% Hh
  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0, dims = c(t_max*K,1)))

  e_h = 1; ht = repmat(h0,t_max,1);
  countout = 0;
  while ( e_h> .01 & count < max_loop){
    # E-step
    invSig = sparseMatrix(i = 1:(t_max*K),j = 1:(t_max*K), x = exp(-ht))
    XinvSig = Matrix::t(bigX) %*% invSig
    HinvSH_theta = Matrix::t(Htheta) %*% invS %*% Htheta;
    Ktheta = HinvSH_theta + XinvSig %*% bigX;
    dtheta = XinvSig %*% Y + HinvSH_theta %*% alptheta;
    thetahat = Matrix::solve(Ktheta,dtheta)
    CKtheta = Matrix::t(Matrix::chol(Ktheta))
    zhat = apply((CKtheta %*% Matrix::solve(bigX) )^2, MARGIN = 2, FUN = sum) + (Y-bigX%*%thetahat)^2;

    # M-step
    e_hj = 1; htt = ht; countin = 0;
    while ( e_hj> .01 & countin < 1000){

      einvhttzhat = exp(-htt)*zhat;
      gQ = -HinvSH_h %*%(htt-alph) -.5*(1-einvhttzhat);
      HQ = -HinvSH_h -.5%*%sparse(1:t_max*K,1:t_max*K,einvhttzhat);
      newhtt = htt - Matrix::solve(HQ,gQ);
      e_hj = mean(abs(newhtt-htt));
      htt = newhtt;
      countin = countin + 1;
      }
      if (countin < 1000){
        e_h = mean(abs(ht-htt));
        ht = htt;
      }

      countout = countout + 1;
    }

    if (countout == 100){
      ht = repmat(h0,t_max,1);
      invSig = sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), exp(-ht))
      XinvSig = Matrix::t(bigX) %*% invSig;
      HinvSH_theta = Matrix::t(Htheta) %*% invS %*% Htheta;
      Ktheta = HinvSH_theta + XinvSig %*% bigX;
      dtheta = XinvSig %*% Y + HinvSH_theta %*% alptheta;
      thetahat = Matrix::solve(Ktheta,dtheta)
      CKtheta = Matrix::t(Matrix::chol(Ktheta))

      zhat = apply((CKtheta %*% solve(bigX)^2),MARGIN = 2, FUN = sum) + (Y-bigX %*% thetahat)^2;
      einvhttzhat = exp(-ht)*zhat;
      HQ = -HinvSH_h -.5*sparseMatrix(i = 1:(t_max*m), j = 1:(t_max*K),einvhttzhat)
    }
    Z = Matrix::t(XinvSig) %*% Matrix::solve(Ktheta,Matrix::t(bigX))
    HH = -.5*Matrix::t(Z)*(Diagonal(K*t_max)-Z);
    Kh = -(HQ+HH);
    Cg = Matrix::t(Matrix::chol(Kh))

    # evaluate the importance weights
    c_pri = -t_max*K/2*log(2*pi) -.5*t_max*sum(log(Sigh))
    c_IS = -t_max*K/2*log(2*pi) + sum(log(diag(Cg)))

    R = 10
    store_llike = rep(0, R)
    for (i in c(1:R)){
      hc = ht + Matrix::solve(Matrix::t(Cg), rnorm(t_max*K))
      llike = -t_max*K*0.5*log(2*pi) - 0.5*sum(hc) - 0.5*t(e) %*% sparse(1:t_max*K,1:t_max*K,exp(-hc)) %*% e

      #shorthc = matrix(hc, ncol = K, nrow = t_max);
            store_llike(i) = deny_h(Y,shorthc,Sigtheta,bigX,theta0) + c_pri -.5*Matrix::t(hc-alph) %*% HinvSH_h %*% (hc-alph) -
              (c_IS - -.5*Matrix::t(hc-ht)%*%Kh%*%(hc-ht) );
    }
    # increase simulation size if the variance of the log-likelihood > 1
    var_llike = var(store_llike)/R

    maxllike = max(store_llike);
    intlike = log(mean(exp(store_llike-maxllike))) + maxllike;
    return(intlike)
}
