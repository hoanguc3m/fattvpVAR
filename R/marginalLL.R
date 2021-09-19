#' @export
intlike_varsv <- function(Yi,thetai,Sig_hi,bigXi,h0i){
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

  alph = Matrix::solve(Hh, sparseMatrix(i = 1:K, j = rep(1,K), x = h0i, dims = c(t_max*K,1)))
  s2 = (Yi-bigXi %*% thetai)^2
  ht = alph + .01*rnorm(t_max)

  e_h = 1;
  count = 0;
  while ( e_h> .01 & count < max_loop){
    einvhts2 = exp(-ht)*s2
    gh = -HinvSH_h %*% (ht-alph) - 0.5*(1-einvhts2);
    Gh = -HinvSH_h -.5*sparseMatrix(i = 1:(t_max*K), j = 1:(t_max*K), x = einvhts2)
    newht = ht - Matrix::solve(Gh,gh)
    e_h = max(abs(newht-ht))
    ht = newht
    count = count + 1
  }

  if (count == max_loop){
    ht = rep(h0i,t_max)
    einvhts2 = exp(-ht)*s2
    Gh = -HinvSH_h -.5*sparseMatrix(i =1:(t_max*K), j = 1:(t_max*K), x = as.numeric(einvhts2))
  }

  Kh = -Gh;
  CKh = Matrix::t(Matrix::chol(Kh))

  # evaluate the importance weights
  c_pri = -t_max*K/2*log(2*pi) -.5*t_max*sum(log(Sig_hi))
  c_IS = -t_max*K/2*log(2*pi) + sum(log( Matrix::diag(CKh)))

  e = Yi - bigXi %*% thetai

  R = 15
  store_llike = rep(0, R)
  for (i in c(1:R)){
        hc = ht + Matrix::solve(Matrix::t(CKh), rnorm(t_max*K))
        llike = -t_max*K*0.5*log(2*pi) - 0.5*sum(hc) - 0.5*Matrix::t(e) %*%
          Matrix::sparseMatrix(i = 1:(t_max*K),
                               j = 1:(t_max*K),
                               x = exp(-hc)) %*% e
        store_llike[i] = as.numeric(llike + c_pri - 0.5*Matrix::t(hc-alph) %*% HinvSH_h %*% (hc-alph) -
                                      (c_IS - 0.5*Matrix::t(hc-ht) %*% Kh %*% (hc-ht) ))
  }
  # increase simulation size if the variance of the log-likelihood > 1
  var_llike = var(store_llike)/R
  maxllike = max(store_llike)
  intlike = log(mean(exp(store_llike-maxllike))) + maxllike
  return(intlike)
}
