#############################################################################################
#' @title D2Q_Obs_Fisher_bet
#' @description Auxiliary function to compute the observed Fisher information matrix for the beta regression.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @return Hessian of the Q-function.

D2Q_Obs_Fisher_bet = function(theta,z,x,v,link.mean,link.precision){
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap]
  lam = theta[-c(1:nkap)]
  nlam = length(lam)
  link_mean = stats::make.link(link.mean)
  link_precision = stats::make.link(link.precision)
  mu = link_mean$linkinv(x%*%kap) # mean parameter.
  phi = link_precision$linkinv(v%*%lam) # phi precision parameter.

  dmudeta = link_mean$mu.eta(x%*%kap)
  dphideta = link_precision$mu.eta(v%*%lam)
  d2mu = d2mudeta2(link.mean,mu)
  d2phi = d2phideta2(link.precision,phi)

  auxKK1 = ( trigamma(mu*phi)+trigamma((1-mu)*phi) ) * phi^2 * dmudeta^2
  auxKK2 = phi*(log(z)-log(1-z)-digamma(mu*phi)+digamma((1-mu)*phi)) * d2mu
  KK = diag(c(auxKK1 - auxKK2))

  auxLL1 = ( (mu^2)*trigamma(mu*phi)+((1-mu)^2)*trigamma((1-mu)*phi)) * dphideta^2
  auxLL2 = ( mu*(log(z)-log(1-z))+digamma(phi)+log(1-z)-mu*digamma(mu*phi)-(1-mu)*digamma((1-mu)*phi)) * d2phi
  LL = diag(c(auxLL1 + auxLL2))

  auxKL = (-log(z)+log(1-z)+digamma(mu*phi)-digamma((1-mu)*phi)+mu*phi*trigamma(mu*phi) -(1-mu)*phi*trigamma((1-mu)*phi))*dmudeta*dphideta
  KL = diag(c(auxKL))

  D2QKappa = (t(x)%*%KK)%*%x
  D2QKL = (t(x)%*%KL)%*%v
  D2QKappa = cbind(D2QKappa, D2QKL)
  D2QLambda = (t(v)%*%LL)%*%v
  D2QLambda = cbind(t(D2QKL), D2QLambda)
  D2Q = rbind(D2QKappa, D2QLambda)

  return(D2Q)
}

#############################################################################################
#' @title DQ2_Obs_Fisher_bet
#' @description Auxiliary function to compute the observed Fisher information matrix for the beta regression.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @return matrix given by the conditional expectation of the gradient of the Q-function and its tranpose.

DQ2_Obs_Fisher_bet = function(theta,z,x,v,link.mean,link.precision){
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap]
  lam = theta[-c(1:nkap)]
  nlam = length(lam)
  link_mean = stats::make.link(link.mean)
  link_precision = stats::make.link(link.precision)
  mu = link_mean$linkinv(x%*%kap) # mean parameter.
  phi = link_precision$linkinv(v%*%lam) # phi precision parameter.

  dmudeta = link_mean$mu.eta(x%*%kap)
  dphideta = link_precision$mu.eta(v%*%lam)
  d2mu = d2mudeta2(link.mean,mu)
  d2phi = d2phideta2(link.precision,phi)

  grad1 = c(( (log(z)-log(1-z)-digamma(mu*phi)+digamma((1-mu)*phi) ) * phi )*dmudeta)
  grad2 = c(( mu*(log(z)-log(1-z)) +digamma(phi) +log(1-z)-mu*digamma(mu*phi)-(1-mu)*digamma((1-mu)*phi)) * dphideta)

  aux_LL = digamma(phi)+log(1-z)+mu*(log(z)-log(1-z))
  aux_LL = aux_LL-mu*digamma(mu*phi)-(1-mu)*digamma((1-mu)*phi)
  aux_LL = (trigamma(phi)+aux_LL^2)
  aux_LL = aux_LL * dphideta^2

  KK_temp = grad1 %*% t(grad1)
  DQ2Kappa = (t(x)%*%KK_temp)%*%x
  KL_temp = grad1 %*% t(grad2)
  DQKL = (t(x)%*%KL_temp)%*%v
  LL_temp = grad2 %*% t(grad2)
  diag(LL_temp) = c(aux_LL)
  DQ2Lambda = (t(v)%*%LL_temp)%*%v

  DQ21 = cbind(DQ2Kappa, DQKL)
  DQ22 = cbind(t(DQKL), DQ2Lambda)
  DQ2 = rbind(DQ21, DQ22)
  return(DQ2)
}

#############################################################################################
#' @title infmat_bet
#' @description Function to compute standard errors based on the Fisher information matrix for the beta regression.
#' This function can also provide the Fisher's information matrix.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @param information optionally, a logical parameter indicating whether the Fisher's information matrix should be returned
#' @return Vector of standard errors or Fisher's information matrix if the parameter 'information' is set to TRUE.

infmat_bet = function(theta,z,x,v,link.mean,link.precision,information=FALSE)
{
  d2.Q = D2Q_Obs_Fisher_bet(theta,z,x,v,link.mean,link.precision)
  dd.Q = DQ2_Obs_Fisher_bet(theta,z,x,v,link.mean,link.precision)
  aux = d2.Q-dd.Q # Fisher Information Matrix.
  if(information){
    out = aux
  } else{
    inv.aux = tryCatch(solve(aux), error = function(e) rep(NA,nrow(aux)))
    out = inv.aux
    if(is.matrix(inv.aux)){ out = sqrt(diag(inv.aux)) } # Standard error.
  }
  return(out)
}

#############################################################################################
#' @title Qf_bet
#' @description Q-function related to the beta model. This function is required in the Expectation-Maximization algorithm.
#' @param theta vector of parameters (all coefficients).
#' @param phiold previous value of the precision parameter (phi).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @return Scalar representing the output of this auxiliary function for the beta case.
Qf_bet = function(theta,phiold,z,x,v,link.mean,link.precision)
{
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap]
  lam = theta[-c(1:nkap)]
  link_mean = stats::make.link(link.mean)
  link_precision = stats::make.link(link.precision)
  mu = link_mean$linkinv(x%*%kap) # mean parameter.
  phi = link_precision$linkinv(v%*%lam) # phi precision parameter.
  #
  mu0phi = mu*phi
  mu1phi = (1-mu)*phi
  mu0phi[which(mu0phi <=0)] = 10^(-10)
  mu1phi[which(mu1phi <=0)] = 10^(-10)
  #
  out = phi*(mu*log(z/(1-z))+digamma(phiold)+log(1-z));
  out = out -lgamma(mu0phi)-lgamma(mu1phi);
  out = out -log(z*(1-z)) -digamma(phiold) -log(1-z) -phiold;
  return(sum(out))
}

#############################################################################################
#' @title gradtheta_bet
#' @description Function to calculate the gradient of the Q-function, which is required for optimization via \code{optim}.
#' This option is related to the beta regression.
#' @param theta vector of parameters (all coefficients).
#' @param phiold previous value of the precision parameter (phi).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @return Scalar representing the output of this auxiliary gradient function for the beta case.
gradtheta_bet = function(theta,phiold,z,x,v,link.mean,link.precision)
{
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap]
  lam = theta[-c(1:nkap)]
  nlam = length(lam)
  link_mean = stats::make.link(link.mean)
  link_precision = stats::make.link(link.precision)
  mu = link_mean$linkinv(x%*%kap) # mean parameter.
  phi = link_precision$linkinv(v%*%lam) # phi precision parameter.
  #
  dmudeta = link_mean$mu.eta(x%*%kap)
  dphideta = link_precision$mu.eta(v%*%lam)
  aux = (log(z)-log(1-z)-digamma(mu*phi)+digamma((1-mu)*phi) ) * phi
  aux = aux*dmudeta
  Ukap = t(x)%*%aux
  aux = mu*(log(z)-log(1-z)) +digamma(phiold)
  aux = aux +log(1-z)-mu*digamma(mu*phi)
  aux = aux -(1-mu)*digamma((1-mu)*phi)
  aux = aux*dphideta
  Ulam = t(v)%*%aux

  return(c(Ukap,Ulam))
}

#############################################################################################
#' @title EMrun_bet
#' @description Function to run the Expectation-Maximization algorithm for the beta regression.
#' @param kap initial values for the coefficients in kappa related to the mean parameter.
#' @param lam initial values for the coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param epsilon tolerance to control the convergence criterion.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @return Vector containing the estimates for kappa and lambda in the beta regression.
EMrun_bet = function(kap,lam,z,x,v,epsilon,link.mean,link.precision){
  n = length(z)
  nkap = ncol(x)
  nlam = ncol(v)
  link_mean = stats::make.link(link.mean)
  link_precision = stats::make.link(link.precision)
  mu = link_mean$linkinv(x%*%kap) # mean parameter.
  phi = link_precision$linkinv(v%*%lam) # phi precision parameter.
  theta = c(kap,lam)
  count = 0

  repeat{
    theta_r = theta
    kap = theta[1:nkap]; lam = theta[-c(1:nkap)]
    mu = link_mean$linkinv(x%*%kap); phi = link_precision$linkinv(v%*%lam)
    phi_r = phi
    ### E step ------------------------------
    ### M step ------------------------------
    M = stats::optim(
      par = theta,
      fn = Qf_bet,
      gr = gradtheta_bet,
      phiold=phi_r,
      z = z ,
      x = x,
      v = v,
      link.mean = link.mean,
      link.precision = link.precision,
      control = list(fnscale = -1),
      method = 'L-BFGS-B'
    )
    theta = M$par
    # Compute Q -----------------------------
    Q_r = Qf_bet(theta_r,phi_r,z,x,v,link.mean,link.precision)
    Q = M$value
    ### Convergence criterion ---------------
    term1 = sqrt(sum((theta-theta_r)^2))
    term2 = abs(Q - Q_r)
    ### -------------------------------------
    count = count+1
    if(max(term1,term2) < epsilon){ break }
    if(count >= 10000){epsilon = 10^(-3)}
    ### -------------------------------------
  }
  
  if(sum(phi <= 0) > 0){
    warning("one or more estimates of precision parameters were negative. Please, 
         consider using another link function for the precision parameter.")
  }
  
  if(sum(is.nan(phi)) > 0){
    warning("one or more estimates of precision parameters were not a number. Please, 
         consider using another link function for the precision parameter.")    
  }

  gphi = 1/(1+phi)
  out = list()
  out[[1]] = c(kap,lam)
  out[[2]] = cbind(mu,gphi)
  out[[3]] = count
  names(out) = c("coeff","mu_gphi","n_iter")
  return(out)
}


#############################################################################################
#' @title simdata_bet
#' @description Function to generate synthetic data from the beta regression.
#' @param kap coefficients kappa related to the mean parameter.
#' @param lam coefficients lambda related to the precision parameter.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param repetitions the number of random draws to be made.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @return a list of response vectors z (with 0 < z_i < 1).
#' @seealso
#' \code{\link{simdata_bes}}, \code{\link{dbessel}}, \code{\link{dbbtest}}
#' @examples
#' n = 100; x = cbind(rbinom(n, 1, 0.5), runif(n, -1, 1)); v = runif(n, -1, 1);
#' z = simdata_bet(kap = c(1, -1, 0.5), lam = c(0.5,- 0.5), x, v, repetitions = 1,
#' link.mean = "logit", link.precision = "log")
#' z = unlist(z)
#' hist(z, xlim = c(0, 1), prob = TRUE)
#' @export
simdata_bet <- function(kap,lam,x,v,repetitions=1,link.mean,link.precision)
{
  ncolx = ncol(x); if(is.null(ncolx)==TRUE){ncolx=1}
  ncolv = ncol(v); if(is.null(ncolv)==TRUE){ncolv=1}
  nkap = length(kap); nlam = length(lam)
  if((nkap-ncolx)==1){x = cbind(1,x)}
  if((nlam-ncolv)==1){v = cbind(1,v)}
  if(abs(nkap-ncolx)>1){stop("check dimension of kappa and x")}
  if(abs(nlam-ncolv)>1){stop("check dimension of lambda and v")}
  #
  link_mean = stats::make.link(link.mean)
  link_precision = stats::make.link(link.precision)
  mu = link_mean$linkinv(x%*%kap)
  phi = link_precision$linkinv(v%*%lam)
  s1 = mu*phi
  s2 = phi*(1-mu)
  n = length(s1)
  Z = rep(0,n)

  Z = lapply(1:repetitions, function(x){  Y =  stats::rbeta(n, shape1=s1, shape2=s2)
  Y[Y < 0.00001] = 0.00001
  Y[Y > 0.99999] = 0.99999
  return(Y) #
  })
  return(Z)
}


#############################################################################################
#' @title envelope_bet
#' @description Function to calculate envelopes based on residuals for the beta regression.
#' @param residual character indicating the type of residual ("pearson", "score" or "quantile").
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param nsim_env number of synthetic data sets to be generated.
#' @param n sample size.
#' @param prob confidence level of the envelope (number between 0 and 1).
#' @param epsilon tolerance parameter used in the Expectation-Maximization algorithm applied to the synthetic data.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @seealso
#' \code{\link{score_residual_bet}}, \code{\link{quantile_residual_bet}}, \code{\link{pred_accuracy_bet}}
#' @return Matrix with dimension 2 x n (1st row = upper bound, second row = lower bound).
envelope_bet <- function(residual,kap,lam,x,v,nsim_env,prob,n,epsilon,link.mean,link.precision)
{
  zsim = simdata_bet(kap,lam,x,v, nsim_env,link.mean,link.precision)
  Res = switch(residual,
               pearson = {
                 Res = pblapply(zsim, function(zs){est = EMrun_bet(kap,lam,z=zs,x,v,epsilon,link.mean,link.precision)$mu_gphi
                 musim = est[,1]
                 gpsim = est[,2]
                 (zs-musim)/(sqrt(gpsim*musim*(1-musim)))})
                 Res
               },
               score = {
                 nkap = length(kap)
                 Res = pblapply(zsim, function(zs){est = EMrun_bet(kap,lam,z=zs,x,v,epsilon,link.mean,link.precision)$coeff
                 kapsim = est[1:nkap]
                 lamsim = est[-(1:nkap)]
                 score_residual_bet(kapsim,lamsim,zs,x,v,link.mean,link.precision)
                 })
                 Res
               },
               quantile = {
                 nkap = length(kap)
                 Res = pblapply(zsim, function(zs){est = EMrun_bet(kap,lam,z=zs,x,v,epsilon,link.mean,link.precision)$coeff
                 kapsim = est[1:nkap]
                 lamsim = est[-(1:nkap)]
                 quantile_residual_bet(kapsim,lamsim,zs,x,v,link.mean,link.precision)
                 })
                 Res
               }
  )
  Res = t(matrix(unlist(Res), n, nsim_env))
  Res = t(apply(Res,1,sort))
  Res = apply(Res,2,sort)
  id1 = max(1,round(nsim_env*(1-prob)/2))
  id2 = round(nsim_env*(1+prob)/2)
  Env = rbind(Res[id2,],apply(Res,2,mean),Res[id1,])
  rownames(Env) = c("upper","mean","lower")
  return(Env)
}

#############################################################################################
#' @title pred_accuracy_bet
#' @description Function to calculate the Residual Sum of Squares for partitions (training and test sets) of
#' the data set. Residuals are calculated here based on the beta regression.
#' @param residual Character indicating the type of residual ("pearson", "score" or "quantile").
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param ntest number of observations in the test set for prediction.
#' @param predict number of partitions (training and test sets) to be evaluated.
#' @param epsilon tolerance parameter used in the Expectation-Maximization algorithm for the training data set.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @seealso
#' \code{\link{score_residual_bet}}, \code{\link{quantile_residual_bet}}, \code{\link{envelope_bet}}
#' @return Vector containing the RSS for each partition of the full data set.
pred_accuracy_bet = function(residual,kap,lam,z,x,v,ntest,predict,epsilon,link.mean,link.precision){
  n = length(z)
  nkap = ncol(x)
  nlam = ncol(v)
  link_mean = stats::make.link(link.mean)
  link_precision = stats::make.link(link.precision)
  RSS_pred = rep(0,predict)
  if(predict>=50){ bar = utils::txtProgressBar(min=0,max=predict,style=3) }
  for(i in 1:predict){
    id_pred = sample(1:n,ntest,replace=FALSE)
    ztr = z[-id_pred]
    xtr = as.matrix(x[-id_pred,])
    vtr = as.matrix(v[-id_pred,])
    zte = z[id_pred]
    xte = as.matrix(x[id_pred,])
    vte = as.matrix(v[id_pred,])
    EMtr = EMrun_bet(kap,lam,ztr,xtr,vtr,epsilon)
    kaptr = EMtr[[1]][1:nkap]
    lamtr = EMtr[[1]][-(1:nkap)]
    mupred = link_mean$linkinv(xte%*%kaptr)
    phipred = link_precision$linkinv(vte%*%lamtr)
    gphi_pred = 1/(1+phipred)
    if(residual=="pearson"){ respred = (zte-mupred)/sqrt(gphi_pred*mupred*(1-mupred)) }
    if(residual=="score"){ respred = score_residual_bet(kaptr,lamtr,zte,xte,vte,link.mean,link.precision) }
    if(residual=="quantile"){ respred = quantile_residual_bet(kaptr,lamtr,zte,xte,vte,link.mean,link.precision) }
    RSS_pred[i] = sum(respred^2)
    if(predict>=50){ utils::setTxtProgressBar(bar,i) }
  }
  return(RSS_pred)
}

#############################################################################################
#' @title score_residual_bet
#' @description Function to calculate the empirical score residuals based on the beta regression.
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param nsim_score number synthetic data sets (default = 100) to be generated as a support to estime mean and s.d. of log(z)-log(1-z).
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @seealso
#' \code{\link{quantile_residual_bet}}
#' @return Vector containing the score residuals.
score_residual_bet = function(kap,lam,z,x,v,nsim_score=100,link.mean,link.precision){
    n = length(z)
    zs = simdata_bet(kap,lam,x,v, nsim_score,link.mean,link.precision)
    zs = lapply(zs, function(x){ log(x) - log(1-x)})
    zs = t(matrix(unlist(zs), n, nsim_score))
    me_zs = apply(zs,2,mean)
    sd_zs = apply(zs,2,stats::sd)
    out = log(z)-log(1-z)
    out = (out-me_zs)/sd_zs
  return(out)
}

#############################################################################################
#' @title quantile_residual_bet
#' @description Function to calculate quantile residuals based on the beta regression. Details about this type of residual can be found in \emph{Pereira (2019)}.
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt", "inverse".
#' @references DOI:10.1080/03610918.2017.1381740 (\href{https://www.tandfonline.com/doi/abs/10.1080/03610918.2017.1381740}{Pereira; 2019})
#' @seealso
#' \code{\link{score_residual_bet}}
#' @return Vector containing the quantile residuals.
quantile_residual_bet = function(kap,lam,z,x,v,link.mean,link.precision){
  n = length(z)
  link_mean = stats::make.link(link.mean)
  link_precision = stats::make.link(link.precision)
  mu = link_mean$linkinv(x%*%kap)
  phi = link_precision$linkinv(v%*%lam)
  out = rep(0,n)
  s1 = mu*phi
  s2 = (1-mu)*phi
  aux = stats::pbeta(z, shape1 = s1, shape2=s2)
  out = stats::qnorm(aux, mean=0,sd=1)
  return(out)
}


