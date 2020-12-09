#############################################################################################
#' @title D2Q_Obs_Fisher_bes
#' @description Auxiliary function to compute the observed Fisher information matrix for the bessel regression.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return Hessian of the Q-function.

D2Q_Obs_Fisher_bes <- function(theta, z, x, v, link.mean, link.precision) {
  n <- length(z)
  nkap <- ncol(x)
  kap <- theta[1:nkap]
  lam <- theta[-c(1:nkap)]
  nlam <- length(lam)
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  mu <- link_mean$linkinv(x %*% kap) # mean parameter.
  phi <- link_precision$linkinv(v %*% lam) # phi precision parameter.
  wz <- Ew1z(z, mu, phi)
  chi <- Ew2z(z, mu, phi)
  xi2 <- 1 + ((z - mu)^2) / (z * (1 - z))
  xit <- (z - mu) / (z * (1 - z))
  dmudeta <- link_mean$mu.eta(x %*% kap)
  dphideta <- link_precision$mu.eta(v %*% lam)
  d2mu <- d2mudeta2(link.mean, mu)
  d2phi <- d2phideta2(link.precision, phi)
  
  auxKK1 <- ((1 - 2 * mu * (1 - mu)) / (mu^2 * ((1 - mu)^2)) + wz * phi^2 / (z * (1 - z))) * dmudeta^2
  auxKK2 <- ((2 * mu - 1) / (mu * (1 - mu)) - wz * phi^2 * xit) * d2mu
  KK <- diag(c(auxKK1 + auxKK2))
  
  auxLL1 <- (2 / (phi^2) + wz * xi2) * dphideta^2
  auxLL2 <- (wz * phi * xi2 - 2 / phi - 1) * d2phi
  LL <- diag(c(auxLL1 + auxLL2))
  
  auxKL <- -2 * phi * wz * xit * dmudeta * dphideta
  KL <- diag(c(auxKL))
  
  D2QKappa <- (t(x) %*% KK) %*% x
  D2QKL <- (t(x) %*% KL) %*% v
  D2QKappa <- cbind(D2QKappa, D2QKL)
  D2QLambda <- (t(v) %*% LL) %*% v
  D2QLambda <- cbind(t(D2QKL), D2QLambda)
  D2Q <- rbind(D2QKappa, D2QLambda)
  
  return(D2Q)
}

#############################################################################################
#' @title DQ2_Obs_Fisher_bes
#' @description Auxiliary function to compute the observed Fisher information matrix for the bessel regression.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return matrix given by the conditional expectation of the gradient of the Q-function and its tranpose.

DQ2_Obs_Fisher_bes <- function(theta, z, x, v, link.mean, link.precision) {
  n <- length(z)
  nkap <- ncol(x)
  kap <- theta[1:nkap]
  lam <- theta[-c(1:nkap)]
  nlam <- length(lam)
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  mu <- link_mean$linkinv(x %*% kap) # mean parameter.
  phi <- link_precision$linkinv(v %*% lam) # phi precision parameter.
  wz <- Ew1z(z, mu, phi)
  chi <- Ew2z(z, mu, phi)
  xi2 <- 1 + ((z - mu)^2) / (z * (1 - z))
  xit <- (z - mu) / (z * (1 - z))
  dmudeta <- link_mean$mu.eta(x %*% kap)
  dphideta <- link_precision$mu.eta(v %*% lam)
  d2mu <- d2mudeta2(link.mean, mu)
  d2phi <- d2phideta2(link.precision, phi)
  
  grad1 <- c(((1 - 2 * mu) / (mu * (1 - mu)) + wz * (phi^2) * ((z - mu) / (z * (1 - z)))) * dmudeta)
  grad2 <- c(((2 / phi) + 1 - wz * phi * xi2) * dphideta)
  
  aux_KK <- (((1 - 2 * mu) / (mu * (1 - mu)))^2 + 2 * xit * wz * phi^2 * (1 - 2 * mu) / (z * (1 - z)) + chi * phi^4 * xit^2) * dmudeta^2
  
  aux_KL <- (((1 - 2 * mu) / (mu * (1 - mu))) * (2 / phi + 1 - wz * phi * xi2) + phi^2 * xit * (wz * (2 / phi + 1) - chi * phi * xi2)) * dmudeta * dphideta
  
  aux_LL <- ((2 / phi + 1)^2 - 2 * wz * phi * (2 / phi + 1) * xi2 + chi * phi^2 * xi2^2) * dphideta^2
  
  KK_temp <- grad1 %*% t(grad1)
  diag(KK_temp) <- c(aux_KK)
  DQ2Kappa <- (t(x) %*% KK_temp) %*% x
  KL_temp <- grad1 %*% t(grad2)
  diag(KL_temp) <- c(aux_KL)
  DQKL <- (t(x) %*% KL_temp) %*% v
  LL_temp <- grad2 %*% t(grad2)
  diag(LL_temp) <- c(aux_LL)
  DQ2Lambda <- (t(v) %*% LL_temp) %*% v
  
  DQ21 <- cbind(DQ2Kappa, DQKL)
  DQ22 <- cbind(t(DQKL), DQ2Lambda)
  DQ2 <- rbind(DQ21, DQ22)
  return(DQ2)
}

#############################################################################################
#' @title infmat_bes
#' @description Function to compute standard errors based on the Fisher information matrix for the bessel regression.
#' This function can also provide the Fisher's information matrix.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @param information optionally, a logical parameter indicating whether the Fisher's information matrix should be returned
#' @return Vector of standard errors or Fisher's information matrix if the parameter 'information' is set to TRUE.

infmat_bes <- function(theta, z, x, v, link.mean, link.precision, information = FALSE) {
  d2.Q <- D2Q_Obs_Fisher_bes(theta, z, x, v, link.mean, link.precision)
  dd.Q <- DQ2_Obs_Fisher_bes(theta, z, x, v, link.mean, link.precision)
  aux <- d2.Q - dd.Q # Fisher Information Matrix.
  if (information) {
    out <- aux
  } else {
    inv.aux <- tryCatch(solve(aux), error = function(e) rep(NA, nrow(aux)))
    out <- inv.aux
    if (is.matrix(inv.aux)) {
      out <- sqrt(diag(inv.aux))
    } # Standard error.
  }
  if(sum(is.nan(out))>0){
    warning("Please, try another optimization method.'")
  }
  return(out)
}

#############################################################################################
#' @title Ew1z
#' @description Auxiliary function required in the Expectation-Maximization algorithm (E-step) and in the calculation of the
#' Fisher information matrix. It represents the conditional expected value E(W_i^s|Z_i), with s = -1; i.e.,
#' latent W_i^(-1) given the observed Z_i.
#' @param z response vector with 0 < z_i < 1.
#' @param mu mean parameter (vector having the same size of z).
#' @param phi precision parameter (vector having the same size of z).
#' @return Vector of expected values.
Ew1z <- function(z, mu, phi) {
  xi <- sqrt(((mu^2) / z) + (((1 - mu)^2) / (1 - z)))
  prod <- phi * xi
  bK.num <- besselK(prod, -2, expon.scaled = TRUE)
  bK.den <- besselK(prod, -1, expon.scaled = TRUE)
  out <- (1 / prod) * (bK.num / bK.den)
  return(out)
}

#############################################################################################
#' @title Ew2z
#' @description Auxiliary function required in the calculation of the
#' Fisher information matrix. It represents the conditional expected value E(W_i^s|Z_i), with s = -2; i.e.,
#' latent W_i^(-2) given the observed Z_i.
#' @param z response vector with 0 < z_i < 1.
#' @param mu mean parameter (vector having the same size of z).
#' @param phi precision parameter (vector having the same size of z).
#' @return vector of expected values.
Ew2z <- function(z, mu, phi) {
  xi <- sqrt(((mu^2) / z) + (((1 - mu)^2) / (1 - z)))
  prod <- phi * xi
  prod2 <- prod^2
  bK.num <- besselK(prod, -3, expon.scaled = TRUE)
  bK.den <- besselK(prod, -1, expon.scaled = TRUE)
  out <- (1 / prod2) * (bK.num / bK.den)
  return(out)
}

#############################################################################################
#' @title Qf_bes
#' @description Q-function related to the bessel model. This function is required in the Expectation-Maximization algorithm.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param wz parameter representing E(1/W_i|Z_i = z_i, theta).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return Scalar representing the output of this auxiliary function for the bessel case.
Qf_bes <- function(theta, wz, z, x, v, link.mean, link.precision) {
  n <- length(z)
  nkap <- ncol(x)
  kap <- theta[1:nkap]
  lam <- theta[-c(1:nkap)]
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  mu <- link_mean$linkinv(x %*% kap) # mean parameter.
  phi <- link_precision$linkinv(v %*% lam) # phi precision parameter.
  #
  out1 <- log(mu) + log(1 - mu) + 2 * log(phi) + phi
  out2 <- 0.5 * wz * (phi^2) * (((mu^2) / z) + (((1 - mu)^2) / (1 - z)))
  
  return(sum(out1 - out2))
}

#############################################################################################
#' @title gradtheta_bes
#' @description Function to calculate the gradient of the Q-function, which is required for optimization via \code{optim}.
#' This option is related to the bessel regression.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param wz parameter representing E(1/W_i|Z_i = z_i, theta).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return vector representing the output of this auxiliary gradient function for the bessel case.
gradtheta_bes <- function(theta, wz, z, x, v, link.mean, link.precision) {
  n <- length(z)
  nkap <- ncol(x)
  kap <- theta[1:nkap]
  lam <- theta[-c(1:nkap)]
  nlam <- length(lam)
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  mu <- link_mean$linkinv(x %*% kap) # mean parameter.
  phi <- link_precision$linkinv(v %*% lam) # phi precision parameter.
  #
  dmudeta <- link_mean$mu.eta(x %*% kap)
  dphideta <- link_precision$mu.eta(v %*% lam)
  aux <- (1 - 2 * mu) / (mu * (1 - mu))
  aux <- aux + wz * (phi^2) * ((z - mu) / (z * (1 - z)))
  aux <- aux * dmudeta
  Ukap <- t(x) %*% aux
  aux <- (2 / phi) + 1 - wz * phi * (1 + ((z - mu)^2) / (z * (1 - z)))
  aux <- aux * dphideta
  Ulam <- t(v) %*% aux
  
  return(c(Ukap, Ulam))
}

#############################################################################################
#' @title EMrun_bes
#' @description Function to run the Expectation-Maximization algorithm for the bessel regression.
#' @param kap initial values for the coefficients in kappa related to the mean parameter.
#' @param lam initial values for the coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param epsilon tolerance to control the convergence criterion.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function. The default is set to "L-BFGS-B".
#' @return Vector containing the estimates for kappa and lambda in the bessel regression.
EMrun_bes <- function(kap, lam, z, x, v, epsilon, link.mean, link.precision, optim_method = "L-BFGS-B") {
  n <- length(z)
  nkap <- ncol(x)
  nlam <- ncol(v)
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  mu <- link_mean$linkinv(x %*% kap) # mean parameter.
  phi <- link_precision$linkinv(v %*% lam) # phi precision parameter.
  theta <- c(kap, lam)
  count <- 0
  
  repeat{
    theta_r <- theta
    kap <- theta[1:nkap]
    lam <- theta[-c(1:nkap)]
    mu <- link_mean$linkinv(x %*% kap)
    phi <- link_precision$linkinv(v %*% lam)
    ### E step ------------------------------
    wz_r <- Ew1z(z, mu, phi)
    ### M step ------------------------------
    M <- tryCatch(stats::optim(
      par = theta,
      fn = Qf_bes,
      gr = gradtheta_bes,
      wz = wz_r,
      z = z,
      x = x,
      v = v,
      link.mean = link.mean,
      link.precision = link.precision,
      control = list(fnscale = -1),
      method = optim_method
    ), error = function(e) {
      "Error"
    })
    
    if (length(M) == 1) {
      warning("Trying with numerical derivatives")
      M <- tryCatch(stats::optim(
        par = theta,
        fn = Qf_bes,
        gr = NULL,
        wz = wz_r,
        z = z,
        x = x,
        v = v,
        link.mean = link.mean,
        link.precision = link.precision,
        control = list(fnscale = -1),
        method = optim_method
      ), error = function(e) {
        "Error"
      })
    }
    
    if (length(M) == 1) {
      warning("Trying with another optimization algorithm")
      if(optim_method == "L-BFGS-B"){
        optim_temp = "Nelder-Mead"
      } else if(optim_method == "Nelder-Mead"){
        optim_temp = "L-BFGS-B"
      } else {
        optim_temp = "Nelder-Mead"
      }
      M <- tryCatch(stats::optim(
        par = theta,
        fn = Qf_bes,
        gr = NULL,
        wz = wz_r,
        z = z,
        x = x,
        v = v,
        link.mean = link.mean,
        link.precision = link.precision,
        control = list(fnscale = -1),
        method = optim_temp
      ), error = function(e) {
        "Error"
      })
    }
    
    if (length(M) == 1) {
      warning("The EM algorithm did not converge.")
      break
    }
    theta <- M$par
    # Compute Q -----------------------------
    Q_r <- Qf_bes(theta_r, wz_r, z, x, v, link.mean, link.precision)
    Q <- M$value
    ### Convergence criterion ---------------
    term1 <- sqrt(sum((theta - theta_r)^2))
    term2 <- abs(Q - Q_r)
    ### -------------------------------------
    count <- count + 1
    if (max(term1, term2) < epsilon) {
      break
    }
    if (count >= 10000) {
      epsilon <- 10^(-3)
    }
    ### -------------------------------------
  }
  
  if (sum(phi <= 0) > 0) {
    stop("one or more estimates of precision parameters were negative. Please, 
         use another link function for the precision parameter.")
  }
  
  gphi <- (1 - phi + (phi^2) * exp(phi) * expint_En(phi, order = 1)) / 2
  out <- list()
  out[[1]] <- c(kap, lam)
  out[[2]] <- cbind(mu, gphi)
  out[[3]] <- count
  names(out) <- c("coeff", "mu_gphi", "n_iter")
  return(out)
}

#############################################################################################
#' @title simdata_bes
#' @description Function to generate synthetic data from the bessel regression.
#' Requires the R package "statmod" generate random numbers from the Inverse Gaussian distribution (\emph{Giner and Smyth, 2016}).
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param repetitions the number of random draws to be made.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return a list of response vectors z (with 0 < z_i < 1).
#' @references DOI:10.32614/RJ-2016-024 (\href{https://journal.r-project.org/archive/2016/RJ-2016-024/index.html}{Giner and Smyth; 2016})
#' @seealso
#' \code{\link{dbessel}}, \code{\link{dbbtest}}, \code{\link{simdata_bet}}
#' @examples
#' n <- 100
#' x <- cbind(rbinom(n, 1, 0.5), runif(n, -1, 1))
#' v <- runif(n, -1, 1)
#' z <- simdata_bes(
#'   kap = c(1, -1, 0.5), lam = c(0.5, -0.5), x, v,
#'   repetitions = 1, link.mean = "logit", link.precision = "log"
#' )
#' z <- unlist(z)
#' hist(z, xlim = c(0, 1), prob = TRUE)
#' @export
simdata_bes <- function(kap, lam, x, v, repetitions = 1, link.mean, link.precision) {
  ncolx <- ncol(x)
  if (is.null(ncolx) == TRUE) {
    ncolx <- 1
  }
  ncolv <- ncol(v)
  if (is.null(ncolv) == TRUE) {
    ncolv <- 1
  }
  nkap <- length(kap)
  nlam <- length(lam)
  if ((nkap - ncolx) == 1) {
    x <- cbind(1, x)
  }
  if ((nlam - ncolv) == 1) {
    v <- cbind(1, v)
  }
  if (abs(nkap - ncolx) > 1) {
    stop("check dimension of kappa and x")
  }
  if (abs(nlam - ncolv) > 1) {
    stop("check dimension of lambda and v")
  }
  #
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  mu <- link_mean$linkinv(x %*% kap)
  phi <- link_precision$linkinv(v %*% lam)
  a <- mu * phi
  b <- phi * (1 - mu)
  n <- length(a)
  
  Z <- lapply(1:repetitions, function(x) {
    Y1 <- rinvgauss(n, mean = a, shape = a^2)
    Y2 <- rinvgauss(n, mean = b, shape = b^2)
    return(Y1 / (Y1 + Y2)) # Z[i] ~ Bessel(mu[i],phi[i])
  })
  return(Z)
}


#############################################################################################
#' @title envelope_bes
#' @description Function to calculate envelopes based on residuals for the bessel regression.
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
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return Matrix with dimension 2 x n (1st row = upper bound, second row = lower bound).
envelope_bes <- function(residual, kap, lam, x, v, nsim_env, prob, n, epsilon, link.mean, link.precision) {
  zsim <- simdata_bes(kap, lam, x, v, nsim_env, link.mean, link.precision)
  Res <- switch(residual,
                pearson = {
                  Res <- pblapply(zsim, function(zs) {
                    est <- EMrun_bes(kap, lam, z = zs, x, v, epsilon, link.mean, link.precision)$mu_gphi
                    musim <- est[, 1]
                    gpsim <- est[, 2]
                    (zs - musim) / (sqrt(gpsim * musim * (1 - musim)))
                  })
                  Res
                },
                score = {
                  nkap <- length(kap)
                  Res <- pblapply(zsim, function(zs) {
                    est <- EMrun_bes(kap, lam, z = zs, x, v, epsilon, link.mean, link.precision)$coeff
                    kapsim <- est[1:nkap]
                    lamsim <- est[-(1:nkap)]
                    score_residual_bes(kapsim, lamsim, zs, x, v, link.mean, link.precision)
                  })
                  Res
                },
                quantile = {
                  nkap <- length(kap)
                  Res <- pblapply(zsim, function(zs) {
                    est <- EMrun_bes(kap, lam, z = zs, x, v, epsilon, link.mean, link.precision)$coeff
                    kapsim <- est[1:nkap]
                    lamsim <- est[-(1:nkap)]
                    quantile_residual_bes(kapsim, lamsim, zs, x, v, link.mean, link.precision)
                  })
                  Res
                }
  )
  Res <- t(matrix(unlist(Res), n, nsim_env))
  Res <- t(apply(Res, 1, sort))
  Res <- apply(Res, 2, sort)
  id1 <- max(1, round(nsim_env * (1 - prob) / 2))
  id2 <- round(nsim_env * (1 + prob) / 2)
  Env <- rbind(Res[id2, ], apply(Res, 2, mean), Res[id1, ])
  rownames(Env) <- c("upper", "mean", "lower")
  return(Env)
}

#############################################################################################
#' @title pred_accuracy_bes
#' @description Function to calculate the Residual Sum of Squares for partitions (training and test sets) of
#' the data set. Residuals are calculated here based on the bessel regression.
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
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return Vector containing the RSS for each partition of the full data set.
pred_accuracy_bes <- function(residual, kap, lam, z, x, v, ntest, predict, epsilon, link.mean, link.precision) {
  n <- length(z)
  nkap <- ncol(x)
  nlam <- ncol(v)
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  RSS_pred <- rep(0, predict)
  if (predict >= 50) {
    bar <- utils::txtProgressBar(min = 0, max = predict, style = 3)
  }
  for (i in 1:predict) {
    id_pred <- sample(1:n, ntest, replace = FALSE)
    ztr <- z[-id_pred]
    xtr <- as.matrix(x[-id_pred, ])
    vtr <- as.matrix(v[-id_pred, ])
    zte <- z[id_pred]
    xte <- as.matrix(x[id_pred, ])
    vte <- as.matrix(v[id_pred, ])
    EMtr <- EMrun_bes(kap, lam, ztr, xtr, vtr, epsilon, link.mean, link.precision)
    kaptr <- EMtr[[1]][1:nkap]
    lamtr <- EMtr[[1]][-(1:nkap)]
    mupred <- link_mean$linkinv(xte %*% kaptr)
    phipred <- link_precision$linkinv(vte %*% lamtr)
    gphi_pred <- (1 - phipred + (phipred^2) * exp(phipred) * expint_En(phipred, order = 1)) / 2
    respred <- switch(residual,
                      pearson = {
                        (zte - mupred) / sqrt(gphi_pred * mupred * (1 - mupred))
                      },
                      score = {
                        score_residual_bes(kaptr, lamtr, zte, xte, vte, link.mean, link.precision)
                      },
                      quantile = {
                        quantile_residual_bes(kaptr, lamtr, zte, xte, vte, link.mean, link.precision)
                      }
    )
    RSS_pred[i] <- sum(respred^2)
    if (predict >= 50) {
      utils::setTxtProgressBar(bar, i)
    }
  }
  return(RSS_pred)
}

#############################################################################################
#' @title score_residual_bes
#' @description Function to calculate the empirical score residuals based on the bessel regression.
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param nsim_score number synthetic data sets (default = 100) to be generated as a support to estime mean and s.d. of log(z)-log(1-z).
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @seealso
#' \code{\link{quantile_residual_bes}}
#' @return Vector containing the score residuals.
score_residual_bes <- function(kap, lam, z, x, v, nsim_score = 100, link.mean, link.precision) {
  n <- length(z)
  zs <- simdata_bes(kap, lam, x, v, nsim_score, link.mean, link.precision)
  zs <- lapply(zs, function(x) {
    log(x) - log(1 - x)
  })
  zs <- t(matrix(unlist(zs), n, nsim_score))
  me_zs <- apply(zs, 2, mean)
  sd_zs <- apply(zs, 2, stats::sd)
  out <- log(z) - log(1 - z)
  out <- (out - me_zs) / sd_zs
  return(out)
}

#############################################################################################
#' @name bessel
#' @title Bessel Distribution
#' @aliases dbessel rbessel pbessel qbessel
#' @description Functions to calculate the cumulative distribution, probability density, 
#' generate bessel random numbers and find quantiles of the bessel distribution
#' @param z vector of numbers in (0,1) for which the p.d.f. is to be evaluated.
#' @param p vector of probabilities.
#' @param q vector of quantiles.
#' @param mu vector of numbers in (0,1) representing the mean parameter.
#' @param phi vector of positive numbers representing the precision parameter.
#' @param n sample size.
#' @return scalar expressing the value of the density at z.
#' @seealso
#' \code{\link{simdata_bes}}, \code{\link{dbbtest}}, \code{\link{simdata_bet}}
#' @examples
#' rbessel(100, mu = 0.2, phi = 10)
#' pbessel(0.4, mu = 0.1)
#' qbessel(0.8, mu = 0.8)
#' plot(dbessel)
#' @rdname bessel
#' @export
dbessel <- function(z, mu = 1/2, phi = 1) {
  dens_bessel <- sapply(z, function(x){
    if(x <= 0 | x >=1){
      return(0)
    } else{
      zeta <- sqrt((x * (1 - 2 * mu) + mu^2) / (x - x^2))
      out <- mu * (1 - mu) * phi * exp(phi) * besselK((phi * zeta), 1)
      out <- out / (zeta * pi * (x* (1 - x))^(3 / 2))
      return(out)
    }
  })
  return(dens_bessel)
}

#' @rdname bessel
#' @export
rbessel <- function(n, mu = 1/2, phi = 1){
  a <- mu * phi
  b <- phi * (1 - mu)
  Y1 <- rinvgauss(n, mean = a, shape = a^2)
  Y2 <- rinvgauss(n, mean = b, shape = b^2)
  return(Y1 / (Y1 + Y2))
}

#' @rdname bessel
#' @export
pbessel <- function(q, mu = 1/2, phi = 1){
  prob_bessel <- sapply(q, function(v){
    if(v<=0){
      return(0)
    } else if(v>=1){
      return(1)
    } else{
      stats::integrate(dbessel,lower = 0, upper = v, mu=mu, phi=phi)$value
    }
  })
  return(prob_bessel)
}

#' @rdname bessel
#' @export
qbessel <- function(p, mu = 1/2, phi = 1){
  quant_bessel <- sapply(p, function(x){
    if(x < 0 | x > 1){
      warn_qbessel <- TRUE
      return(NaN)} else{
        return(stats::uniroot(function(y) {
          pbessel(y, mu=mu,phi=phi) - x
        },lower = 0, upper = 1)$root)
      }
  }
  )
  if(any(is.nan(quant_bessel))){
    warning("NaNs produced", call. = TRUE, domain = "R")
  }
  return(quant_bessel)
}

#############################################################################################
#' @title quantile_residual_bes
#' @description Function to calculate quantile residuals based on the bessel regression. Details about this type of residual can be found in \emph{Pereira (2019)}.
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @references DOI:10.1080/03610918.2017.1381740 (\href{https://www.tandfonline.com/doi/abs/10.1080/03610918.2017.1381740}{Pereira; 2019})
#' @seealso
#' \code{\link{score_residual_bes}}
#' @return Vector containing the quantile residuals.
quantile_residual_bes <- function(kap, lam, z, x, v, link.mean, link.precision) {
  n <- length(z)
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  mu <- link_mean$linkinv(x %*% kap)
  phi <- link_precision$linkinv(v %*% lam)
  out <- rep(0, n)
  for (i in 1:n) {
    aux <- stats::integrate(f = dbessel, lower = 0, upper = z[i], mu[i], phi[i])
    aux <- aux$value
    if (aux > 0.99999) {
      aux <- 0.99999
    }
    if (aux < 0.00001) {
      aux <- 0.00001
    }
    out[i] <- stats::qnorm(aux, mean = 0, sd = 1)
  }
  return(out)
}
