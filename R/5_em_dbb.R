#############################################################################################
#' @title dbbtest
#' @description Function to run the discrimination test between beta and bessel regressions (DBB).
#' @param formula symbolic description of the model (set: z ~ x or z ~ x | v); see details below.
#' @param data arguments considered in the formula description. This is usually a data frame composed by:
#' (i) the response with bounded continuous observations (0 < z_i < 1),
#' (ii) covariates for the mean submodel (columns of matrix x) and
#' (iii) covariates for the precision submodel (columns of matrix v).
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @param em_controls a list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm, the default is set to 5000; 
#' \code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm, the default is set to 10^(-5).
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function. The default is set to "L-BFGS-B".
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @return Object of class dbbtest, which is a list containing two elements. The 1st one is a table of terms
#' considered in the decision rule of the test; they are sum(z2/n) = sum_{i=1}^{n}(z_i^2)/n, sum(quasi_mu) = sum_{i=1}^{n}(tilde{mu_i}^2 + tilde{mu_i}(1-tilde{mu_i})/2)
#' |D_bessel| and |D_beta| as indicated in the main reference. The 2nd term of the list is the name of the selected model (bessel or beta).
#' @seealso
#' \code{\link{simdata_bes}}, \code{\link{dbessel}}, \code{\link{simdata_bet}}
#' @examples
#' # Illustration using the Weather task data set available in the bbreg package.
#' dbbtest(agreement ~ priming + eliciting,
#'   data = WT,
#'   link.mean = "logit", link.precision = "identity"
#' )
#' @export
dbbtest <- function(formula, data, link.mean, link.precision, em_controls = list(maxit = 5000, em_tol = 10^(-5)), optim_method = "L-BFGS-B", optim_controls = list()) {
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  ## Processing call
  eps <- match.call()
  if (missing(data)) {
    data <- environment(formula)
  }
  MF <- match.call(expand.dots = FALSE)
  arg_names <- match(c("formula", "data"), names(MF), nomatch = 0)
  MF <- MF[c(1, arg_names)]
  MF$drop.unused.levels <- TRUE
  ## Processing formula
  Fo <- as.Formula(formula)
  if (length(Fo)[2] < 2) {
    Fo <- as.Formula(stats::formula(Fo), ~1)
  } else {
    if (length(Fo)[2] > 2) {
      Fo <- Formula(stats::formula(Fo, rhs = 1:2))
    }
  }
  MF$formula <- Fo
  MF[[1]] <- as.name("model.frame")
  MF <- eval(MF, parent.frame())
  MTerms_x <- stats::terms(Fo, data = data, rhs = 1)
  MTerms_v <- stats::delete.response(stats::terms(Fo, data = data, rhs = 2))
  z <- stats::model.response(MF, "numeric")
  x <- stats::model.matrix(MTerms_x, MF)
  v <- stats::model.matrix(MTerms_v, MF)
  
  out <- dbbtest.fit(z = z,x = x,v = v, 
                     link.mean = link.mean, 
                     link.precision = link.precision, em_controls = em_controls, optim_method = optim_method, optim_controls = optim_controls)
  return(out)
}


#############################################################################################
#' @title dbbtest.fit
#' @description Function to run the discrimination test between beta and bessel regressions (DBB).
#' @param z vector of response variables with length \code{n}. Each coordinate must belong to the standard unit interval (0,1). 
#' @param x matrix of covariates with respect to the mean with dimension \code{(n,nkap)}.
#' @param v matrix of covariates with respect to the precision parameter. The default is \code{NULL}. If not \code{NULL} must be of dimension \code{(n,nlam)}.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @param em_controls a list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm, the default is set to 5000; 
#' \code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm, the default is set to 10^(-5).
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function. The default is set to "L-BFGS-B".
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @return Object of class dbbtest, which is a list containing two elements. The 1st one is a table of terms
#' considered in the decision rule of the test; they are sum(z2/n) = sum_{i=1}^{n}(z_i^2)/n, sum(quasi_mu) = sum_{i=1}^{n}(tilde{mu_i}^2 + tilde{mu_i}(1-tilde{mu_i})/2)
#' |D_bessel| and |D_beta| as indicated in the main reference. The 2nd term of the list is the name of the selected model (bessel or beta).
#' @seealso
#' \code{\link{simdata_bes}}, \code{\link{dbbtest}}, \code{\link{simdata_bet}}
#' @examples
#' # Illustration using the Weather task data set available in the bbreg package.
#' n <- 100
#' x <- cbind(rbinom(n, 1, 0.5), runif(n, -1, 1))
#' v <- runif(n, -1, 1)
#' z <- simdata_bes(
#'   kap = c(1, -1, 0.5), lam = c(0.5, -0.5), x, v,
#'   repetition = 1, link.mean = "logit", link.precision = "log"
#' )
#' z <- unlist(z)
#' dbbtest.fit(z = z, x = x ,v = v, link.mean = "logit", link.precision = "identity")
#' @export
dbbtest.fit <- function(z,x,v = NULL, link.mean, link.precision, em_controls = list(maxit = 5000, em_tol = 10^(-5)), optim_method = "L-BFGS-B", optim_controls = list()) {
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  
  x <- as.matrix(x)
  
  n <- length(z)
  if(is.null(v)){
    v = matrix(rep(1,n), nrow = n)
  } else{
    v = as.matrix(v)
  }
  nkap <- ncol(x)
  nlam <- ncol(v)
  if (nkap == 0) {
    x <- cbind(rep(1, n), x)
    nkap <- 1
  }
  if (nlam == 0) {
    v <- cbind(rep(1, n), v)
    nlam <- 1
  }
  #
  idx <- which(apply(x == 1, 2, sum) != n) # find non-intercept columns in x.
  idv <- which(apply(v == 1, 2, sum) != n) # find non-intercept columns in v.
  #
  if (length(idx) > 0 & length(idx) == nkap) {
    outquasi <- stats::glm(z ~ 0 + x[, idx], family = stats::quasi(variance = "mu(1-mu)", link = link.mean), start = rep(0, nkap))
  }
  if (length(idx) > 0 & length(idx) < nkap) {
    outquasi <- stats::glm(z ~ x[, idx], family = stats::quasi(variance = "mu(1-mu)", link = link.mean), start = rep(0, nkap))
  }
  if (length(idx) == 0) {
    outquasi <- stats::glm(z ~ 1, family = stats::quasi(variance = "mu(1-mu)", link = link.mean), start = rep(0, nkap))
  }
  kapquasi <- outquasi$coefficients
  muquasi <- link_mean$linkinv(x %*% kapquasi)
  sumz2 <- sum(z^2) / n
  sumquasi <- sum(muquasi * (1 - muquasi) / 2 + muquasi^2)
  selection <- 0 # symbol: 0 = beta and 1 = bessel
  if (sumz2 < sumquasi) {
    # Set initial values for bessel.
    lam <- startvalues(z, x, v, link.mean, link.precision, "bessel")
    lam <- lam[[2]]
    #
    EM <- EMrun_bes_dbb(lam, z, v, mu = muquasi, link.precision, em_controls, optim_method, optim_controls)
    phi <- link_precision$linkinv(v %*% EM)
    gphi <- (1 - phi + (phi^2) * exp(phi) * expint_En(phi, order = 1)) / 2
    Wbes <- gphi
    # Set initial values for beta
    lam <- startvalues(z, x, v, link.mean, link.precision, "beta")
    lam <- lam[[2]]
    #
    EM <- EMrun_bet_dbb(lam, z, v, mu = muquasi, link.precision, em_controls, optim_method, optim_controls)
    phi <- link_precision$linkinv(v %*% EM)
    gphi <- 1 / (1 + phi)
    Wbet <- gphi
    #
    Dbes <- sumz2 - sum(muquasi * (1 - muquasi) * Wbes + muquasi^2) / n
    Dbet <- sumz2 - sum(muquasi * (1 - muquasi) * Wbet + muquasi^2) / n
    if (abs(Dbes) <= abs(Dbet)) {
      selection <- 1
    }
  }
  out <- list()
  tab <- c(sumz2, sumquasi, abs(Dbes), abs(Dbet))
  names(tab) <- c("sum(z2/n)", "sum(quasi_mu)", "|D_bessel|", "|D_beta|")
  out[[1]] <- tab
  out[[2]] <- if (selection == 0) {
    "beta"
  } else {
    "bessel"
  }
  names(out) <- c("terms", "model")
  class(out) <- "dbbtest"
  return(out)
}





#############################################################################################
#' @title Qf_bes_dbb
#' @description Q-function related to the bessel model. This function was adapted for the discrimination test between bessel and beta (DBB) required in the Expectation-Maximization algorithm.
#' @param lam coefficients in lambda related to the covariates in v.
#' @param wz parameter wz representing E(1/W_i|Z_i = z_i, theta).
#' @param z response vector with 0 < z_i < 1.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param mu mean parameter (vector having the same size of z).
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return Scalar representing the output of this auxiliary function for the bessel case.
Qf_bes_dbb <- function(lam, wz, z, v, mu, link.precision) {
  link_precision <- stats::make.link(link.precision)
  phi <- link_precision$linkinv(v %*% lam) # precision parameter
  out1 <- log(mu) + log(1 - mu) + 2 * log(phi) + phi
  out2 <- 0.5 * wz * (phi^2) * (((mu^2) / z) + (((1 - mu)^2) / (1 - z)))
  return(sum(out1 - out2))
}

#############################################################################################
#' @title Qf_bet_dbb
#' @description Q-function related to the beta model. This function was adapted for the discrimination test between bessel and beta (DBB) required in the Expectation-Maximization algorithm.
#' @param lam coefficients in lambda related to the covariates in v.
#' @param phiold previous value of the precision parameter (phi).
#' @param z response vector with 0 < z_i < 1.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param mu mean parameter (vector having the same size of z).
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return Scalar representing the output of this auxiliary function for the beta case.
Qf_bet_dbb <- function(lam, phiold, z, v, mu, link.precision) {
  link_precision <- stats::make.link(link.precision)
  phi <- link_precision$linkinv(v %*% lam) # precision parameter
  mu0phi <- mu * phi
  mu1phi <- (1 - mu) * phi
  mu0phi[which(mu0phi <= 0)] <- 10^(-10)
  mu1phi[which(mu1phi <= 0)] <- 10^(-10)
  out <- phi * (mu * log(z / (1 - z)) + digamma(phiold) + log(1 - z))
  out <- out - lgamma(mu0phi) - lgamma(mu1phi)
  out <- out - log(z * (1 - z)) - digamma(phiold) - log(1 - z) - phiold
  return(sum(out))
}

#############################################################################################
#' @title gradlam_bes_dbb
#' @description Gradient of the Q-function (adapted for the discrimination test between bessel and beta - DBB) to calculate the gradient required for optimization via \code{optim}.
#' This option is related to the bessel regression.
#' @param lam coefficients in lambda related to the covariates in v.
#' @param wz parameter wz representing E(1/W_i|Z_i = z_i, theta).
#' @param z response vector with 0 < z_i < 1.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param mu mean parameter (vector having the same size of z).
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return Scalar representing the output of this auxiliary gradient function for the bessel case.
gradlam_bes_dbb <- function(lam, wz, z, v, mu, link.precision) {
  link_precision <- stats::make.link(link.precision)
  dphideta <- link_precision$mu.eta(v %*% lam)
  phi <- link_precision$linkinv(v %*% lam)
  aux <- (2 / phi) + 1 - wz * phi * (1 + ((z - mu)^2) / (z * (1 - z)))
  aux <- aux * dphideta
  Ulam <- t(v) %*% aux
  return(Ulam)
}

#############################################################################################
#' @title gradlam_bet
#' @description Gradient of the Q-function (adapted for the discrimination test between bessel and beta - DBB) to calculate the gradient required for optimization via \code{optim}.
#' This option is related to the beta regression.
#' @param lam coefficients in lambda related to the covariates in v.
#' @param phiold previous value of the precision parameter (phi).
#' @param z response vector with 0 < z_i < 1.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param mu mean parameter (vector having the same size of z).
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @return Scalar representing the output of this auxiliary gradient function for the beta case.
gradlam_bet_dbb <- function(lam, phiold, z, v, mu, link.precision) {
  link_precision <- stats::make.link(link.precision)
  phi <- link_precision$linkinv(v %*% lam)
  dphideta <- link_precision$mu.eta(v %*% lam)
  aux <- mu * (log(z) - log(1 - z)) + digamma(phiold)
  aux <- aux + log(1 - z) - mu * digamma(mu * phi)
  aux <- aux - (1 - mu) * digamma((1 - mu) * phi)
  aux <- aux * dphideta
  Ulam <- t(v) %*% aux
  return(Ulam)
}

#############################################################################################
#' @title EMrun_bes_dbb
#' @description Function (adapted for the discrimination test between bessel and beta - DBB) to run the Expectation-Maximization algorithm for the bessel regression.
#' @param lam initial values for the coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param mu mean parameter (vector having the same size of z).
#' @param link.precision a string containing the link function the precision parameter.
#' @param em_controls a list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm, the default is set to 5000; 
#' \code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm, the default is set to 10^(-5).
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function. The default is set to "L-BFGS-B".
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @return Vector containing the estimates for lam in the bessel regression.
EMrun_bes_dbb <- function(lam, z, v, mu, link.precision, em_controls = list(maxit = 5000, em_tol = 10^(-5)), optim_method = "L-BFGS-B", optim_controls = list()) {
  link_precision <- stats::make.link(link.precision)
  phi <- link_precision$linkinv(v %*% lam) # phi precision parameter.
  count <- 0
  
  maxit <- em_controls$maxit
  epsilon <- em_controls$em_tol
  
  repeat{
    lam_r <- lam
    phi_r <- phi
    ### E step ------------------------------
    wz_r <- Ew1z(z, mu, phi_r)
    ### M step ------------------------------
    M <- tryCatch(stats::optim(
      par = lam, fn = Qf_bes_dbb, gr = gradlam_bes_dbb, wz = wz_r, z = z, v = v,
      mu = mu, link.precision = link.precision, control = c(fnscale = -1, optim_controls), method = optim_method
    ), error = function(e) {
      "Error"
    })
    
    if (length(M) == 1) {
      warning("Trying with numerical derivatives")
      M <- tryCatch(stats::optim(
        par = lam, fn = Qf_bes_dbb, gr = NULL, wz = wz_r, z = z, v = v,
        mu = mu, link.precision = link.precision, control = c(fnscale = -1, optim_controls), method = optim_method
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
        par = lam, fn = Qf_bes_dbb, gr = NULL, wz = wz_r, z = z, v = v,
        mu = mu, link.precision = link.precision, control = c(fnscale = -1, optim_controls), method = optim_temp
      ), error = function(e) {
        "Error"
      })
    }
    if (length(M) == 1) {
      warning("The EM algorithm did not converge.")
      break
    }
    lam <- M$par
    phi <- link_precision$linkinv(v %*% lam)
    # Compute Q -----------------------------
    Q_r <- Qf_bes_dbb(lam_r, wz_r, z, v, mu, link.precision)
    Q <- Qf_bes_dbb(lam, wz_r, z, v, mu, link.precision)
    ### Convergence criterion ---------------
    term1 <- sqrt(sum((lam - lam_r)^2))
    term2 <- abs(Q - Q_r)
    ### -------------------------------------
    count <- count + 1
    if (max(term1, term2) < epsilon) {
      break
    }
    if (count >= maxit) {
      warning("The EM algorithm did not converge.")
      break
    }
    ### -------------------------------------
  }
  return(lam)
}

#############################################################################################
#' @title EMrun_bet_dbb
#' @description Function (adapted for the discrimination test between bessel and beta - DBB) to run the Expectation-Maximization algorithm for the beta regression.
#' @param lam initial values for the coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param mu mean parameter (vector having the same size of z).
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @param em_controls a list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm, the default is set to 5000; 
#' \code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm, the default is set to 10^(-5).
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function. The default is set to "L-BFGS-B".
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @return Vector containing the estimates for lam in the beta regression.
EMrun_bet_dbb <- function(lam, z, v, mu, link.precision, em_controls = list(maxit = 5000, em_tol = 10^(-5)), optim_method = "L-BFGS-B", optim_controls = list()) {
  link_precision <- stats::make.link(link.precision)
  phi <- link_precision$linkinv(v %*% lam) # phi precision parameter.
  count <- 0
  
  maxit <- em_controls$maxit
  epsilon <- em_controls$em_tol

  repeat{
    lam_r <- lam
    phi_r <- phi
    ### E step ------------------------------
    ### M step ------------------------------
    M <- tryCatch(stats::optim(
      par = lam, fn = Qf_bet_dbb, gr = gradlam_bet_dbb, phiold = phi_r, z = z, v = v, mu = mu,
      link.precision = link.precision, control = c(fnscale = -1, optim_controls), method = optim_method
    ), error = function(e) {
      "Error"
    })
    
    if (length(M) == 1) {
      warning("Trying with numerical derivatives")
      M <- tryCatch(stats::optim(
        par = lam, fn = Qf_bet_dbb, gr = NULL, phiold = phi_r, z = z, v = v, mu = mu,
        link.precision = link.precision, control = c(fnscale = -1, optim_controls), method = optim_method
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
        par = lam, fn = Qf_bet_dbb, gr = NULL, phiold = phi_r, z = z, v = v, mu = mu,
        link.precision = link.precision, control = c(fnscale = -1, optim_controls), method = optim_temp
      ), error = function(e) {
        "Error"
      })
    }
    if (length(M) == 1) {
      warning("The EM algorithm did not converge.")
      break
    }
    lam <- M$par
    phi <- link_precision$linkinv(v %*% lam)
    # Compute Q -----------------------------
    Q_r <- Qf_bet_dbb(lam_r, phi_r, z, v, mu, link.precision)
    Q <- Qf_bet_dbb(lam, phi_r, z, v, mu, link.precision)
    ### Convergence criterion ---------------
    term1 <- sqrt(sum((lam - lam_r)^2))
    term2 <- abs(Q - Q_r)
    ### -------------------------------------
    count <- count + 1
    if (max(term1, term2) < epsilon) {
      break
    }
    if (count >= maxit) {
      warning("The EM algorithm did not converge.")
      break
    }
    ### -------------------------------------
  }
  return(lam)
}
