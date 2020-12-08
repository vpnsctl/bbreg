
#############################################################################################
#' @title startvalues
#' @description Function providing initial values for the Expectation-Maximization algorithm.
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @param model The model to which the initial values will be applied, "bessel" or "beta".
startvalues <- function(z, x, v, link.mean, link.precision, model) {
  nkap <- ncol(x)
  nlam <- ncol(v)
  n <- length(z)
  link_mean <- stats::make.link(link.mean)
  link_precision <- stats::make.link(link.precision)
  
  fit_aux <- stats::glm.fit(x = x, y = z, family = stats::quasibinomial(link = link.mean))
  kap_start <- fit_aux$coefficients
  
  mu_est <- link_mean$linkinv(x %*% kap_start)
  
  g_phi <- sum((z - mu_est)^2 / (mu_est * (1 - mu_est))) / (n - nkap)
  
  phi_est <- g_inv(g_phi, model)
  
  lam_start <- c(link_precision$linkfun(phi_est), rep(0, nlam - 1))
  
  out <- list()
  
  out[[1]] <- kap_start
  out[[2]] <- lam_start
  return(out)
}

#############################################################################################
#' @title d2mudeta2
#' @description Function to obtain the second derivatives of the mean parameter with respect to the linear predictor.
#' @param link.mean a string containing the link function for the mean.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param mu mean parameter.

d2mudeta2 <- function(link.mean, mu) {
  d2mu <- switch(link.mean,
                 "logit" = {
                   mu * (1 - mu) * (1 - 2 * mu)
                 },
                 "probit" = {
                   (-mu / sqrt(2 * pi)) * exp(-mu^2 / 2)
                 },
                 "cloglog" = {
                   -(1 - mu) * log(1 - mu) * (1 + log(1 - mu))
                 },
                 "cauchit" = {
                   (-2 / pi) * tan(pi * (mu - 1 / 2)) / ((1 + tan(pi * (mu - 1 / 2))^2)^2)
                 }
  )
  return(d2mu)
}

#############################################################################################
#' @title d2phideta2
#' @description Function to obtain the second derivatives of the precision parameter with respect to the linear predictor.
#' @param link.precision a string containing the link function the precision parameter.
#' The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @param phi precision parameter.

d2phideta2 <- function(link.precision, phi) {
  d2phi <- switch(link.precision,
                  "identity" = {
                    0
                  },
                  "log" = {
                    phi
                  },
                  "sqrt" = {
                    2
                  },
                  "inverse" = {
                    2 * phi^3
                  },
                  "1/precision^2" = {
                    3 * phi^5 / 4
                  }
  )
  return(d2phi)
}

#############################################################################################
#' @title g_phi
#' @description Function to obtain the function g() that relates the precision parameter to the dispersion parameter for bessel and beta distributions.
#' @param phi precision parameter.
#' @param model "bessel" or "beta"

g_phi <- function(phi, model) {
  if (model == "beta") {
    return(1 / (1 + phi))
  } else if (model == "bessel") {
    return((1 - phi + (phi^2) * exp(phi) * expint_En(phi, order = 1)) / 2)
  }
}

#############################################################################################
#' @title g_inv
#' @description Function to obtain an approximation of the inverse of the function g() that relates the precision parameter to the dispersion parameter for bessel and beta distributions.
#' @param sigma2 dispersion parameter.
#' @param model "bessel" or "beta"

g_inv <- function(sigma2, model) {
  if (model == "beta") {
    return(1 / sigma2 - 1)
  } else if (model == "bessel") {
    if (sigma2 >= 0.3) {
      return(0.5)
    } else if (sigma2 >= 0.15) {
      return(2)
    } else if (sigma2 >= 0.07) {
      return(8)
    } else if (sigma2 >= 0.044) {
      return(15)
    } else if (sigma2 >= 0.023) {
      return(30)
    } else if (sigma2 >= 0.01) {
      return(80)
    } else if (sigma2 >= 0.005) {
      return(150)
    } else if (sigma2 >= 0.0025) {
      return(300)
    } else if (sigma2 >= 0.0015) {
      return(525)
    } else {
      return(650)
    }
  }
}
