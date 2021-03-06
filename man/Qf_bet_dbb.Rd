% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5_em_dbb.R
\name{Qf_bet_dbb}
\alias{Qf_bet_dbb}
\title{Qf_bet_dbb}
\usage{
Qf_bet_dbb(lam, phiold, z, v, mu, link.precision)
}
\arguments{
\item{lam}{coefficients in lambda related to the covariates in v.}

\item{phiold}{previous value of the precision parameter (phi).}

\item{z}{response vector with 0 < z_i < 1.}

\item{v}{matrix containing the covariates for the precision submodel. Each column is a different covariate.}

\item{mu}{mean parameter (vector having the same size of z).}

\item{link.precision}{a string containing the link function the precision parameter.
The possible link functions for the precision parameter are "identity", "log", "sqrt".}
}
\value{
Scalar representing the output of this auxiliary function for the beta case.
}
\description{
Q-function related to the beta model. This function was adapted for the discrimination test between bessel and beta (DBB) required in the Expectation-Maximization algorithm.
}
