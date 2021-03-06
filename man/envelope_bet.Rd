% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_em_beta.R
\name{envelope_bet}
\alias{envelope_bet}
\title{envelope_bet}
\usage{
envelope_bet(
  residual,
  kap,
  lam,
  x,
  v,
  nsim_env,
  prob,
  n,
  link.mean,
  link.precision,
  em_controls,
  optim_method,
  optim_controls
)
}
\arguments{
\item{residual}{character indicating the type of residual ("pearson", "score" or "quantile").}

\item{kap}{coefficients in kappa related to the mean parameter.}

\item{lam}{coefficients in lambda related to the precision parameter.}

\item{x}{matrix containing the covariates for the mean submodel. Each column is a different covariate.}

\item{v}{matrix containing the covariates for the precision submodel. Each column is a different covariate.}

\item{nsim_env}{number of synthetic data sets to be generated.}

\item{prob}{confidence level of the envelope (number between 0 and 1).}

\item{n}{sample size.}

\item{link.mean}{a string containing the link function for the mean.
The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".}

\item{link.precision}{a string containing the link function the precision parameter.
The possible link functions for the precision parameter are "identity", "log", "sqrt".}

\item{em_controls}{a list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm;
\code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm.}

\item{optim_method}{main optimization algorithm to be used. The available methods are the same as those of \code{optim} function.}

\item{optim_controls}{a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.}
}
\value{
Matrix with dimension 2 x n (1st row = upper bound, second row = lower bound).
}
\description{
Function to calculate envelopes based on residuals for the beta regression.
}
\seealso{
\code{\link{score_residual_bet}}, \code{\link{quantile_residual_bet}}, \code{\link{pred_accuracy_bet}}
}
