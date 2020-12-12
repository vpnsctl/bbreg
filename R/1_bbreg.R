#############################################################################################
#' @name bbreg
#' @title Bessel and Beta Regression Models
#' @aliases bbreg bbreg.fit
#' @description Function to fit, via Expectation-Maximization (EM) algorithm, the bessel or the beta regression to a given data set with a bounded continuous response variable.
#' @param formula symbolic description of the model (examples: \code{z ~ x1 + x2} and \code{z ~ x1 + x2 | v1 + v2}); see details below.
#' @param data elements expressed in formula. This is usually a data frame composed by:
#' (i) the bounded continuous observations in \code{z} (0 < z_i < 1),
#' (ii) covariates for the mean submodel (columns \code{x1} and \code{x2}) and
#' (iii) covariates for the precision submodel (columns \code{v1} and \code{v2}).
#' @param model character ("bessel" or "beta") indicating the type of model to be fitted. The default is NULL,
#' meaning that a discrimination test must be applied to select the model.
#' @param z vector of response variables with length \code{n}. Each coordinate must belong to the standard unit interval (0,1). 
#' @param x matrix of covariates with respect to the mean with dimension \code{(n,nkap)}.
#' @param v matrix of covariates with respect to the precision parameter. The default is \code{NULL}. If not \code{NULL} must be of dimension \code{(n,nlam)}.
#' @param residual character indicating the type of residual to be evaluated ("pearson", "score" or "quantile"). The default is "pearson".
#' @param envelope number of simulations (synthetic data sets) to build envelopes for residuals (with \code{100*prob\%} confidence level).
#' The default \code{envelope = 0} dismisses the envelope analysis.
#' @param prob probability indicating the confidence level for the envelopes (default: \code{prob} = 0.95).
#' If \code{envelope} = 0, \code{prob} is ignored.
#' @param predict number of partitions (training set to fit the model and a test set to calculate residuals) to be evaluated in a predictive accuracy
#' study related to the \code{RSS_pred} (residual sum of squares for the partition "test set"). The default \code{predict} = 0 dismisses the \code{RSS_pred} analysis.
#' The partitions are randomly defined when predict is set as a positive integer.
#' @param ptest proportion of the sample size to be considered in the test set for the \code{RSS_pred} analysis (default = 0.25 = 25\% of the sample size).
#' If predict = 0, ptest is ignored.
#' @param link.mean optionally, a string containing the link function for the mean. If omitted, the 'logit' link function will be used.
#' The possible link functions for the mean are "logit","probit", "cauchit", "cloglog".
#' @param link.precision optionally, a string containing the link function the precision parameter. If omitted and the only precision
#' covariate is the intercept, the 'identity' link function will be used, if omitted and there is a precision covariate other than the
#' intercept, the 'log' link function will be used. The possible link functions for the precision parameter are "identity", "log", "sqrt".
#' @param em_controls a list containing two elements: \code{maxit} that contains the maximum number of iterations of the EM algorithm, the default is set to 5000; 
#' \code{em_tol} that defines the tolerance value to control the convergence criterion in the EM-algorithm, the default is set to 10^(-5).
#' @param optim_method main optimization algorithm to be used. The available methods are the same as those of \code{optim} function. The default is set to "L-BFGS-B".
#' @param optim_controls a list of control arguments to be passed to the \code{optim} function in the optimization of the model. For the control options, see
#' the 'Details' in the help of \code{\link[stats]{optim}} for the possible arguments.
#' @return \code{bbreg} returns an object of class "bbreg". The function \code{bbreg.fit} returns an object of class "bbreg_fit".
#' The objects of classes \code{bbreg} and \code{bbreg_fit} return lists containing:
#' 
#'  \itemize{
#'   \item \code{coefficients} - a list with elements "mean" and "precision" containing the estimated coefficients of the model;
#'   \item \code{call} - the formula used by the model. If using \code{bbreg.fit}, this returns \code{NULL}.
#'   \item \code{modelname} - the fitted model, bessel or beta;
#'   \item \code{message} - message to be displayed, regarding the usage or not of the discrimination criterion;
#'   \item \code{residualname} - the name of the chosen residual in the call;
#'   \item \code{niter} - number of iterations of the EM algorithm;
#'   \item \code{start} - the initial guesses of the parameters
#'   \item \code{intercept} - vector indicating if the intercept is present in the mean and/or in the precision regressions;
#'   \item \code{link.mean} - link function of the mean;
#'   \item \code{link.precision} - link function of the precision parameter;
#'   \item \code{kappa} - vector containing the estimates of the mean-related coefficients;
#'   \item \code{lambda} - vector containing the estimates of the precision-related coefficients;
#'   \item \code{mu} - estimated means;
#'   \item \code{fitted.values} - the fitted values in the response scale;
#'   \item \code{efron.pseudo.r2} - Efron's pseudo R^2: the squared correlation between the response variables and the predicted values;
#'   \item \code{x} - the covariates related to the mean;
#'   \item \code{v} - the covariates related to the precision parameter;
#'   \item \code{z} - the response variables;
#'   \item \code{gphi} - the estimated dispersion parameters;
#'   \item \code{residuals} - the values of the chosen residual in the call;
#'   \item \code{std_errors} - the standard errors of the estimated parameters;
#'   \item \code{RSS} - sum of squared residuals;
#'   \item \code{RSS_pred} - sum of squared residuals from the predictions (when prediction is performed);
#'   \item \code{DBB} - vector containing the discrimination statistics (when discrimination is performed);
#'   \item \code{envelope} - the numerical envelopes used to build the Q-Q plot with simulated envelopes;
#'   \item \code{terms} - (only for \code{bbreg})the \code{terms} object used;
#'   \item \code{levels} - (where relevant, only for \code{bbreg}) the levels of the factors used;
#'   \item \code{contrasts} - (where relevant, only for \code{bbreg}) the contrasts used.
#' }
#'
#'
#' @details The bessel regression originates from a class of normalized inverse-Gaussian (N-IG) process introduced in \emph{Lijoi et al. (2005)}
#' as an alternative to the widely used Dirichlet process in the Bayesian context. These authors consider a ratio of inverse-Gaussian random variables
#' to define the new process. In the particular univariate case, the N-IG is obtained from the representation "Z = Y1/(Y1+Y2)", with "Y1" and "Y2" being
#' independent inverse-Gaussian random variables having scale = 1 and shape parameters "a1 > 0" and "a2 > 0", respectively.
#' Denote "Y1 ~ IG(a1)" and "Y2 ~ IG(a2)". The density of "Z" has support in the interval (0,1) and it depends on the modified Bessel function of third
#' kind with order 1, named here as "\emph{K1(-)}". The presence of "\emph{K1(-)}" in the structure of the p.d.f. establishes the name of the new distribution;
#' consider Z ~ Bessel(a1,a2). Note that the name "beta distribution" is also an analogy to the presence of a function (the beta function)
#' in its p.d.f. structure. The bessel regression model is defined by assuming "Z_1,...,Z_n" as a random sample of continuous bounded responses with
#' "Z_i ~ Bessel(mu_i,phi_i)" for "i = 1,...,n". Using this parameterization, one can write: "E(Z_i) = mu_i" and "Var(Z_i) = mu_i(1-mu_i) g(phi_i)",
#' where "\emph{g(-)}" is a function depending on the exponential integral of "phi_i". The following link functions are assumed "logit(mu_i) = x_i^T kappa" and
#' "log(phi_i) = v_i^T lambda", where "kappa' = (kappa_1,...,kappa_p)" and "lambda' = (lambda_1,...,lambda_q)" are real valued vectors.
#' The terms "x_i^T" and "v_i^T" represent, respectively, the i-th row of the matrices "x" (\emph{nxp}) and "v" (\emph{nxq}) containing covariates in their columns
#' ("x_{i,1}" and "v_{i,1}" may be 1 to handle intercepts). As it can be seen, this regression model has two levels with covariates explaining the mean
#' "mu_i" and the parameter "phi_i". For more details about the bessel regression see \emph{Barreto-Souza, Mayrink and Simas (2020)}.
#'
#' This package implements an Expectation Maximization (EM) algorithm to fit the bessel regression. The full EM approach proposed in \emph{Barreto-Souza and Simas (2017)} for the beta
#' regression is also available here. Fitting the beta regression via EM-algorithm is a major difference between the present package \pkg{bbreg} and the
#' well known \code{betareg} created by Alexandre B. Simas and currently maintained by Achim Zeileis. The estimation procedure on the \code{betareg} packages
#' is given by maximizing the beta model likelihood via \code{\link[stats]{optim}}.
#' In terms of initial values, \pkg{bbreg} uses quasi-likelihood estimates as the starting points for
#' the EM-algorithms. The formulation of the target model also has the same structure as in the standard functions \code{lm}, \code{glm} and \code{betareg},
#' with also the same structure as the latter when precision covariates are being used. The user is supposed to
#' write a formula object describing elements of the regression (response, covariates for the mean submodel,
#' covariates for the precision submodel, presence of intercepts, and interactions). As an example, the description
#' "z ~ x" indicates: "response = z" (continuous and bounded by 0 and 1), "covariates = columns of x" (mean submodel) and
#' precision submodel having only an intercept. On the other hand, the configuration "z ~ x | v" establishes that the covariates given
#' in the columns of "v" must be used in the precision submodel. Intercepts may be removed by setting
#' "z ~ 0 + x | 0 + v" or "z ~ x - 1|v - 1". Absence of intercept and covariates is not allowed in any submodel.
#' The type of model to be fitted ("bessel" or "beta") can be specified through the argument "model" of
#' \pkg{bbreg}. If the user does not specify the model, the package will automatically apply a discrimination
#' test (DBB - Discrimination between Bessel and Beta),
#' developed in \emph{Barreto-Souza, Mayrink and Simas (2020)}, to select the most appropriate model for the given
#' data set. In this case, some quantities related to the DBB are included in the final output; they are:
#' "sum(Z2/n)" = mean of z_i^2, "sum(quasi_mu)" = sum (for i = 1,...,n) of muq_i + muq_i(1-muq_i)/2,
#' with muq_i being the quasi-likelihood estimator of mu_i and, finally, the quantities "|D_bessel|" and
#' "|D_beta|" depending on muq_i and the EM-estimates of phi_i under bessel or beta.
#'
#' In the current version, three types of residuals are available for analysis ("Pearson", "Score" and "Quantile").
#' The user may choose one of them via the argument "residual". The score residual is computed empirically, based
#' on 100 artificial data sets generated from the fitted model. The sample size
#' is the same of the original data and the simulations are used to estimate the mean and s.d. required in the score
#' residual formulation. The user
#' may also choose to build envelopes for the residuals with confidence level in "prob". This feature also requires simulations of synthetic data
#' ("envelope" is the number of replications). Residuals are obtained for each data set and confronted against the quantiles of the N(0,1). Predictive
#' accuracy of the fitted model is also explored by setting "predict" as a positive integer (this value represents the number of random partitions to be evaluated).
#' In this case, the full data set is separated in a training (partition to fit the model) and a test set (to evaluate residuals) for which the
#' RSS (Residual Sum of Squares) is computed. The default partition is 75\% (training) and 25\% (test); this can be modified by choosing the
#' proportion \code{ptest} for the test set (large \code{ptest} is not recommended).
#'
#' @references
#' arXiv:2003.05157 (\href{https://arxiv.org/abs/2003.05157}{Barreto-Souza, Mayrink and Simas; 2020})
#'
#' DOI:10.1080/00949655.2017.1350679 (\href{https://www.tandfonline.com/doi/abs/10.1080/00949655.2017.1350679?journalCode=gscs20}{Barreto-Souza and Simas; 2017})
#'
#' DOI:10.18637/jss.v034.i02 (\href{https://www.jstatsoft.org/article/view/v034i02}{Cribari-Neto and Zeileis; 2010})
#'
#' DOI:10.1198/016214505000000132 (\href{https://www.tandfonline.com/doi/abs/10.1198/016214505000000132}{Lijoi et al.; 2005})
#'
#' @seealso
#' \code{\link{summary.bbreg}}, \code{\link{plot.bbreg}}, \code{\link{simdata_bes}}, \code{\link{dbessel}}, \code{\link{dbbtest}}, \code{\link{simdata_bet}}, \code{\link[Formula]{Formula}}
#'
#' @examples
#' # Example with artificial data.
#' n <- 100
#' x <- cbind(rbinom(n, 1, 0.5), runif(n, -1, 1))
#' v <- runif(n, -1, 1)
#' z <- simdata_bes(
#'   kap = c(1, -1, 0.5), lam = c(0.5, -0.5), x, v,
#'   repetition = 1, link.mean = "logit", link.precision = "log"
#' )
#' z <- unlist(z)
#' fit1 <- bbreg(z ~ x | v)
#' summary(fit1)
#' plot(fit1)
#'
#' # Examples using the Weather Task (WT) data available in bbreg.
#' \donttest{
#' fit2 <- bbreg(agreement ~ priming + eliciting, data = WT)
#' summary(fit2)
#' }
#' \donttest{
#' fit3 <- bbreg(agreement ~ priming + eliciting, envelope = 30, predict = 10, data = WT)
#' summary(fit3)
#' }
#' # Example with precision covariates
#' \donttest{
#' fit4 <- bbreg(agreement ~ priming + eliciting | eliciting, data = WT)
#' summary(fit4)
#' }
#' # Example with different link functions:
#' \donttest{
#' fit5 <- bbreg(agreement ~ priming + eliciting | eliciting,
#'   data = WT,
#'   link.mean = "cloglog", link.precision = "sqrt"
#' )
#' summary(fit5)
#' }
#' 
#' @rdname bbreg
#' @export
bbreg <- function(formula, data, link.mean = c("logit", "probit", "cauchit", "cloglog"),
                  link.precision = c("identity", "log", "sqrt"),
                  model = NULL, residual = NULL, envelope = 0, prob = 0.95, predict = 0, 
                  ptest = 0.25, em_controls = list(maxit = 5000, em_tol = 10^(-5)),
                  optim_method = "L-BFGS-B", optim_controls = list()) {
  ## Processing call
  # If data is not provided, verify the current R workspace
  
  #Processing em_controls:
  if(is.null(em_controls$maxit)){
    em_controls$maxit = 5000
  }
  
  if(is.null(em_controls$em_tol)){
    em_controls$em_tol = 10^(-5)
  }
  
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
      warning("right hand side of formula should not have > 2 parts (ignoring 3rd and higher cases)")
    }
  }
  MF$formula <- Fo
  
  ## Model.frame: Convert MF into a data matrix.
  MF[[1]] <- as.name("model.frame")
  MF <- eval(MF, parent.frame())
  
  ## Extract terms (covariate matrices and response vector)
  MTerms_x <- stats::terms(Fo, data = data, rhs = 1)
  MTerms_v <- stats::delete.response(stats::terms(Fo, data = data, rhs = 2))
  z <- stats::model.response(MF, "numeric")
  x <- stats::model.matrix(MTerms_x, MF)
  v <- stats::model.matrix(MTerms_v, MF)
  
  
  object <- bbreg.fit(z=z,x=x,v=v,link.mean = link.mean, link.precision = link.precision,
                      model = model, residual = residual, envelope = envelope, prob = prob,
                      predict = predict, ptest = ptest, em_controls = em_controls, optim_method = optim_method, optim_controls = optim_controls)
  
  object$call <- Fo
  object$terms <- list(mean = MTerms_x, precision = MTerms_v)
  object$levels <- list(mean = stats::.getXlevels(MTerms_x, MF), precision = stats::.getXlevels(MTerms_v, MF))
  object$contrasts <- list(mean = attr(x, "contrasts"), precision = attr(v, "contrasts"))
  
  class(object) <- "bbreg"
  
  return(object)
  
}


#' @rdname bbreg
#' @export

bbreg.fit <- function(z, x, v = NULL, link.mean = c("logit", "probit", "cauchit", "cloglog"),
                      link.precision = c("identity", "log", "sqrt"),
                      model = NULL, residual = NULL, envelope = 0, prob = 0.95, predict = 0, 
                      ptest = 0.25, em_controls = list(maxit = 5000, em_tol = 10^(-5)), optim_method = "L-BFGS-B", optim_controls = list()) {
  
  #Processing em_controls:
  if(is.null(em_controls$maxit)){
    em_controls$maxit = 5000
  }
  
  if(is.null(em_controls$em_tol)){
    em_controls$em_tol = 10^(-5)
  }
  
  
  n <- length(z)
  x = as.matrix(x)
  if(is.null(v)){
    v = matrix(rep(1,n), nrow = n)
  } else{
    v = as.matrix(v)
  }
  nkap <- ncol(x)
  nlam <- ncol(v)
  
  ## Check for intercepts:
  intercept_x <- sum(colSums(x==1)==n) > 0
  intercept_v <- sum(colSums(v==1)==n) > 0
  intercept <- c(intercept_x, intercept_v)
  
  
  if (nkap == 0) {
    stop("empty matrix x is not allowed")
  }
  if (nlam == 0) {
    stop("empty matrix v is not allowed")
  }
  
  ## Validation of input variables and arguments.
  if (n < 1) {
    stop("response variable is not provided")
  }
  if (min(z) <= 0 & max(z) >= 1) {
    stop("requirement 0 < z_i < 1 is not true for all responses")
  }
  #
  if (is.null(residual) == TRUE) {
    residual <- "pearson"
  }
  if (is.character(residual) == FALSE) {
    stop("residual must be a character (pearson, score or quantile)")
  }
  if (is.character(residual) == TRUE) {
    residual <- match.arg(residual, c("pearson", "score", "quantile"))
  }
  #
  if (is.numeric(envelope) == FALSE) {
    stop("argument envelope must be numeric (0 or positive interger)")
  }
  if (envelope < 0 | envelope %% 1 != 0) {
    stop("consider 0 or positive integers for envelope")
  }
  if (envelope > 0 & envelope < 10) {
    stop("number of simulations to build envelopes is too small (try at least envelope = 10)")
  }
  if (is.numeric(prob) == FALSE) {
    stop("prob must be numeric in the interval (0,1)")
  }
  if (prob <= 0 | prob >= 1) {
    stop("prob must be in the interval (0,1)")
  }
  if (envelope == 0 & prob != 0.95) {
    warning(paste0("prob = ", prob, " is ignored since envelope = 0"))
  }
  #
  if (is.numeric(predict) == FALSE) {
    stop("argument predict must be numeric (0 or positive integer)")
  }
  if (predict < 0 | predict %% 1 != 0) {
    stop("consider 0 or positive integers for predict")
  }
  if (is.numeric(ptest) == FALSE) {
    stop("ptest must be numeric in the interval (0,1)")
  }
  if (ptest <= 0 | ptest >= 1) {
    stop("ptest must be in the interval (0,1)")
  }
  if (predict == 0 & ptest != 0.25) {
    warning(paste0("ptest = ", ptest, " is ignored since predict = 0"))
  }
  #

  ## Set ntest (size of the test set) when "prediction > 0"
  if (predict > 0) {
    ntest <- round(ptest * n)
    if (ntest == 0 | ntest == n) {
      predict == 0
      warning("size of the test set is near 0 or n (RSS_pred is not calculated)")
    }
  }
  
  ## Checking and processing link functions
  if (length(link.mean) > 1) {
    link.mean <- link.mean[1]
  }
  possible_link_mean <- c("logit", "probit", "cauchit", "cloglog")
  
  if (!(link.mean %in% possible_link_mean)) {
    stop(paste0("link function for the mean must be one of ", possible_link_mean))
  }
  
  if (length(link.precision) > 1) {
    if (nlam == 1 & intercept_v) {
      link.precision <- link.precision[1]
    } else {
      link.precision <- link.precision[2]
    }
  }
  possible_link_precision <- c("identity", "log", "sqrt")
  
  if (!(link.precision %in% possible_link_precision)) {
    stop(paste0("link function for the precision parameter must be one of ", possible_link_precision))
  }
  
  ## Set initial values.
  
  ## EM algorithm
  if (is.character(model) == TRUE) {
    aux <- match.arg(model, c("bessel", "beta"))
    if (aux == "bessel") {
      start <- startvalues(z, x, v, link.mean, link.precision, "bessel")
      kap <- start[[1]]
      lam <- start[[2]]
      start <- c(kap, lam)
      #names(start) <- c(paste0("kappa[", 1:nkap, "]"), paste0("lambda[", 1:nlam, "]"))
      names(start) = c(colnames(x), paste(colnames(v),".precision", sep = ""))
      
      modelname <- "Bessel regression"
      message <- paste0(modelname, " via EM - Ignoring the Discrimination test (DBB)")
      inits <- list(kap, lam)
      EM <- EMrun_bes(kap, lam, z, x, v, link.mean, link.precision, em_controls, optim_method, optim_controls)
      
      niter <- EM[[3]] # number of iterations of the EM algorithm
      Est <- EM[[1]] # coefficients
      kap <- Est[1:nkap]
      lam <- Est[-(1:nkap)]
      mu <- as.numeric(EM[[2]][, 1]) # mu
      gphi <- as.numeric(EM[[2]][, 2]) # g(phi)
      if (residual == "pearson") {
        res <- (z - mu) / sqrt(gphi * mu * (1 - mu))
      }
      if (residual == "score") {
        res <- score_residual_bes(kap, lam, z, x, v, link.mean, link.precision)
      }
      if (residual == "quantile") {
        res <- quantile_residual_bes(kap, lam, z, x, v, link.mean, link.precision)
      }
      RSS <- sum(res^2)
      SEr <- infmat_bes(Est, z, x, v, link.mean, link.precision) # standard errors
      DBB <- NULL
      Env <- NULL
      if (envelope > 0) {
        Env <- envelope_bes(residual, Est[1:nkap], Est[-(1:nkap)], x, v, envelope, prob, n, link.mean, link.precision, em_controls, optim_method, optim_controls)
        # % of residuals inside the envelope
        Env_prop <- 100 * sum(sort(res) < Env[1, ] & sort(res) > Env[3, ]) / n
      }
      RSS_pred <- NULL
      if (predict > 0) {
        RSS_pred <- pred_accuracy_bes(residual, kap, lam, z, x, v, ntest, predict, link.mean, link.precision, em_controls, optim_method, optim_controls)
      }
      rm("EM", "Est")
    }
    if (aux == "beta") {
      start <- startvalues(z, x, v, link.mean, link.precision, "beta")
      kap <- start[[1]]
      lam <- start[[2]]
      start <- c(kap, lam)
      names(start) = c(colnames(x), paste(colnames(v),".precision", sep = ""))
      #names(start) <- c(paste0("kappa[", 1:nkap, "]"), paste0("lambda[", 1:nlam, "]"))
      
      modelname <- "Beta regression"
      message <- paste0(modelname, " via EM - Ignoring the Discrimination test (DBB)")
      EM <- EMrun_bet(kap, lam, z, x, v, link.mean, link.precision, em_controls, optim_method, optim_controls)
      niter <- EM[[3]] # number of iterations of the EM algorithm
      Est <- EM[[1]] # coefficients
      kap <- Est[1:nkap]
      lam <- Est[-(1:nkap)]
      mu <- as.numeric(EM[[2]][, 1]) # mu
      gphi <- as.numeric(EM[[2]][, 2]) # g(phi)
      if (residual == "pearson") {
        res <- (z - mu) / sqrt(gphi * mu * (1 - mu))
      }
      if (residual == "score") {
        res <- score_residual_bet(kap, lam, z, x, v, link.mean, link.precision)
      }
      if (residual == "quantile") {
        res <- quantile_residual_bet(kap, lam, z, x, v, link.mean, link.precision)
      }
      RSS <- sum(res^2)
      SEr <- infmat_bet(Est, z, x, v, link.mean, link.precision) # standard errors
      DBB <- NULL
      Env <- NULL
      if (envelope > 0) {
        Env <- envelope_bet(residual, kap, lam, x, v, envelope, prob, n, link.mean, link.precision, em_controls, optim_method, optim_controls)
        # % of residuals inside the envelope
        Env_prop <- 100 * sum(sort(res) < Env[1, ] & sort(res) > Env[3, ]) / n
      }
      RSS_pred <- NULL
      if (predict > 0) {
        RSS_pred <- pred_accuracy_bet(residual, kap, lam, z, x, v, ntest, predict, link.mean, link.precision, em_controls, optim_method, optim_controls)
      }
      rm("EM", "Est")
    }
  }
  if (is.null(model)) {
    # run the discrimination
    aux <- dbbtest.fit(z, x, v, link.mean, link.precision, em_controls, optim_method, optim_controls)
    # fit the chosen model
    if (aux[[2]] == "bessel") {
      start <- startvalues(z, x, v, link.mean, link.precision, "bessel")
      kap <- start[[1]]
      lam <- start[[2]]
      start <- c(kap, lam)
      names(start) = c(colnames(x), paste(colnames(v),".precision", sep = ""))
      #names(start) <- c(paste0("kappa[", 1:nkap, "]"), paste0("lambda[", 1:nlam, "]"))
      
      modelname <- "Bessel regression"
      message <- paste0(modelname, " via EM - Model selected via Discrimination test (DBB)")
      EM <- EMrun_bes(kap, lam, z, x, v, link.mean, link.precision, em_controls, optim_method, optim_controls)
      niter <- EM[[3]] # number of iterations of the EM algorithm
      Est <- EM[[1]] # coefficients
      kap <- Est[1:nkap]
      lam <- Est[-(1:nkap)]
      mu <- as.numeric(EM[[2]][, 1]) # mu
      gphi <- as.numeric(EM[[2]][, 2]) # g(phi)
      if (residual == "pearson") {
        res <- (z - mu) / sqrt(gphi * mu * (1 - mu))
      }
      if (residual == "score") {
        res <- score_residual_bes(kap, lam, z, x, v, link.mean, link.precision)
      }
      if (residual == "quantile") {
        res <- quantile_residual_bes(kap, lam, z, x, v, link.mean, link.precision)
      }
      RSS <- sum(res^2)
      SEr <- infmat_bes(Est, z, x, v, link.mean, link.precision) # standard errors
      DBB <- aux[[1]]
      Env <- NULL
      if (envelope > 0) {
        Env <- envelope_bes(residual, kap, lam, x, v, envelope, prob, n, link.mean, link.precision, em_controls, optim_method, optim_controls)
        # % of residuals inside the envelope
        Env_prop <- 100 * sum(sort(res) < Env[1, ] & sort(res) > Env[3, ]) / n
      }
      RSS_pred <- NULL
      if (predict > 0) {
        RSS_pred <- pred_accuracy_bes(residual, kap, lam, z, x, v, ntest, predict, link.mean, link.precision, em_controls, optim_method, optim_controls)
      }
      rm("EM", "Est")
    }
    if (aux[[2]] == "beta") {
      start <- startvalues(z, x, v, link.mean, link.precision, "beta")
      kap <- start[[1]]
      lam <- start[[2]]
      start <- c(kap, lam)
      names(start) = c(colnames(x), paste(colnames(v),".precision", sep = ""))
      #names(start) <- c(paste0("kappa[", 1:nkap, "]"), paste0("lambda[", 1:nlam, "]"))
      
      modelname <- "Beta regression"
      message <- paste0(modelname, " via EM - Model selected via Discrimination test (DBB)")
      EM <- EMrun_bet(kap, lam, z, x, v, link.mean, link.precision, em_controls, optim_method, optim_controls)
      niter <- EM[[3]] # number of iterations of the EM algorithm
      Est <- EM[[1]] # coefficients
      kap <- Est[1:nkap]
      lam <- Est[-(1:nkap)]
      mu <- as.numeric(EM[[2]][, 1]) # mu
      gphi <- as.numeric(EM[[2]][, 2]) # g(phi)
      if (residual == "pearson") {
        res <- (z - mu) / sqrt(gphi * mu * (1 - mu))
      }
      if (residual == "score") {
        res <- score_residual_bet(kap, lam, z, x, v, link.mean, link.precision)
      }
      if (residual == "quantile") {
        res <- quantile_residual_bet(kap, lam, z, x, v, link.mean, link.precision)
      }
      RSS <- sum(res^2)
      SEr <- infmat_bet(Est, z, x, v, link.mean, link.precision) # standard errors
      DBB <- aux[[1]]
      Env <- NULL
      if (envelope > 0) {
        Env <- envelope_bet(residual, kap, lam, x, v, envelope, prob, n, link.mean, link.precision, em_controls, optim_method, optim_controls)
        # % of residuals inside the envelope
        Env_prop <- 100 * sum(sort(res) < Env[1, ] & sort(res) > Env[3, ]) / n
      }
      RSS_pred <- NULL
      if (predict > 0) {
        RSS_pred <- pred_accuracy_bet(residual, kap, lam, z, x, v, ntest, predict, link.mean, link.precision, em_controls, optim_method, optim_controls)
      }
      rm("EM", "Est")
    }
  }
  
  #Efron pseudo R^2: the squared linear correlation between the response variable and the estimated mean
  efron.pseudo.r2 <- stats::cor(z, mu)^2
  
  #Naming variables accordingly
  names(kap) <- colnames(x)
  names(lam) <- colnames(v)
  names(mu) <- 1:length(mu)
  names(gphi) <- 1:length(gphi)
  
  ###############
  object <- list()
  object$call <- NULL
  object$modelname <- modelname
  object$message <- message
  object$residualname <- residual
  object$niter <- niter
  object$start <- start
  object$intercept <- intercept
  object$link.mean <- link.mean
  object$link.precision <- link.precision
  object$kappa <- kap
  object$lambda <- lam
  object$mu <- mu
  object$efron.pseudo.r2 <- efron.pseudo.r2
  object$fitted.values <- mu
  object$efron.pseudo.r2 <- efron.pseudo.r2
  object$coefficients <- list(mean = kap, precision = lam)
  object$start <- start
  object$x <- x
  object$v <- v
  object$z <- z
  object$gphi <- gphi
  object$residuals <- as.numeric(res)
  object$RSS <- RSS
  object$std_errors <- SEr
  object$DBB <- DBB
  object$envelope <- Env
  if (is.null(Env) == FALSE) {
    object$envelope_prop <- Env_prop
  }
  object$RSS_pred <- RSS_pred
  ###############
  
  class(object) <- "bbreg_fit"
  return(object)
}
