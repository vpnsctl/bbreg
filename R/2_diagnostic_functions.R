## Import key packages
#' @import Formula
#' @import expint
#' @import pbapply
#' @import statmod


#############################################################################################
#' @title plot.bbreg
#' @description Function to build useful plots for bounded regression models.
#' @param x object of class "bbreg" containing results from the fitted model.
#' If the model is fitted with envelope = 0, the Q-Q plot will be produced without envelopes.
#' @param which a number of a vector of numbers between 1 and 4. Plot 1:
#' Residuals vs. Index; Plot 2: Q-Q Plot (if the fit contains simulated envelopes,
#' the plot will be with the simulated envelopes); Plot 3: Fitted means vs. Response;
#' Plot 4: Residuals vs. Fitted means.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param main character; title to be placed at each plot additionally (and above) all captions.
#' @param qqline logical; if \code{TRUE} and the fit does *not* contain simulated
#' envelopes, a qqline will be added to the normal Q-Q plot.
#' @param ... graphical parameters to be passed.
#' @seealso
#' \code{\link{summary.bbreg}}, \code{\link{coef.bbreg}}, \code{\link{vcov.bbreg}}, \code{\link{fitted.bbreg}}, \code{\link{predict.bbreg}}
#' @examples
#' \donttest{
#' n <- 100
#' x <- cbind(rbinom(n, 1, 0.5), runif(n, -1, 1))
#' v <- runif(n, -1, 1)
#' z <- simdata_bes(
#'   kap = c(1, 1, -0.5), lam = c(0.5, -0.5), x, v, repetitions = 1,
#'   link.mean = "logit", link.precision = "log"
#' )
#' z <- unlist(z)
#' fit <- bbreg(z ~ x | v, envelope = 10)
#' plot(fit)
#' plot(fit, which = 2)
#' plot(fit, which = c(1, 4), ask = FALSE)
#' }
#' @export
plot.bbreg <- function(x, which = c(1, 2, 3, 4), ask = TRUE, main = "", qqline = TRUE, ...) {
  if (length(which) == 1) {
    ask <- FALSE
  }
  grDevices::devAskNewPage(ask = ask)
  res <- x$residuals
  call_mod <- deparse(x$call)
  name <- x$modelname
  residualname <- paste0(toupper(substring(x$residualname, 1, 1)), substring(x$residualname, 2))
  residualname <- paste0(residualname, " residuals")
  mu_est <- stats::fitted(x, type = "response")
  
  if (1 %in% which) {
    # First plot (residuals vs index)
    ylab <- residualname
    xlab <- paste0("Index\n bbreg(", call_mod, ")")
    title_1 <- paste0(residualname, " vs Index - ", name)
    graphics::plot(res, xlab = xlab, ylab = ylab, main = main, ...)
    graphics::abline(0, 0, lty = 3)
    graphics::mtext(title_1, side = 3)
  }
  
  # Second plot (QQ-plot)
  if (2 %in% which) {
    env <- x$envelope
    n <- length(res)
    residualname <- paste0(toupper(substring(x$residualname, 1, 1)), substring(x$residualname, 2))
    ylim <- range(res, env)
    xlim <- c(stats::qnorm(0.5 / n), stats::qnorm(1 - 0.5 / n))
    xlab <- paste0("Theoretical quantiles\n bbreg(", call_mod, ")")
    ylab <- paste0(residualname, " residuals")
    if (is.null(env) == FALSE) {
      title_2 <- paste0("Q-Q Plot with simulated envelopes - ", name)
    } else {
      title_2 <- paste0("Q-Q Plot - ", name)
    }
    RR <- stats::qqnorm(res, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, main = main, ...)
    if (is.null(env) == FALSE) {
      aux <- sort(RR$x)
      graphics::lines(aux, env[1, ], col = grDevices::rgb(0.7, 0.7, 0.7))
      graphics::lines(aux, env[3, ], col = grDevices::rgb(0.7, 0.7, 0.7))
      graphics::polygon(c(aux, rev(aux)), c(env[3, ], rev(env[1, ])), col = grDevices::rgb(0.7, 0.7, 0.7), border = NA)
      graphics::lines(aux, env[2, ], lty = 2, lwd = 2)
    } else {
      if (qqline) {
        stats::qqline(res, lty = 3)
      }
    }
    graphics::points(RR$x, RR$y, ...)
    graphics::mtext(title_2, side = 3)
  }
  
  # Third plot (fitted vs response)
  
  if (3 %in% which) {
    obs <- x$z
    
    title_3 <- paste0("Response vs Fitted means - ", name)
    ylab <- "Response"
    xlab <- paste0("Predicted values\n bbreg(", call_mod, ")")
    graphics::plot(mu_est, obs, xlab = xlab, ylab = ylab, main = main, ...)
    graphics::abline(0, 1, lty = 3)
    graphics::mtext(title_3, side = 3)
  }
  
  
  # Fourth plot
  
  if (4 %in% which) {
    title_4 <- paste0(residualname, " vs Fitted means - ", name)
    ylab <- paste0(residualname, " residuals")
    xlab <- paste0("Predicted values\n bbreg(", call_mod, ")")
    graphics::plot(mu_est, res, xlab = xlab, ylab = ylab, main = main, ...)
    graphics::abline(0, 0, lty = 3)
    graphics::mtext(title_4, side = 3)
  }
}


#############################################################################################
#' @title fitted.bbreg
#' @description Function providing the fitted means for the model (bessel or beta).
#' @param object object of class "bbreg" containing results from the fitted model.
#' @param type the type of variable to get the fitted values. The default is the "response" type, which provided the estimated values for the means. The type "link" provides the estimates for the linear predictor of the mean. The type "precision" provides estimates for the precision parameters whereas the type "variance" provides estimates for the variances.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{predict.bbreg}}, \code{\link{summary.bbreg}}, \code{\link{coef.bbreg}}, \code{\link{vcov.bbreg}}, \code{\link{plot.bbreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting, data = WT)
#' fitted(fit)
#' fitted(fit, type = "precision")
#' }
#' @export
fitted.bbreg <- function(object, type = c("response", "link", "precision", "variance"), ...) {
  fit <- object
  if (length(type) > 1) {
    type <- type[1]
  }
  possible_types <- c("response", "link", "precision", "variance")
  if (!(type %in% c("response", "link", "precision", "variance"))) {
    stop(paste0("type must be one of ", possible_types))
  }
  
  link_precision <- stats::make.link(fit$link.precision)
  
  fitted_values <- switch(type,
                          "response" = {
                            mu <- fit$mu
                            names(mu) <- 1:length(mu)
                            mu
                          },
                          "link" = {
                            link_fitted <- c(fit$x %*% fit$kappa)
                            names(link_fitted) <- 1:length(link_fitted)
                            link_fitted
                          },
                          "precision" = {
                            fitted_prec <- c(link_precision$linkinv(fit$v %*% fit$lambda))
                            names(fitted_prec) <- 1:length(fitted_prec)
                            fitted_prec
                          },
                          "variance" = {
                            variance_fitted <- c(fit$mu * (1 - fit$mu) * fit$gphi)
                            names(variance_fitted) <- 1:length(variance_fitted)
                            variance_fitted
                          }
  )
  return(fitted_values)
}


#############################################################################################
#' @title predict.bbreg
#' @description Function to obtain various predictions based on the fitted model (bessel or beta).
#' @param object object of class "bbreg" containing results from the fitted model.
#' @param newdata optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted response values will be provided.
#' @param type the type of prediction. The default is the "response" type, which provided the estimated values for the means. The type "link" provides the estimates for the linear predictor. The type "precision" provides estimates for the precision parameters whereas the type "variance" provides estimates for the variances.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{fitted.bbreg}}, \code{\link{summary.bbreg}}, \code{\link{coef.bbreg}}, \code{\link{vcov.bbreg}}, \code{\link{plot.bbreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting, data = WT)
#' predict(fit)
#' new_data_example <- data.frame(priming = c(0, 0, 1), eliciting = c(0, 1, 1))
#' predict(fit, new_data = new_data_example)
#' predict(fit, new_data = new_data_example, type = "precision")
#' }
#' @export
predict.bbreg <- function(object, newdata = NULL, type = c("response", "link", "precision", "variance"), ...) {
  fit <- object
  if (length(type) > 1) {
    type <- type[1]
  }
  
  possible_types <- c("response", "link", "precision", "variance")
  if (!(type %in% c("response", "link", "precision", "variance"))) {
    stop(paste0("type must be one of ", possible_types))
  }
  
  
  if (missing(newdata)) {
    predictions <- stats::fitted(fit, type)
  } else {
    formula_temp <- Formula(fit$call)
    matrix_temp_x <- stats::model.matrix(object = formula_temp, data = newdata, rhs = 1)
    matrix_temp_v <- stats::model.matrix(object = formula_temp, data = newdata, rhs = 2)
    
    kappa_est <- fit$kappa
    lambda_est <- fit$lambda
    
    link_mean <- stats::make.link(fit$link.mean)
    link_precision <- stats::make.link(fit$link.precision)
    
    mu_est <- link_mean$linkinv(matrix_temp_x %*% kappa_est)
    mu_est <- c(mu_est)
    names(mu_est) <- 1:length(mu_est)
    
    phi_est <- c(link_precision$linkinv(matrix_temp_v %*% fit$lambda))
    names(phi_est) <- 1:length(phi_est)
    
    if (fit$modelname == "Bessel regression") {
      gphi_est <- c((1 - phi_est + (phi_est^2) * exp(phi_est) * expint_En(phi_est, order = 1)) / 2)
    } else {
      gphi_est <- c(1 / (1 + phi_est))
    }
    
    predictions <- switch(type,
                          "response" = {
                            mu_est
                          },
                          "link" = {
                            link_predict <- c(matrix_temp_x %*% kappa_est)
                            names(link_predict) <- 1:length(link_predict)
                            link_predict
                          },
                          "precision" = {
                            phi_est
                          },
                          "variance" = {
                            variance_fitted <- c(mu_est * (1 - mu_est) * gphi_est)
                            names(variance_fitted) <- 1:length(variance_fitted)
                            variance_fitted
                          }
    )
  }
  return(predictions)
}


#############################################################################################
#' @title print.bbreg
#' @description Function providing a brief description of results related to the regression model (bessel or beta).
#' @param x object of class "bbreg" containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{fitted.bbreg}}, \code{\link{summary.bbreg}}, \code{\link{coef.bbreg}}, \code{\link{vcov.bbreg}}, \code{\link{plot.bbreg}}, \code{\link{predict.bbreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting, data = WT)
#' fit
#' }
#' @export
print.bbreg <- function(x, ...) {
  nkap <- length(x$kappa)
  nlam <- length(x$lambda)
  #
  coeff_kappa <- x$kappa
  names(coeff_kappa) <- colnames(x$x)
  coeff_lambda <- x$lambda
  names(coeff_lambda) <- colnames(x$v)
  cat("\n")
  cat(x$message)
  cat("\n\n")
  cat("Call:", "\n")
  call_model <- deparse(x$call)
  cat(paste0("bbreg(", call_model, ")"))
  cat("\n\n")
  cat(paste0("Coefficients modeling the mean (with ", x$link.mean, " link):", "\n"))
  print(coeff_kappa)
  cat(paste0("Coefficients modeling the precision (with ", x$link.precision, " link):", "\n"))
  print(coeff_lambda)
}


#############################################################################################
#' @title coef.bbreg
#' @description Function to extract the coefficients of a fitted regression model (bessel or beta).
#' @param object object of class "bbreg" containing results from the fitted model.
#' @param parameters a string to determine which coefficients should be extracted: 'all' extracts all coefficients, 'mean' extracts the coefficients of the mean parameters and 'precision' extracts coefficients of the precision parameters.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{fitted.bbreg}}, \code{\link{summary.bbreg}}, \code{\link{vcov.bbreg}}, \code{\link{plot.bbreg}}, \code{\link{predict.bbreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting, data = WT)
#' coef(fit)
#' coef(fit, parameters = "precision")
#' }
#' @export
coef.bbreg <- function(object, parameters = c("all", "mean", "precision"), ...) {
  if(length(parameters)>1){
    parameters = parameters[1]
  }
  coef_ext <- switch(parameters,
                     "all" = {
                       object$coefficients
                     }, "mean" = {
                       object$coefficients$mean
                     }, "precision" = {
                       object$coefficients$precision
                     }
  )
  return(coef_ext)
}

#############################################################################################
#' @title vcov.bbreg
#' @description Function to extract the variance-covariance matrix of the parameters of the fitted regression model (bessel or beta).
#' @param object an object of class "bbreg" containing results from the fitted model.
#' @param parameters a string to determine which coefficients should be extracted: 'all' extracts all coefficients, 'mean' extracts the coefficients of the mean parameters and 'precision' extracts coefficients of the precision parameters.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{infmat_bes}}, \code{\link{infmat_bet}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting | priming, data = WT)
#' vcov(fit)
#' vcov(fit, parameters = "precision")
#' }
#' @export
vcov.bbreg <- function(object, parameters = c("all", "mean", "precision"), ...) {
  fit <- object
  nkap <- length(fit$kappa)
  nlam <- length(fit$lambda)
  #
  coeff_kappa <- fit$coefficients$mean
  coeff_lambda <- fit$coefficients$precision
  all_coeff <- c(coeff_kappa, coeff_lambda)
  coeff_names = c(names(coeff_kappa),paste(names(coeff_lambda), ".precision", sep = ""))
  
  if (length(parameters) > 1) {
    parameters <- parameters[1]
  }
  model_bb <- fit$modelname
  vcov_bb_complete <- switch(model_bb,
                             "Bessel regression" = {
                               infmat_bes(all_coeff, fit$z, fit$x, fit$v, fit$link.mean, fit$link.precision, information = TRUE)
                             },
                             "Beta regression" = {
                               infmat_bet(all_coeff, fit$z, fit$x, fit$v, fit$link.mean, fit$link.precision, information = TRUE)
                             }
  )
  vcov_bb_complete <- tryCatch(solve(vcov_bb_complete), error = function(e) rep(NA, nrow(vcov_bb_complete)))
  vcov_bb <- switch(parameters,
                    "all" = {
                      colnames(vcov_bb_complete) <- rownames(vcov_bb_complete) <- coeff_names
                      vcov_bb_complete
                    },
                    "mean" = {
                      vcov_bb_complete[1:nkap, 1:nkap]
                    },
                    "precision" = {
                      vcov_bb_complete[(nkap + 1):(nkap + nlam), (nkap + 1):(nkap + nlam)]
                    }
  )
  
  return(vcov_bb)
}




#############################################################################################
#' @title summary.bbreg
#' @description Function providing a summary of results related to the regression model (bessel or beta).
#' @param object an object of class "bbreg" containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @seealso
#' \code{\link{fitted.bbreg}}, \code{\link{plot.bbreg}}, \code{\link{predict.bbreg}}
#' @examples
#' \donttest{
#' fit <- bbreg(agreement ~ priming + eliciting | priming, data = WT)
#' summary(fit)
#' }
#' @export
summary.bbreg <- function(object, ...) {
  M <- object
  nkap <- length(M$kappa)
  nlam <- length(M$lambda)
  #
  Est <- c(M$kappa, M$lambda)
  SEr <- M$std_errors
  tab <- cbind(Est, SEr, Est / SEr, 2 * stats::pnorm(-abs(Est / SEr)))
  colnames(tab) <- c("Estimate", "Std.error", "z-value", "Pr(>|z|)")
  rownames(tab) <- c(colnames(M$x), colnames(M$v))
  tab <- list(mean = tab[seq.int(length.out = nkap), , drop = FALSE], precision = tab[seq.int(length.out = nlam) + nkap, , drop = FALSE])
  #
  digits <- max(3, getOption("digits") - 3)
  #
  message(sprintf("\n%s:\n", M$message))
  message("Call:\n", paste0("bbreg(", deparse(M$call, width.cutoff = floor(getOption("width") * 0.85)), ")"), "", sep = "\n")
  message(sprintf("%s\n", paste0("Number of iterations of the EM algorithm = ", M$niter)))
  #
  if (is.null(M$DBB) == FALSE) {
    cat(sprintf("\n %s:\n", "Results of the discrimination test DBB"))
    print(structure(round(M$DBB, digits = digits)))
  }
  #
  residualname <- paste0(toupper(substring(M$residualname, 1, 1)), substring(M$residualname, 2))
  cat(sprintf("\n %s:\n", paste0(residualname, " residuals")))
  print(structure(round(c(M$RSS, as.vector(stats::quantile(M$residuals))), digits = digits), .Names = c("RSS", "Min", "1Q", "Median", "3Q", "Max")))
  #
  if (NROW(tab$mean)) {
    cat(paste0("\n Coefficients modeling the mean (with ", M$link.mean, " link):\n"))
    stats::printCoefmat(tab$mean, digits = digits, signif.legend = FALSE)
  } else {
    message("\n No coefficients modeling the mean. \n")
  }
  #
  if (NROW(tab$precision)) {
    cat(paste0("\n Coefficients modeling the precision (with ", M$link.precision, " link):\n"))
    stats::printCoefmat(tab$precision, digits = digits, signif.legend = FALSE)
  } else {
    message("\n No coefficients modeling the precision. \n")
  }
  #
  gp <- unique(M$gphi)
  if (length(gp) == 1) {
    cat(sprintf("%s\n", paste0("g(phi) = ", round(gp, digits = digits))))
  }
  #
  if (getOption("show.signif.stars")) {
    cat("---\n Signif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n\n")
  }
  #
  if (is.null(M$RSS_pred) == FALSE) {
    message(sprintf("%s\n", paste0("Average RSS from the predictive accuracy simulation = ", round(mean(M$RSS_pred), digits = digits))))
  }
  #
  if (is.null(M$envelope) == FALSE) {
    message(sprintf("%s\n", paste0("Percentage of residual within the envelope = ", round(M$envelope_prop, digits = digits))))
  }
}


