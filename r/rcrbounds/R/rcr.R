#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


# global reference to rcrpy (will be initialized in .onLoad)
rcrpy <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  rcrpy <<- reticulate::import("rcrbounds", delay_load = TRUE)
}


#' Install rcrbounds Python package
#'
#' The rcrbounds R package is a wrapper to a Python
#' package, also called "rcrbounds". This function
#' installs the required Python package from the
#' Python Package Index (PyPI), and may need to be
#' run once before using the functions in the R
#' package.
#'
#' For most users, this function can be run once with the default
#' arguments and then ignored.
#'
#' @inheritParams reticulate::py_install
#' @seealso [rcr()]
#'          [reticulate::py_install()] which this function wraps.
#' @export
install_rcrpy <- function(method = "auto", conda = "auto") {
  reticulate::py_install("rcrbounds",
    pip = TRUE,
    method = method,
    conda = conda
  )
}

#' Estimate an RCR model
#'
#' `rcr()` estimates the linear casual effect of one explanatory
#' variable on an outcome, under a relative correlation restriction
#' of the form described in Krauth (2016).
#'
#' Most of the calculation is done in an external Python function which
#' is called by `rcr()`. You may need to use the [install_rcrpy()] function
#' to install the package the first time you wish to use this function.
#'
#' @param formula An object of class "formula" (or one that
#'        can be coerced to that class): a symbolic description
#'        of the model to be fitted. The details of model specification
#'        are given under 'Details'.
#' @param data an optional data frame, list or environment (or object
#'        coercible by as.data.frame to a data frame) containing the
#'        variables in the model. If not found in data, the variables
#'        are taken from environment(formula), typically the environment
#'        from which rcr is called.
#' @param subset an optional vector specifying a subset of observations
#'        to be used in the fitting process.
#' @param weights an optional vector of weights to be used in the fitting
#'        process. Should be NULL or a numeric vector.
#' @param na.action a function which indicates what should happen when
#'        the data contain NAs. The default is set by the na.action setting
#'        of options, and is na.fail if that is unset. The ‘factory-fresh’
#'        default is na.omit. Another possible value is NULL, no action.
#'        Value na.exclude can be useful.
#' @param model description goes here.
#' @param cluster An optional vector defining groups for cluster-robust
#'        covariance matrix estimation. Should be NULL or a numeric
#'        vector. See also ‘Details’,
#' @param rc_range An optional numeric vector of length two indicating
#'        the assumed range for the relative correlation parameter
#'        lambda.`-Inf` is allowed for the lower bound, and `Inf` is
#'        allowed for the upper bound.
#' @param vceadj An optional adjustment factor to perform degrees-of-freedom
#'        adjustments for the estimated covariance matrix.  That is,
#'        the estimated covariance matrix will be multiplied by `vceadj`.
#' @returns `rcr()` returns an object of class "rcr", which is
#' a list containing the following components:
#' \describe{
#'   \item{`coefficients`}{a named vector of coefficients}
#'   \item{`cov.unscaled`}{the estimated covariance matrix for `coefficients`}
#'   \item{`na.action`}{(where relevant) information returned by model.frame on the special handling of NAs..}
#'   \item{`xlevels`}{(where relevant) a record of the levels of the factors used in fitting..}
#'   \item{`call`}{the matched call.}
#'   \item{`terms`}{the terms object used.}
#'   \item{`model`}{if requested (the default), the model frame used.}
#'   \item{`weights`}{(where relevant) the specified weights.}
#'   \item{`cluster`}{(where relevant) the specified cluster identifier.}
#'   \item{`pyobj`}{a pointer to the Python object returned
#'          by [rcr.fit()]. }
#' }
#' @seealso [rcr.fit()],
#' @examples
#' f1 <- 1
#' @export
rcr <- function(formula,
                data,
                subset,
                weights,
                na.action,
                model = TRUE,
                cluster = NULL,
                rc_range = c(0, 1),
                vceadj = 1.0) {
  contrasts <- NULL
  cl <- match.call()
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  # offset and contrasts are not arguments to the function but are kept
  # here to keep the code just like the code for lm.
  m <- match(c(
    "formula", "data", "subset", "weights", "na.action",
    "offset","cluster"
  ), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  formula <- Formula::as.Formula(formula)
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- stats::model.response(mf, "numeric")
  mt <- stats::terms(formula, data = data)
  mtX <- stats::terms(stats::update(formula, ~ . - 1), data = data, rhs = 1)
  X <- stats::model.matrix(mtX, mf, contrasts)
  mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2))
  exog <- stats::model.matrix(mtZ, mf, contrasts)
  endog <- cbind(Y, X)
  weights <- stats::model.weights(mf)
  if (is.null(cluster)){
    cov_type <- "nonrobust"
  } else {
    cov_type <- "cluster"
  }
  result <- list()
  pyobj <- rcr.fit(endog, exog,
    vceadj = vceadj,
    rc_range = rc_range,
    weights = weights,
    cov_type = cov_type,
    groupvar = cluster)
  result$coefficients <- pyobj$params
  names(result$coefficients) <- pyobj$param_names
  result$cov.unscaled <- pyobj$cov_params
  rownames(result$cov.unscaled) <- pyobj$param_names
  colnames(result$cov.unscaled) <- pyobj$param_names
  result$na.action <- attr(mf, "na.action")
  result$xlevels <- stats::.getXlevels(mt, mf)
  result$call <- cl
  result$terms <- mt
  if (model) {
    result$model <- mf
  }
  result$weights <- weights
  result$cluster <- cluster
  result$pyobj <- pyobj
  class(result) <- c("rcr")
  result
}

#' Fitting function for RCR model
#'
#' `rcr.fit()` is a bare-bones wrapper function to call the Python
#' function that estimates the RCR model. It is not intended for
#' most users; you should use [rcr()] unless you wish to interact
#' more directly with the Python code.
#'
#' The object returned by `rcr.fit()` is a pointer to a Python object of
#' class "`rcrbounds.RCRResults`" and not a normal R object. All attributes
#' *and methods* of the Python object are available within R. See the
#' documentation for the Python rcrbounds module for details on these
#' attributes and methods.
#'
#' Note that the Python object ceases to exist when the R session ends,
#' and cannot be saved in .RData or any other file. The pointer
#' to the Python object then becomes a NULL pointer.
#'
#' @inheritParams rcr
#' @param endog a matrix or `data.frame` representing the
#'              endogenous (outcome and treatment) variables
#'              in the model
#' @param exog a matix or `data.frame` representing
#'              the exogenous (control) variables in the model.
#'              Its first column must be an intercept.
#' @param groupvar an optional cluster ID variable for
#'                 cluster-robust covariance matrix estimates,
#'                 or NULL
#' @param cov_type the type of covariance matrix to
#'                 estimate, either "`nonrobust`" (the
#'                 default) or "`cluster`"
#' @param citype the confidence interval type, either
#'               "`conservative`" (the default), "`upper`",
#'               "`lower`", or "`Imbens-Manski`"
#' @param cilevel the desired confidence level,
#'                on a 0-100 scale.
#' @returns `rcr.fit()` returns a pointer to the Python object of class
#'          "`rcrbounds.RCRResults`" as returned by the Python call.
#' @seealso [rcr()]
#' @export
rcr.fit <- function(endog,
                    exog,
                    rc_range = c(0.0, 1.0),
                    weights = NULL,
                    groupvar = NULL,
                    cov_type = "nonrobust",
                    vceadj = 1.0,
                    citype = "conservative",
                    cilevel = 95) {
  if (!is.null(weights)) {
    weights <- reticulate::np_array(weights)
  }
  rcrpy$RCR(
    endog,
    exog,
    rc_range,
    cov_type,
    vceadj,
    citype,
    cilevel,
    weights,
    groupvar
  )$fit()
}

#' Print method for RCR results
#'
#' `print.rcr()` prints the main results of estimating an RCR model
#'
#' @param x an object of class "`rcr`", usually the result of a call
#'        to `rcr`
#' @param ... Additional arguments to be passed on to `print()`,
#'            for example `digits=`.
#' @returns `print.rcr()` returns its input object `x` invisibly.
#' @seealso [rcr()], [print.summary.rcr()]
#' @export
print.rcr <- function(x,...){
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(stats::coefficients(x,...),
    ...)
  cat("\n")
  invisible(x)
}

#' Calculate covariance matrix for RCR coefficients
#'
#' `vcov()` Returns the variance-covariance matrix of the main parameters
#'  of the rcr model object.
#'
#' @param object an object of class "rcr", usually the result of a call
#'        to `rcr`.
#' @param ... Additional arguments (not used).
#' @returns A symmetric 5-by-5 matrix with row and column names corresponding
#'          to the parameter names in `coef(object)`.
#' @export
vcov.rcr <- function(object, ...) {
  object$cov.unscaled
}

#' Confidence intervals for RCR results
#'
#' Computes confidence intervals for one or more parameters in an RCR model.
#'
#' @param object an object of class "rcr", usually the result of a call
#'        to `rcr`.
#' @param parm a specification of which parameters are to be given
#'        confidence intervals, either a vector of numbers or a vector
#'        of names. Available parameters include "`rcInf`", "`effectInf`",
#'        "`rc0`", "`effectL`", "`effectH`", and "`effect`".
#'        If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param citype the confidence interval type: "`conservative`" (the default)
#'        "`upper`", "`lower`" or "`Imbens-Manski`".
#' @param ... additional arguments (not used).
#' @returns `confint.rcr()` returns a matrix with columns giving
#'          upper and lower confidence limits for each parameter
#' @seealso [rcr()]
#' @export
confint.rcr <- function(object,
                        parm,
                        level = 0.95,
                        citype = "conservative",
                        ...) {
  cf <- object$coefficients
  ses <- sqrt(diag(object$cov.unscaled))
  pnames <- names(ses)
  if (missing(parm)) {
    parm <- c(pnames, "effect")
  }
  a <- NULL
  if (citype == "conservative" | citype == "Imbens-Manski") {
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
  } else {
    if (citype == "lower") {
      a <- c(0, level)
    }
    if (citype == "upper") {
      a <- c(1 - level, 1)
    }
  }
  if (is.null(a)){
    msg <- paste("invalid citype",citype)
    stop(msg)
  }
  fac <- stats::qnorm(a)
  pct <- paste(format(100 * a, 3), " %")
  ci <- array(NA_real_,
    dim = c(length(parm), 2L),
    dimnames = list(parm, pct)
  )
  ci[] <- cbind(cf[parm], cf[parm]) + ses[parm] %o% fac
  effects <- (parm == "effect")
  if (any(effects)) {
    ci[effects, ] <- matrix(effect_ci(object,
      level = level,
      citype = citype
    ),
    nrow = sum(effects),
    ncol = 2,
    byrow = TRUE
    )
  }
  ci
}

#' Summary method for RCR results
#'
#' `summary.rcr()` produces a tabular summary of rcr model
#' results estimated by the [rcr()] function.
#'
#' @inheritParams confint.rcr
#' @param x an object of class "summary.rcr", usually the result
#'        of a call to `summary.rcr`.
#' @param ... Additional arguments (e.g., `digits`) to be
#'        passed to `print()`
#' @returns `summary.rcr()` returns an rcr.summary object, which
#'          is a list with the following components:
#' @seealso [rcr()].
#' @export
summary.rcr <- function(object,
                        level = 0.95,
                        citype = "conservative",
                        ...) {
  cf <- stats::coef(object)
  ses <- sqrt(diag(stats::vcov(object)))
  tvalue <- cf / ses
  pvalue <- 2 * stats::pnorm(-abs(tvalue))
  result <- list()
  result$call <- object$call
  result$terms <- object$terms
  result$coefficients <- cbind(cf, ses, tvalue, pvalue)
  colnames(result$coefficients) <- c(
    "Estimate",
    "Std. Error",
    "t value",
    "Pr(>|t|)"
  )
  result$cov.unscaled <- object$cov.unscaled
  result$effect_ci <- confint.rcr(object,
    "effect",
    citype = citype,
    level = level
  )
  result$citype <- citype
  result$level <- level
  class(result) <- "summary.rcr"
  result
}

#' @rdname summary.rcr
#' @export
print.summary.rcr <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(stats::coefficients(x), ...)
  cat("---\n")
  cat(
    x$citype,
    "confidence interval:\n"
  )
  print(x$effect_ci, ...)
  invisible(x)
}

effect_ci <- function(object,
                      citype = "conservative",
                      level = 0.95) {
  cf <- object$coefficients[4:5]
  vcf <- object$cov.unscaled[4:5, 4:5]
  ses <- sqrt(diag(vcf))
  ci <- NULL
  if (citype == "conservative") {
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    fac <- stats::qnorm(a)
    ci <- cf + fac * ses
  }
  if (citype == "lower") {
    a <- c(0, level)
    fac <- stats::qnorm(a)
    ci <- cf + fac * ses
  }
  if (citype == "upper") {
    a <- c(1 - level, 1)
    fac <- stats::qnorm(a)
    ci <- cf + fac * ses
  }
  if (citype == "Imbens-Manski") {
    cv_min <- stats::qnorm(1 - (1 - level))
    cv_max <- stats::qnorm(1 - (1 - level) / 2)
    cv_mid <- cv_min
    delta <- (cf[2] - cf[1]) / max(ses)
    if (is.finite(delta)) {
      while ((cv_max - cv_min) > 0.000001) {
        cv_mid <- (cv_min + cv_max) / 2.0
        if (stats::pnorm(cv_mid + delta) -
          stats::pnorm(-cv_mid) < level) {
          cv_min <- cv_mid
        } else {
          cv_max <- cv_mid
        }
      }
    }
    ci <- cf
    if (ses[1] > 0.0) {
      ci[1] <- cf[1] - cv_mid * ses[1]
    } else {
      ci[1] <- -Inf
    }
    if (ses[2] > 0.0) {
      ci[2] <- cf[2] + cv_mid * ses[2]
    } else {
      ci[2] <- Inf
    }
  }
  if (is.null(ci)) {
    msg <- paste("Invalid citype:", citype)
    stop(msg)
  }
  ci
}

#' Hypothesis testing for RCR causal effect
#'
#' `effect_test()` allows a user to test a point null hypothesis on the
#' (interval identified) causal effect parameter in the RCR model.
#'
#' The hypothesis test is based on inverting the Imbens-Manski confidence
#' interval rather than on a specific test statistic or set of critical
#' values.  For example, we reject the null hypothesis at 5% significance
#' if the value under the null is outside of the 95% Imbens-Manski
#' confidence interval.
#'
#' @param object an object of class "rcr", usually the result of a call
#'        to `rcr`.
#' @param h0 the value of the parameter under the null hypothesis.  Default
#'           is zero.
#' @returns `effect_test()` returns the p-value associated with the
#'          specified null hypothesis.
#' @seealso [rcr()], [confint.rcr()]
#' @export
effect_test <- function(object, h0 = 0.0) {
  pmin <- 0.0
  pmid <- 0.5
  pmax <- 1.0
  if (h0 >= object$coefficients[4] & h0 <= object$coefficients[5]) {
    pvalue <- 1.0
  } else {
    while (pmax - pmin > 0.0000001) {
      pmid <- (pmax + pmin)/2.0
      ci <- effect_ci(object,
                      citype = "Imbens-Manski",
                      level = pmid
      )
      if (h0 >= ci[1] & h0 <= ci[2]) {
        pmax <- pmid
      } else {
        pmin <- pmid
      }
    }
    pvalue <- 1.0 - pmin
  }
  pvalue
}
