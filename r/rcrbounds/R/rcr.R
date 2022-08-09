#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


# global reference to rcrpy (will be initialized in .onLoad)
rcrpy <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to rcrbounds
  rcrpy <<- reticulate::import("rcrbounds", delay_load = TRUE)
}


#' Install rcrbounds Python package
#'
#' `install_rcrpy()` is a utility function to install the rcrbounds
#' Python package from the Python Package Index (PyPI). This function
#' may need to be run once before using the rcrbounds R package.
#' For most users, this function can be run once with the default
#' arguments and then ignored.
#'
#' @inheritParams reticulate::py_install
#' @returns `install_rcrpy()` returns whatever is returned by the
#'           call to [reticulate::py_install], usually `NULL`.
#' @seealso [rcr()] to estimate an RCR model,
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
#' `rcr()` estimates the linear casual effect of a scalar explanatory
#' variable on a scalar outcome, under a relative correlation restriction
#' of the form described in Krauth (2016).
#'
#' Most of the calculations are done in an external Python function which
#' is called by `rcr()`. You may need to call the [install_rcrpy()] function
#' to install the package the first time you wish to use this function.
#'
#' The `formula` argument should be in the form
#' `outcome ~ treatment | c1 + c2 + c3`,
#' where `outcome` is the outcome variable, `treatment` is
#' the treatment variable, and `c1`, `c2`, etc. are the control
#' variables.
#'
#' @param formula An object of class "`formula`" (or one that
#'        can be coerced to that class): a symbolic description
#'        of the model to be fitted. The details of model specification
#'        are given under 'Details'.
#' @param data An optional `data frame`, `list` or environment (or object
#'        coercible by `as.data.frame` to a data frame) containing the
#'        variables in the model. If not found in data, the variables
#'        are taken from `environment(formula)`, typically the environment
#'        from which `rcr` is called.
#' @param rc_range An optional numeric vector of length two indicating
#'        the assumed range for the relative correlation parameter
#'        lambda.`-Inf` is allowed for the lower bound, and `Inf` is
#'        allowed for the upper bound.
#' @param subset An optional vector specifying a subset of observations
#'        to be used in the fitting process.
#' @param weights An optional vector of weights to be used in the fitting
#'        process.
#' @param cluster An optional vector defining groups for cluster-robust
#'        covariance matrix estimation.
#' @param na.action An optional function which indicates what should happen when
#'        the data contain `NA`s. The default is set by the `na.action` setting
#'        of options, and is `na.fail` if that is unset. The ‘factory-fresh’
#'        default is `na.omit`. Another possible value is `NULL`, no action.
#'        Value `na.exclude` can be useful.
#' @param model,pyobj Optional logicals. If `TRUE`, the corresponding
#'        intermediate components of the fit (the model frame and the Python
#'        object) are returned.
#' @param vceadj An optional adjustment factor to perform degrees-of-freedom
#'        adjustments for the estimated covariance matrix.  That is,
#'        the estimated covariance matrix will be multiplied by the value of
#'        `vceadj`.
#' @returns `rcr()` returns an object of class "`rcr`", which is
#' a list containing the following components:
#' \describe{
#'   \item{`coefficients`}{a named vector of coefficients}
#'   \item{`cov.unscaled`}{the estimated covariance matrix for `coefficients`}
#'   \item{`na.action`}{(where relevant) information returned by model.frame on
#'     the special handling of NAs..}
#'   \item{`xlevels`}{(where relevant) a record of the levels of the factors
#'     used in fitting..}
#'   \item{`call`}{the matched call.}
#'   \item{`terms`}{the terms object used.}
#'   \item{`model`}{if requested (the default), the model frame used.}
#'   \item{`weights`}{(where relevant) the specified weights.}
#'   \item{`cluster`}{(where relevant) the specified cluster identifier.}
#'   \item{`pyobj`}{if requested (the default), a pointer to the Python object
#'          returned by [rcr_fit()]. }
#' }
#' @seealso [rcr_fit()] the lower-level model fitting command,
#'          [coefficients()] to retrieve coefficient estimates,
#'          [vcov.rcr()] to retrieve the covariance matrix,
#'          [summary.rcr()] to produce a table summarizing results,
#'          [confint.rcr()] to produce confidence intervals,
#'          [effect_test()] to test null hypotheses, and
#'          [plot.rcr()] to plot the results.
#' @references Krauth, B. V. (2016). "Bounding a linear causal effect using
#'             relative correlation restrictions"
#'             *Journal of Econometric Methods* 5(1): 117-141.
#' @examplesIf reticulate::py_module_available("rcrbounds")
#' # A simple example with default options
#' rcr(weight ~ Time | Diet, ChickWeight)
#' # Use rc_range to change the range of values for the
#' # relative correlation (lambda) parameter
#' rcr(weight ~ Time | Diet, ChickWeight, rc_range = c(0, 2))
#' @export
rcr <- function(formula,
                data,
                rc_range = c(0, 1),
                subset,
                weights,
                cluster = NULL,
                na.action,
                model = TRUE,
                pyobj = TRUE,
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
    "offset", "cluster"
  ), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  formula <- Formula::as.Formula(formula)
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  y <- stats::model.response(mf, "numeric")
  mt <- stats::terms(formula, data = data)
  mt_x <- stats::terms(stats::update(formula, ~ . - 1), data = data, rhs = 1)
  x <- stats::model.matrix(mt_x, mf, contrasts)
  mt_z <- stats::delete.response(stats::terms(formula, data = data, rhs = 2))
  exog <- stats::model.matrix(mt_z, mf, contrasts)
  endog <- cbind(y, x)
  weights <- stats::model.weights(mf)
  if (is.null(cluster)) {
    cov_type <- "nonrobust"
  } else {
    cov_type <- "cluster"
  }
  result <- list()
  python_result <- rcr_fit(endog, exog,
    vceadj = vceadj,
    rc_range = rc_range,
    weights = weights,
    cov_type = cov_type,
    groupvar = cluster
  )
  result$coefficients <- python_result$params
  names(result$coefficients) <- python_result$param_names
  result$cov.unscaled <- python_result$cov_params
  rownames(result$cov.unscaled) <- python_result$param_names
  colnames(result$cov.unscaled) <- python_result$param_names
  result$na.action <- attr(mf, "na.action")
  result$xlevels <- stats::.getXlevels(mt, mf)
  result$call <- cl
  result$terms <- mt
  if (model) {
    result$model <- mf
  }
  if (pyobj) {
    result$pyobj <- python_result
  }
  result$weights <- weights
  result$cluster <- cluster
  class(result) <- c("rcr")
  result
}

#' Fitting function for RCR model
#'
#' `rcr_fit()` is a bare-bones wrapper function to call the Python
#' code to estimate the RCR model. It is not intended for
#' most users; you should use [rcr()] unless you wish to interact
#' more directly with the Python code.
#'
#' The object returned by `rcr_fit()` is a pointer to a Python object of
#' class "`rcrbounds.RCRResults`" and not a normal R object. All attributes
#' *and methods* of the Python object are available within R. See the
#' documentation for the Python `rcrbounds` module for details on these
#' attributes and methods.
#'
#' Note that the Python object ceases to exist when the R session ends,
#' and cannot be saved in `.RData` or any other file. The pointer
#' to the Python object then becomes a null pointer.
#'
#' @inheritParams rcr
#' @param endog A matrix or `data.frame` representing the
#'              endogenous (outcome and treatment) variables
#'              in the model
#' @param exog A matix or `data.frame` representing
#'              the exogenous (control) variables in the model.
#'              Its first column must be an intercept.
#' @param groupvar An optional cluster ID variable for
#'                 cluster-robust covariance matrix estimates,
#'                 or NULL
#' @param cov_type The type of covariance matrix to
#'                 estimate, either "`nonrobust`" (the
#'                 default) or "`cluster`"
#' @param citype The confidence interval type, either
#'               "`conservative`" (the default), "`upper`",
#'               "`lower`", or "`Imbens-Manski`"
#' @param cilevel The desired confidence level,
#'                on a 0-100 scale.
#' @returns `rcr_fit()` returns a pointer to the Python object of class
#'          "`rcrbounds.RCRResults`" as returned by the Python call.
#' @seealso [rcr()] the high-level version of this function which
#'          should be used by most users.
#' @export
rcr_fit <- function(endog,
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
#' `print.rcr()` prints the main results of estimating an RCR model.
#'
#' @param x An object of class "`rcr`", usually the result of a call
#'        to `rcr`
#' @param ... Additional arguments to be passed on to `print()`,
#'            for example `digits=`.
#' @returns `print.rcr()` returns its input object `x` invisibly.
#' @seealso [rcr()] to estimate the model
#' @export
print.rcr <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(
    stats::coefficients(x, ...),
    ...
  )
  cat("\n")
  invisible(x)
}

#' Calculate covariance matrix for RCR coefficients
#'
#' `vcov()` returns the variance-covariance matrix of the main parameters
#'  of an rcr model object.
#'
#' @param object An object of class "`rcr`", usually the result of a call
#'        to [rcr()].
#' @param ... Additional arguments (not used).
#' @returns `vcov()` returns the variance-covariance matrix of the main
#'           point-identified parameters of the rcr model object.
#' @seealso [rcr()] to estimate the model, [coef()] to retrieve coefficient
#'          estimates
#' @examplesIf reticulate::py_module_available("rcrbounds")
#' # Estimate the model
#' result <- rcr(weight ~ Time | Diet, ChickWeight)
#' # Use coef() to recover the coefficients
#' coef(result)
#' # Use vcov() to recover the covariance matrix
#' vcov(result)
#' @export
vcov.rcr <- function(object, ...) {
  object$cov.unscaled
}

#' Confidence intervals for RCR results
#'
#' `confint.rcr()` computes confidence intervals for one or more parameters
#' in an RCR model.
#'
#' @param object An object of class "`rcr`", usually the result of a call
#'        to [rcr()].
#' @param parm An optional specification of which parameters are to be given
#'        confidence intervals, either a vector of numbers or a vector
#'        of names. Available parameters include "`rcInf`", "`effectInf`",
#'        "`rc0`", "`effectL`", "`effectH`", and "`effect`".
#'        If missing, all parameters are considered.
#' @param level The confidence level required, on a scale
#'        of 0 to 1. Default is 0.95.
#' @param citype The confidence interval type: "`conservative`" (the default)
#'        "`upper`", "`lower`" or "`Imbens-Manski`". The "`conservative`"
#'        option produces a two-tailed confidence interval for "`effect`"
#'        that covers the *entire* identified set with the specified
#'        asymptotic probability, while the "`Imbens-Manski`" option produces
#'        a two-tailed confidence interval for "`effect`" that covers
#'        *each value* in the identified set with the specified asymptotic
#'        probability. Both options produce a conventional two-tailed
#'        confidence interval for all other parameters. The "`upper`" and
#'        "`lower`" options produce one-tailed confidence intervals for all
#'        parameters.
#' @param ... Additional arguments (not used).
#' @returns `confint.rcr()` returns a matrix with rows corresponding to
#'          the elements of `parm` and columns giving upper and lower
#'          confidence limits for each parameter
#' @seealso [rcr()] to estimate the model, [effect_test()] to perform
#'          hypothesis tests on the causal effect.
#' @examplesIf reticulate::py_module_available("rcrbounds")
#' # Estimate the model
#' result <- rcr(weight ~ Time | Diet, ChickWeight)
#' # Use confint() to produce the confidence intervals
#' confint(result)
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
  if (citype == "conservative" || citype == "Imbens-Manski") {
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
  if (is.null(a)) {
    msg <- paste("invalid citype", citype)
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
#' @param x An object of class "`summary.rcr`", usually the result
#'        of a call to `summary.rcr`.
#' @param ... Additional arguments (e.g., `digits`) to be
#'        passed to `print()`
#' @returns `summary.rcr()` returns an object of class "`rcr.summary`", which
#' is a list containing the following components:
#' \describe{
#'   \item{`call`}{the matched call.}
#'   \item{`terms`}{the terms object used.}
#'   \item{`coefficients`}{a matrix of coefficient estimates, standard errors,
#'         t-statistics and p-values}
#'   \item{`cov.unscaled`}{the estimated covariance matrix for `coefficients`}
#'   \item{`effect_ci`}{the confidence interval for the causal effect}
#'   \item{`citype`}{the confidence interval type used for `effect_ci`}
#'   \item{`level`}{the confidence level used for `effect_ci`}
#' }
#' `print.summary.rcr()` returns its own `x` argument invisibly.
#' @seealso [rcr()] to estimate the model,
#'          [confint.rcr()] to construct confidence intervals for other
#'          parameters,
#'          [effect_test()] to test null hypotheses on the causal effect
#' @examplesIf reticulate::py_module_available("rcrbounds")
#' # Estimate the model
#' result <- rcr(weight ~ Time | Diet, ChickWeight)
#' # Use summary() to produce the summary
#' summary(result)
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
#' The test is based on inverting the Imbens-Manski confidence
#' interval rather than on a specific test statistic or set of critical
#' values.  For example, we reject the null hypothesis at 5% significance
#' if the value under the null is outside of the 95% Imbens-Manski
#' confidence interval.  As a result, this function reports a p-value but
#' there is no test statistic.
#'
#' @param object An object of class "`rcr`", usually the result of a call
#'        to [rcr()].
#' @param h0 The value of the parameter under the null hypothesis.  Default
#'           is zero.
#' @returns `effect_test()` returns the p-value associated with the
#'          specified null hypothesis.
#' @seealso [rcr()] to estimate the model, [confint.rcr()] to construct
#'          confidence intervals.
#' @examplesIf reticulate::py_module_available("rcrbounds")
#' # Estimate the model
#' result <- rcr(weight ~ Time | Diet, ChickWeight)
#' # Use effect_test() with no options to test the null of zero effect
#' effect_test(result)
#' # This tests the null that the effect is equal to 8.0
#' effect_test(result, 8.0)
#' @export
effect_test <- function(object, h0 = 0.0) {
  pmin <- 0.0
  pmid <- 0.5
  pmax <- 1.0
  if (h0 >= object$coefficients[4] && h0 <= object$coefficients[5]) {
    pvalue <- 1.0
  } else {
    while (pmax - pmin > 0.0000001) {
      pmid <- (pmax + pmin) / 2.0
      ci <- effect_ci(object,
        citype = "Imbens-Manski",
        level = pmid
      )
      if (h0 >= ci[1] && h0 <= ci[2]) {
        pmax <- pmid
      } else {
        pmin <- pmid
      }
    }
    pvalue <- 1.0 - pmin
  }
  pvalue
}


#' Plot method for RCR results
#'
#' `plot.rcr()` plots the main results of estimating an RCR model. The plot
#' will display three elements: the `lambda(beta_x)` function, the value of
#' `effectInf`, and the value of `rcInf`.
#'
#' It is recommended that you create an initial plot with the default
#' values of `xlim` and `ylim`, and then adjust these values as needed
#' to zero in on the portion of the graph that is of interest to you.
#'
#' It is also recommended that you use `plot.rcr()` during the same R
#' session in which you created the `rcr` object using the [rcr()]
#' command, and ideally shortly after you create the object. The
#' reason for this is that it uses Python code embedded in the
#' `rcr` object, and this code will typically not be available if the
#' current Python session ends (e.g., by an end in the R session).
#'
#' @param x An object of class "`rcr`", usually the result of a call
#'        to `rcr`
#' @param xlab,ylab Axis labels for x-axis and y-axis, respectively
#' @param xlim,ylim Range of values for x-axis and y-axis, respectively.
#' @param col A vector of colors to be used for the model elements.  `col[1]`
#'            will be the color of the `lambda(beta_x)` function, `col[2]` will
#'            be the color of the vertical line depicting `effectInf`, and
#'            `col[3]`will be the color of the horizontal line depicting
#'            `rcInf`. If fewer than 3 colors are provided, `plot.rcr` will
#'            cycle through the provided colors. Colors can be provided as
#'            simple names, or in any other form described in the "Color
#'            Specification" section of the documentation for [par()].
#' @param lty A vector of line types to be used for the model elements. `lty[1]`
#'            will be the type of the `lambda(beta_x)` function, `lty[2]` will
#'            be the type of the vertical line depicting `effectInf`, and
#'            `lty[3]`will be the type of the horizontal line depicting `rcInf`.
#'            If fewer than 3 line types are provided, `plot.rcr` will cycle
#'            through the provided line types. Line types can be provided as
#'            simple names, or in any other form described in the "Line Type
#'            Specification" section of the documentation for [par()].
#' @param ... Additional graphical parameters to be passed on to
#'            [plot.default()].
#' @returns `plot.rcr()` returns its input object `x` invisibly.
#' @seealso [rcr()] to estimate the model
#' @export
plot.rcr <- function(x,
                     xlab = expression(beta[x]),
                     ylab = expression(lambda),
                     xlim = c(-50, 50),
                     ylim = NULL,
                     col = "black",
                     lty = c("solid", "dotted", "dashed"),
                     ...) {
  # Regenerate Python object if needed
  if (reticulate::py_is_null_xptr(x$pyobj)) {
    x <- eval(x$call)
  }
  # Extend colors and line types
  colors <- rep_len(col, 3)
  ltypes <- rep_len(lty, 3)
  # Get function values
  xvals <- seq(xlim[1], xlim[2], length.out = 100)
  rcvals <- x$pyobj$model$rcvals(xvals, TRUE)
  # Plot function values
  plot(rcvals[[2]], rcvals[[1]],
    type = "l",
    xlab = xlab,
    ylab = ylab,
    col = colors[1],
    lty = ltypes[1],
    ylim = ylim,
    ...
  )
  # Add vertical line at effectInf
  effect_inf <- x$coefficients[2]
  if (effect_inf >= min(xlim) && effect_inf <= max(xlim)) {
    graphics::abline(
      v = effect_inf,
      lty = ltypes[2],
      col = colors[2]
    )
    graphics::mtext(expression(beta[x]^infinity),
      side = 3,
      at = effect_inf,
      col = colors[2]
    )
  }
  # Add horizontal line at rcInf
  rc_inf <- x$coefficients[1]
  plotrange <- graphics::par("usr")
  if (rc_inf >= plotrange[3] && rc_inf <= plotrange[4]) {
    graphics::abline(
      h = rc_inf,
      lty = ltypes[3],
      col = colors[3]
    )
    graphics::mtext(expression(lambda^infinity),
      side = 4,
      at = rc_inf,
      col = colors[3]
    )
  }
  # Return original object invisibly
  invisible(x)
}
