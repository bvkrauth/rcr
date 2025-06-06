% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcr.R
\name{effect_test}
\alias{effect_test}
\title{Hypothesis testing for RCR causal effect}
\usage{
effect_test(object, h0 = 0)
}
\arguments{
\item{object}{An object of class "\code{rcr}", usually the result of a call
to \code{\link[=rcr]{rcr()}}.}

\item{h0}{The value of the parameter under the null hypothesis.  Default
is zero.}
}
\value{
\code{effect_test()} returns the p-value associated with the
specified null hypothesis.
}
\description{
\code{effect_test()} allows a user to test a point null hypothesis on the
(interval identified) causal effect parameter in the RCR model.
}
\details{
The test is based on inverting the Imbens-Manski confidence
interval rather than on a specific test statistic or set of critical
values.  For example, we reject the null hypothesis at 5\% significance
if the value under the null is outside of the 95\% Imbens-Manski
confidence interval.  As a result, this function reports a p-value but
there is no test statistic.
}
\examples{
\dontshow{if (reticulate::py_module_available("rcrbounds")) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Estimate the model
result <- rcr(weight ~ Time | Diet, ChickWeight)
# Use effect_test() with no options to test the null of zero effect
effect_test(result)
# This tests the null that the effect is equal to 8.0
effect_test(result, 8.0)
\dontshow{\}) # examplesIf}
}
\seealso{
\code{\link[=rcr]{rcr()}} to estimate the model, \code{\link[=confint.rcr]{confint.rcr()}} to construct
confidence intervals.
}
