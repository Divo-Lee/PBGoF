\name{sn.fit.robust}
\alias{sn.fit.robust}
\title{Numerically Robust Fitting of the Skew-Normal Distribution}
\usage{
sn.fit.robust(data = NULL, para_form = c("DP", "CP"))
}
\arguments{
\item{data}{A numeric vector containing at least 10 finite observations.}
\item{para_form}{Parameterization of the result: \code{"DP"} for direct
parameters or \code{"CP"} for centered parameters.}
}
\value{
A named numeric vector containing parameter estimates and standard errors.
If all fitting procedures fail, all entries are \code{NA}.
}
\description{
Fits a skew-normal distribution by MLE, then falls back to default penalized
MLE and matching-prior penalized MLE when necessary.
}
\details{
Robustness here concerns numerical fitting failures, not resistance to
outliers or model contamination.
}
\examples{
\donttest{
set.seed(123)
x <- sn::rsn(100, xi = 0, omega = 1, alpha = 5)
sn.fit.robust(x, "DP")
sn.fit.robust(x, "CP")
}
}
