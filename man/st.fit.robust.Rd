\name{st.fit.robust}
\alias{st.fit.robust}
\title{Numerically Robust Fitting of a Univariate Skew-t Distribution}
\usage{st.fit.robust(data = NULL, para_form = c("DP", "CP"))}
\arguments{
\item{data}{A numeric vector with at least 10 finite observations.}
\item{para_form}{Parameterization to return: \code{"DP"} or \code{"CP"}.}
}
\value{A named numeric vector containing four estimates and their standard errors.}
\description{Fits a univariate skew-t distribution using ordinary maximum likelihood, followed by Q-penalized maximum likelihood with alternative optimizers when needed.}
\details{DP contains \code{xi}, \code{omega}, \code{alpha}, and \code{nu}. CP contains the mean, standard deviation, skewness, and kurtosis when these moments exist. Robust refers to numerical fallback, not resistance to outliers.}
\examples{
set.seed(123)
x <- sn::rst(100, alpha=3, nu=8)
st.fit.robust(x, "DP")
}
