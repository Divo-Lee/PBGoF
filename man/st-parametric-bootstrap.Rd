\name{st-parametric-bootstrap}
\alias{st.para.bootstrap.ks.test}
\alias{st.para.bootstrap.cvm.test}
\title{Parametric Bootstrap Goodness-of-Fit Tests for a Skew-t Distribution}
\usage{
st.para.bootstrap.ks.test(data = NULL, B = 1000L, seed = 103L, verbose = FALSE)

st.para.bootstrap.cvm.test(data = NULL, B = 1000L, seed = 103L, verbose = FALSE)
}
\arguments{
\item{data}{A numeric vector with at least 10 finite observations.}
\item{B}{Positive integer number of bootstrap samples.}
\item{seed}{Non-negative integer or \code{NULL}.}
\item{verbose}{Logical; print a short summary.}
}
\value{A numeric bootstrap p-value with attributes \code{statistic}, \code{B}, \code{valid}, \code{failed}, and \code{fit}.}
\description{Tests a fitted univariate skew-t distribution using KS or Cramer--von Mises statistics. All four skew-t parameters are re-estimated in each bootstrap sample with \code{st.fit.robust()}.}
\details{The KS statistic is \code{sqrt(n) D}. The CvM statistic is \eqn{1/(12n)+\sum_i[U_i-(2i-1)/(2n)]^2}. Failed fits are excluded, and the finite-simulation correction is used.}
\examples{
\donttest{
set.seed(123)
x <- sn::rst(100, alpha=3, nu=8)
st.para.bootstrap.ks.test(x, B=100)
st.para.bootstrap.cvm.test(x, B=100)
}
}
