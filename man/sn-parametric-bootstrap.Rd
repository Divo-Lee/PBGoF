\name{sn-parametric-bootstrap}
\alias{sn.para.bootstrap.ks.test}
\alias{sn.para.bootstrap.cvm.test}
\title{Parametric Bootstrap Goodness-of-Fit Tests for a Skew-Normal Distribution}
\usage{
sn.para.bootstrap.ks.test(data = NULL, B = 1000L, seed = 103L,
  verbose = FALSE)

sn.para.bootstrap.cvm.test(data = NULL, B = 1000L, seed = 103L,
  verbose = FALSE)
}
\arguments{
\item{data}{A numeric vector containing at least 10 finite observations.}
\item{B}{A positive integer giving the number of bootstrap samples.}
\item{seed}{A single integer, or \code{NULL} to use the current random-number
stream. The caller's random-number state is restored on exit.}
\item{verbose}{Logical; whether to print a short summary.}
}
\value{
A single numeric bootstrap p-value with attributes \code{statistic},
\code{B}, \code{valid}, and \code{failed}.
}
\description{
Tests skew-normal goodness of fit using a parametric bootstrap, re-estimating
the model in each bootstrap sample with \code{sn.fit.robust()}.
}
\details{
The KS function uses \eqn{\sqrt{n}D}; the CvM function uses the
Cramér--von Mises statistic. Failed fits are excluded and the p-value uses
the finite-simulation correction \eqn{(1 + \sum I(T_b >= T_0))/(B_{valid}+1)}.
The names ending in \code{.1} are compatibility aliases.
}
\examples{
\donttest{
set.seed(123)
x <- sn::rsn(50, xi = 0, omega = 1, alpha = 3)
sn.para.bootstrap.ks.test(x, B = 99)
sn.para.bootstrap.cvm.test(x, B = 99)
}
}
