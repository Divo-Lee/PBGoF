\name{cdf.plot.check}
\alias{cdf.plot.check}
\title{Graphical Check Using Empirical and Fitted Distribution Functions}
\usage{
cdf.plot.check(data, model = "SN", curve_col = NULL, lwd = 2,
  main = NULL, xlab = "Values", ylab = "Distribution function",
  legend = FALSE, ...)
}
\arguments{
\item{data}{A numeric vector containing finite observations.}
\item{model}{Either \code{"SN"} or \code{"ST"}.}
\item{curve_col}{Color of the fitted CDF. The default is \code{"#6CC6C6"} for both SN and ST.}
\item{lwd}{Line width of the fitted CDF. The ECDF uses \code{lwd - 0.5}.}
\item{main}{Main plot title. The default is \code{NULL}.}
\item{xlab}{Label for the horizontal axis.}
\item{ylab}{Label for the vertical axis.}
\item{legend}{Logical; if \code{TRUE}, add a legend for the ECDF and fitted CDF.}
\item{...}{Additional graphical arguments passed to \code{graphics::plot()}.}
}
\value{
Invisibly returns a list containing the model, fitted parameters, grid values, fitted CDF values, and ECDF object.
}
\description{
Plots the empirical cumulative distribution function together with the CDF of a fitted skew-normal or skew-t distribution.
}
\details{
The fitted CDF is overlaid on the ECDF with the standard horizontal reference lines at \code{y = 0} and \code{y = 1}.
}
\examples{
set.seed(123)
x <- sn::rsn(200, xi = 0, omega = 1, alpha = 5)
cdf.plot.check(x, model = "SN")
}
