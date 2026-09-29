\name{QQ.plot}
\alias{QQ.plot}
\title{QQ Plot for a Fitted Skew-Normal or Skew-t Distribution}
\usage{
QQ.plot(data, model = "SN", main = NULL,
  ylab = "Observed quantiles", pch = 16, cex = 0.65, lwd = 2)
}
\arguments{
\item{data}{A numeric vector containing finite observations.}
\item{model}{Either \code{"SN"} or \code{"ST"}. Other models are not supported.}
\item{main}{Main plot title. The default is \code{NULL}.}
\item{ylab}{Label for the vertical axis.}
\item{pch}{Plotting character for the observed quantiles.}
\item{cex}{Point size.}
\item{lwd}{Line width of the reference line.}
}
\value{
Invisibly returns a list containing the model, fitted parameters, theoretical quantiles, observed quantiles, and reference-line color.
}
\description{
Draws a QQ plot for a fitted skew-normal or skew-t distribution.
}
\details{
This function is intended only for checking skew-normal and skew-t models. The reference line uses color \code{"#D55E00"}. The horizontal-axis label is automatically set to \code{Theoretical SN quantiles} or \code{Theoretical ST quantiles}.
}
\examples{
set.seed(123)
x <- sn::rsn(200, xi = 0, omega = 1, alpha = 5)
QQ.plot(x, model = "SN")
}
