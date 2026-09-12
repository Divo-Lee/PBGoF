\name{sn.plot.check}
\alias{sn.plot.check}
\title{Graphical Check of a Fitted Skew-Normal Distribution}
\usage{
sn.plot.check(
  data,
  breaks = "FD",
  col = "grey85",
  border = "white",
  curve_col = "#D55E00",
  lwd = 2,
  main = "Skew-normal fit check",
  xlab = deparse(substitute(data)),
  add_rug = TRUE,
  ...
)
}
\arguments{
\item{data}{A numeric vector containing at least 10 finite observations.}
\item{breaks}{Histogram breaks passed to \code{graphics::hist()}.}
\item{col}{Fill color for the histogram.}
\item{border}{Border color for histogram bars.}
\item{curve_col}{Color of the fitted skew-normal density curve.}
\item{lwd}{Line width of the fitted density curve.}
\item{main}{Main plot title.}
\item{xlab}{Label for the horizontal axis.}
\item{add_rug}{Logical; whether to add observations as a rug plot.}
\item{...}{Additional graphical arguments passed to \code{graphics::hist()}.}
}
\value{
Invisibly returns a list containing the fitted DP parameters, histogram object,
and coordinates of the fitted density curve.
}
\description{
Draws a density-scale histogram of the observed data and overlays a
skew-normal density fitted by \code{sn.fit.robust()}.
}
\details{
This visual model check complements, but does not replace, the formal
goodness-of-fit tests provided by PBGoF.
}
\examples{
set.seed(123)
x <- sn::rsn(100, xi = 0, omega = 1, alpha = 4)
sn.plot.check(x)
}
