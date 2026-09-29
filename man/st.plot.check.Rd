\name{st.plot.check}
\alias{st.plot.check}
\title{Graphical Check of a Fitted Skew-t Distribution}
\usage{
st.plot.check(data, breaks = "FD", col = "grey83", border = "grey83",
  curve_col = "#6CC6C6", lwd = 2, main = NULL, xlab = "Values",
  add_rug = TRUE, ...)
}
\arguments{
\item{data}{A numeric vector with at least 10 finite observations.}
\item{breaks}{Histogram breaks passed to \code{hist}.}
\item{col}{Histogram fill color.}
\item{border}{Histogram border color.}
\item{curve_col}{Fitted density color. The default is \code{"#6CC6C6"}.}
\item{lwd}{Fitted density line width.}
\item{main}{Main title.}
\item{xlab}{Horizontal-axis label.}
\item{add_rug}{Logical; add observed values as a rug plot.}
\item{...}{Additional graphical arguments passed to \code{hist}.}
}
\value{
Invisibly returns fitted DP parameters, histogram information, and fitted
curve coordinates.
}
\description{
Draws a density-scale histogram and overlays the skew-t density fitted by
\code{st.fit.robust()}.
}
\examples{
set.seed(123)
x <- sn::rst(100, alpha=3, nu=8)
st.plot.check(x)
}
