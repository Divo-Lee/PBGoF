#' Graphical check of a fitted skew-t distribution
#' @param data Numeric vector.
#' @param breaks Histogram breaks.
#' @param col Histogram fill color.
#' @param border Histogram border color.
#' @param curve_col Fitted density color.
#' @param lwd Fitted density line width.
#' @param main Plot title.
#' @param xlab x-axis label.
#' @param add_rug Logical; add observed values as a rug plot.
#' @param ... Additional arguments to [graphics::hist()].
#' @return Invisibly returns fitted DP parameters, the histogram, and fitted
#'   curve coordinates.
#' @examples
#' set.seed(123)
#' x <- sn::rst(100, alpha=3, nu=8)
#' st.plot.check(x)
#' @export
st.plot.check <- function(data, breaks="FD", col="grey83",
                          border="grey83",
                          curve_col="#6CC6C6", lwd=2, main=NULL,
                          xlab="Values", add_rug=TRUE, ...) {
  .validate_st(data, B=1L, seed=NULL, verbose=FALSE)
  if (length(add_rug)!=1L || !is.logical(add_rug) || is.na(add_rug)) stop("`add_rug` must be TRUE or FALSE.",call.=FALSE)
  if (length(lwd)!=1L || !is.numeric(lwd) || !is.finite(lwd) || lwd<=0) stop("`lwd` must be positive and finite.",call.=FALSE)
  fit <- .st_fit(data)
  if (is.null(fit)) stop("Skew-t fitting failed; the plot cannot be produced.",call.=FALSE)
  histogram <- graphics::hist(data, breaks=breaks, plot=FALSE)
  x_grid <- seq(min(histogram$breaks), max(histogram$breaks), length.out=512L)
  fitted_density <- sn::dst(x_grid, dp=fit)
  y_max <- max(c(histogram$density, fitted_density))
  graphics::plot(histogram, freq=FALSE, col=col, border=border, main=main,
                 xlab=xlab, ylim=c(0,1.08*y_max), ...)
  graphics::lines(x_grid, fitted_density, col=curve_col, lwd=lwd)
  if (add_rug) graphics::rug(data, col=grDevices::adjustcolor(curve_col,alpha.f=.45))
  invisible(list(parameters=fit, histogram=histogram, x=x_grid,
                 density=fitted_density))
}
