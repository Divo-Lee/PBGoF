#' Graphical check of a fitted skew-normal distribution
#'
#' @description
#' Draw a density-scale histogram of the observed data and overlay the density
#' of a skew-normal distribution fitted by [sn.fit.robust()].
#'
#' @param data A numeric vector containing at least 10 finite observations.
#' @param breaks Histogram breaks passed to [graphics::hist()].
#' @param col Fill color for the histogram.
#' @param border Border color for histogram bars.
#' @param curve_col Color of the fitted skew-normal density curve.
#' @param lwd Line width of the fitted density curve.
#' @param main Main plot title.
#' @param xlab Label for the horizontal axis.
#' @param add_rug Logical; if TRUE, add observed values as a rug plot.
#' @param ... Additional graphical arguments passed to [graphics::hist()].
#'
#' @return Invisibly returns a list with the fitted DP parameters, histogram
#'   object, and the coordinates of the fitted curve.
#'
#' @details
#' This function is intended as a visual model check. It complements, but does
#' not replace, the formal goodness-of-fit tests provided by PBGoF.
#'
#' @examples
#' set.seed(123)
#' x <- sn::rsn(200, xi = 0, omega = 1, alpha = 5)
#' sn.plot.check(x)
#'
#' @export
sn.plot.check <- function(data,
                          breaks = "FD",
                          col = "grey85",
                          border = "white",
                          curve_col = "#D55E00",
                          lwd = 2,
                          main = "Skew-normal fit check",
                          xlab = deparse(substitute(data)),
                          add_rug = TRUE,
                          ...) {
  fit <- sn.fit.robust(data, para_form = "DP")
  parameters <- fit[c("xi", "omega", "alpha")]
  if (any(!is.finite(parameters))) {
    stop("Skew-normal fitting failed; the plot cannot be produced.",
         call. = FALSE)
  }

  if (length(add_rug) != 1L || !is.logical(add_rug) || is.na(add_rug)) {
    stop("add_rug must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(lwd) != 1L || !is.numeric(lwd) || !is.finite(lwd) || lwd <= 0) {
    stop("lwd must be a positive finite number.", call. = FALSE)
  }

  histogram <- graphics::hist(
    data,
    breaks = breaks,
    plot = FALSE
  )

  x_grid <- seq(
    min(histogram$breaks),
    max(histogram$breaks),
    length.out = 512L
  )
  fitted_density <- sn::dsn(
    x_grid,
    xi = parameters[["xi"]],
    omega = parameters[["omega"]],
    alpha = parameters[["alpha"]]
  )

  y_max <- max(c(histogram$density, fitted_density))
  graphics::plot(
    histogram,
    freq = FALSE,
    col = col,
    border = border,
    main = main,
    xlab = xlab,
    ylim = c(0, 1.08 * y_max),
    ...
  )
  graphics::lines(
    x_grid,
    fitted_density,
    col = curve_col,
    lwd = lwd
  )
  if (add_rug) {
    graphics::rug(data, col = grDevices::adjustcolor(curve_col, alpha.f = 0.45))
  }
  graphics::legend(
    "topright",
    legend = "Fitted skew-normal density",
    col = curve_col,
    lwd = lwd,
    bty = "n"
  )

  invisible(list(
    parameters = parameters,
    histogram = histogram,
    x = x_grid,
    density = fitted_density
  ))
}
