#' Graphical check using the empirical and fitted distribution functions
#'
#' @description
#' Plots the empirical cumulative distribution function (ECDF) together with
#' the CDF of a fitted skew-normal or skew-t distribution.
#'
#' @param data A numeric vector containing finite observations.
#' @param model Character string specifying the fitted model: `"SN"` or
#'   `"ST"`. The default is `"SN"`.
#' @param curve_col Color of the fitted CDF. By default, `"#6CC6C6"` is used
#'   for both SN and ST, matching [sn.plot.check()] and [st.plot.check()].
#' @param lwd Line width of the fitted CDF. The ECDF is drawn with line width
#'   `lwd - 0.5`.
#' @param main Main plot title. The default is `NULL`.
#' @param xlab Label for the horizontal axis.
#' @param ylab Label for the vertical axis.
#' @param legend Logical; if `TRUE`, add a legend for the ECDF and fitted CDF.
#' @param ... Additional graphical arguments passed to [graphics::plot()].
#'
#' @return Invisibly returns a list containing the model, fitted parameters,
#'   grid values, fitted CDF values, and ECDF object.
#'
#' @details
#' The fitted CDF is overlaid on the ECDF with the standard horizontal
#' reference lines at $y=0$ and $y=1$.
#'
#' @examples
#' set.seed(123)
#' x <- sn::rsn(200, xi = 0, omega = 1, alpha = 5)
#' cdf.plot.check(x, model = "SN")
#'
#' @export
cdf.plot.check <- function(data,
                           model = "SN",
                           curve_col = NULL,
                           lwd = 2,
                           main = NULL,
                           xlab = "Values",
                           ylab = "Distribution function",
                           legend = FALSE,
                           ...) {
  if (!is.numeric(data) || is.matrix(data) || is.array(data) ||
      length(data) < 2L || any(!is.finite(data))) {
    stop("data must be a numeric vector containing at least two finite observations.",
         call. = FALSE)
  }
  if (length(model) != 1L || is.na(model) ||
      !toupper(model) %in% c("SN", "ST")) {
    stop("model must be either 'SN' or 'ST'.", call. = FALSE)
  }
  model <- toupper(model)
  if (length(lwd) != 1L || !is.numeric(lwd) || !is.finite(lwd) || lwd <= 0.5) {
    stop("lwd must be greater than 0.5.", call. = FALSE)
  }
  if (length(legend) != 1L || !is.logical(legend) || is.na(legend)) {
    stop("legend must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(curve_col)) {
    curve_col <- "#6CC6C6"
  }
  if (length(curve_col) != 1L || is.na(curve_col)) {
    stop("curve_col must be a single non-missing color specification.",
         call. = FALSE)
  }

  if (model == "SN") {
    fit <- sn.fit.robust(data, para_form = "DP")
    parameters <- fit[c("xi", "omega", "alpha")]
    if (any(!is.finite(parameters))) {
      stop("Skew-normal fitting failed; the plot cannot be produced.",
           call. = FALSE)
    }
    fitted_cdf <- function(x) {
      sn::psn(x, xi = parameters[["xi"]],
              omega = parameters[["omega"]],
              alpha = parameters[["alpha"]])
    }
  } else {
    .validate_st(data, B = 1L, seed = NULL, verbose = FALSE)
    parameters <- .st_fit(data)
    if (is.null(parameters) || any(!is.finite(parameters))) {
      stop("Skew-t fitting failed; the plot cannot be produced.",
           call. = FALSE)
    }
    fitted_cdf <- function(x) sn::pst(x, dp = parameters)
  }

  x_grid <- seq(min(data), max(data), length.out = 512L)
  fitted_values <- fitted_cdf(x_grid)
  empirical <- stats::ecdf(data)

  graphics::plot(
    empirical,
    verticals = TRUE,
    do.points = FALSE,
    col.01line = "gray70",
    lwd = lwd - 0.5,
    main = main,
    xlab = xlab,
    ylab = ylab,
    ...
  )
  graphics::lines(x_grid, fitted_values, col = curve_col, lwd = lwd)
  if (legend) {
    graphics::legend(
      "topleft",
      legend = c("ECDF", paste0("Fitted ", model, " CDF")),
      col = c("black", curve_col),
      lwd = c(lwd - 0.5, lwd),
      bty = "n"
    )
  }

  invisible(list(
    model = model,
    parameters = parameters,
    x = x_grid,
    cdf = fitted_values,
    ecdf = empirical,
    curve_col = curve_col,
    fitted_lwd = lwd,
    ecdf_lwd = lwd - 0.5,
    legend = legend
  ))
}
