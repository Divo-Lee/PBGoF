#' Quantile-quantile plot for a fitted skew-normal or skew-t distribution
#'
#' @description
#' Draws a QQ plot comparing the observed sample quantiles with theoretical
#' quantiles from a fitted skew-normal or skew-t distribution.
#'
#' @param data A numeric vector containing finite observations.
#' @param model Character string specifying `"SN"` or `"ST"`.
#' @param main Main plot title. The default is `NULL`.
#' @param ylab Label for the vertical axis.
#' @param pch Plotting character for the observed quantiles.
#' @param cex Point size.
#' @param lwd Line width of the reference line.
#'
#' @return Invisibly returns a list containing the model, fitted parameters,
#'   theoretical quantiles, observed quantiles, and reference-line color.
#'
#' @details
#' This function is intended only for checking skew-normal and skew-t models.
#' The reference line is drawn with color `"#D55E00"`.
#'
#' @examples
#' set.seed(123)
#' x <- sn::rsn(200, xi = 0, omega = 1, alpha = 5)
#' QQ.plot(x, model = "SN")
#'
#' @export
QQ.plot <- function(data,
                    model = "SN",
                    main = NULL,
                    ylab = "Observed quantiles",
                    pch = 16,
                    cex = 0.65,
                    lwd = 2) {
  if (!is.numeric(data) || is.matrix(data) || is.array(data) ||
      length(data) < 2L || any(!is.finite(data))) {
    stop("data must be a numeric vector containing at least two finite observations.",
         call. = FALSE)
  }
  if (length(model) != 1L || is.na(model) ||
      !toupper(model) %in% c("SN", "ST")) {
    stop("QQ.plot() only supports model = 'SN' or model = 'ST'.",
         call. = FALSE)
  }
  model <- toupper(model)
  if (length(cex) != 1L || !is.numeric(cex) || !is.finite(cex) || cex <= 0) {
    stop("cex must be a positive finite number.", call. = FALSE)
  }
  if (length(lwd) != 1L || !is.numeric(lwd) || !is.finite(lwd) || lwd <= 0) {
    stop("lwd must be a positive finite number.", call. = FALSE)
  }

  if (model == "SN") {
    fit <- sn.fit.robust(data, para_form = "DP")
    parameters <- fit[c("xi", "omega", "alpha")]
    if (any(!is.finite(parameters))) {
      stop("Skew-normal fitting failed; the QQ plot cannot be produced.",
           call. = FALSE)
    }
    theoretical <- sn::qsn(
      stats::ppoints(length(data)),
      xi = parameters[["xi"]],
      omega = parameters[["omega"]],
      alpha = parameters[["alpha"]]
    )
    xlab <- "Theoretical SN quantiles"
  } else {
    .validate_st(data, B = 1L, seed = NULL, verbose = FALSE)
    parameters <- .st_fit(data)
    if (is.null(parameters) || any(!is.finite(parameters))) {
      stop("Skew-t fitting failed; the QQ plot cannot be produced.",
           call. = FALSE)
    }
    theoretical <- sn::qst(stats::ppoints(length(data)), dp = parameters)
    xlab <- "Theoretical ST quantiles"
  }

  observed <- sort(data)
  limits <- range(c(theoretical, observed), finite = TRUE)
  reference_col <- "#D55E00"
  graphics::plot(
    theoretical,
    observed,
    pch = pch,
    cex = cex,
    xlim = limits,
    ylim = limits,
    main = main,
    xlab = xlab,
    ylab = ylab
  )
  graphics::abline(0, 1, col = reference_col, lwd = lwd, lty = 2)

  invisible(list(
    model = model,
    parameters = parameters,
    theoretical = theoretical,
    observed = observed,
    reference_col = reference_col
  ))
}
