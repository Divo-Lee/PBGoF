#' Robust fitting of the skew-normal distribution
#'
#' @description
#' Estimate parameters of the skew-normal (SN) distribution using a
#' numerically robust sequential estimation strategy. Ordinary maximum
#' likelihood estimation
#' (MLE) is first attempted. If numerical instability occurs, maximum
#' penalized likelihood estimation (MPLE) is used, followed by MPLE with
#' modified penalty (\code{MPpenalty}) as a final fallback.
#'
#' This function is designed for parametric bootstrap goodness-of-fit
#' procedures, where stable repeated estimation of skew-normal parameters is
#' required.
#'
#' @param data
#' A numeric vector containing observations.
#'
#' @param para_form
#' Parameterization of the returned estimates.
#'
#' \code{"DP"} returns direct parameters:
#' \code{xi}, \code{omega}, and \code{alpha}.
#'
#' \code{"CP"} returns centered parameters:
#' \code{mean}, \code{sd}, and \code{gamma1}.
#'
#' Default is \code{"DP"}.
#'
#' @return
#' A named numeric vector containing parameter estimates and their standard
#' errors.
#'
#' For \code{para_form="DP"}:
#'
#' \itemize{
#'   \item \code{xi}: location parameter estimate.
#'   \item \code{se.xi}: standard error of \code{xi}.
#'   \item \code{omega}: scale parameter estimate.
#'   \item \code{se.omega}: standard error of \code{omega}.
#'   \item \code{alpha}: shape parameter estimate.
#'   \item \code{se.alpha}: standard error of \code{alpha}.
#' }
#'
#' For \code{para_form="CP"}:
#'
#' \itemize{
#'   \item \code{mean}: mean parameter estimate.
#'   \item \code{se.mean}: standard error of \code{mean}.
#'   \item \code{sd}: standard deviation estimate.
#'   \item \code{se.sd}: standard error of \code{sd}.
#'   \item \code{gamma1}: skewness parameter estimate.
#'   \item \code{se.gamma1}: standard error of \code{gamma1}.
#' }
#'
#' If all estimation procedures fail, a vector of \code{NA} values is returned.
#'
#' @details
#' The estimation is performed using the \code{\link[sn]{selm}} function from
#' the \pkg{sn} package.
#'
#' The estimation procedure follows:
#'
#' \enumerate{
#'   \item ordinary maximum likelihood estimation (MLE);
#'   \item maximum penalized likelihood estimation (MPLE);
#'   \item MPLE with modified penalty (\code{MPpenalty}).
#' }
#'
#' Here, "robust" refers to robustness against numerical fitting failures; it
#' does not mean robustness against outliers or model contamination.
#'
#' This strategy improves numerical stability for highly skewed samples or
#' small sample sizes where the shape parameter may become difficult to
#' estimate.
#'
#' @references
#'
#' Azzalini, A. (1985).
#' A class of distributions which includes the normal ones.
#' \emph{Scandinavian Journal of Statistics}, 12(2), 171--178.
#'
#' Azzalini, A. and Capitanio, A. (2014).
#' \emph{The Skew-Normal and Related Families}.
#' Cambridge University Press.
#'
#' Azzalini, A. and Arellano-Valle, R. B. (2013).
#' Maximum penalized likelihood estimation for skew-normal and skew-t
#' distributions.
#' \emph{Journal of Statistical Planning and Inference}, 143, 419--433.
#'
#' @importFrom sn selm
#'
#' @examples
#'
#' set.seed(123)
#'
#' x <- sn::rsn(
#'   100,
#'   xi = 0,
#'   omega = 1,
#'   alpha = 5
#' )
#'
#' sn.fit.robust(x, para_form="DP")
#'
#' sn.fit.robust(x, para_form="CP")
#'
#' @export
sn.fit.robust <- function(data = NULL,
                          para_form = c("DP", "CP")) {
  para_form <- match.arg(para_form)

  # Data checking ---------------------------------------------------------
  if (is.null(data)) {
    stop("`data` is missing.", call. = FALSE)
  }
  if (!typeof(data) %in% c("integer", "double") || !is.null(dim(data))) {
    stop("`data` must be a numeric vector.", call. = FALSE)
  }
  if (length(data) < 10L) {
    stop(
      "`data` must contain at least 10 observations.",
      call. = FALSE
    )
  }
  if (any(!is.finite(data))) {
    stop("`data` must contain only finite values (no NA, NaN, or Inf).",
         call. = FALSE)
  }
  if (diff(range(data)) == 0) {
    stop("`data` must contain at least two distinct values.", call. = FALSE)
  }

  # Extract and validate the requested parameterization ------------------
  extract_par <- function(fit) {
    parameter_type <- if (para_form == "DP") "dp" else "cp"
    parameter_sets <- tryCatch(
      methods::slot(fit, "param"),
      error = function(e) NULL
    )
    variance_sets <- tryCatch(
      methods::slot(fit, "param.var"),
      error = function(e) NULL
    )

    estimates <- parameter_sets[[parameter_type]]
    variance <- variance_sets[[parameter_type]]
    if (is.null(estimates) || length(estimates) < 3L ||
        is.null(variance) || NROW(variance) < 3L || NCOL(variance) < 3L) {
      return(NULL)
    }

    est <- unname(estimates[seq_len(3L)])
    variance_diagonal <- diag(variance)[seq_len(3L)]
    # Tiny negative values can arise solely from floating-point roundoff.
    tolerance <- 100 * .Machine$double.eps *
      max(1, max(abs(variance_diagonal), na.rm = TRUE))
    if (any(!is.finite(variance_diagonal)) ||
        any(variance_diagonal < -tolerance)) {
      return(NULL)
    }
    variance_diagonal[
      variance_diagonal < 0 & variance_diagonal >= -tolerance
    ] <- 0
    se <- unname(sqrt(variance_diagonal))
    if (!all(is.finite(c(est, se))) || any(se < 0)) {
      return(NULL)
    }

    # The second component is a scale in both parameterizations.
    if (est[2L] <= 0) {
      return(NULL)
    }

    if (para_form == "DP") {
      return(c(
        xi = est[1L], se.xi = se[1L],
        omega = est[2L], se.omega = se[2L],
        alpha = est[3L], se.alpha = se[3L]
      ))
    }

    c(
      mean = est[1L], se.mean = se[1L],
      sd = est[2L], se.sd = se[2L],
      gamma1 = est[3L], se.gamma1 = se[3L]
    )
  }

  attempt <- function(...) {
    fit <- tryCatch(
      sn::selm(data ~ 1, family = "SN", ...),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    extract_par(fit)
  }

  # MLE, default Q-penalized MPLE, then matching-prior MPLE ---------------
  result <- attempt()
  if (!is.null(result)) {
    return(result)
  }

  result <- attempt(method = "MPLE")
  if (!is.null(result)) {
    return(result)
  }

  result <- attempt(method = "MPLE", penalty = "MPpenalty")
  if (!is.null(result)) {
    return(result)
  }

  warning("All skew-normal estimation procedures failed.", call. = FALSE)

  if (para_form == "DP") {
    return(c(
      xi = NA_real_, se.xi = NA_real_,
      omega = NA_real_, se.omega = NA_real_,
      alpha = NA_real_, se.alpha = NA_real_
    ))
  }

  c(
    mean = NA_real_, se.mean = NA_real_,
    sd = NA_real_, se.sd = NA_real_,
    gamma1 = NA_real_, se.gamma1 = NA_real_
  )
}
