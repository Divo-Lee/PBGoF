#' Unified parametric-bootstrap goodness-of-fit test
#'
#' Perform a parametric-bootstrap goodness-of-fit test for either a
#' skew-normal or a skew-t distribution. The fitted model is re-estimated in
#' every bootstrap sample.
#'
#' @param data A numeric vector with at least 10 finite observations.
#' @param model The fitted model, either `"SN"` or `"ST"`.
#' @param statistic The goodness-of-fit statistic: `"KS"` or `"CvM"`.
#' @param B A positive integer giving the number of bootstrap samples.
#' @param seed A single integer used to initialize the bootstrap random-number
#'   stream, or `NULL` to use the current stream.
#' @param verbose Logical; if `TRUE`, print a short summary.
#'
#' @return A numeric bootstrap p-value with attributes containing the observed
#'   statistic and bootstrap diagnostics. For the skew-t model, the fitted
#'   parameters are also returned in the `fit` attribute.
#'
#' @details
#' This function dispatches to the existing model- and statistic-specific
#' implementations, so the estimation and bootstrap calibration are unchanged.
#' The legacy functions `sn.para.bootstrap.ks.test()`,
#' `sn.para.bootstrap.cvm.test()`, `st.para.bootstrap.ks.test()`, and
#' `st.para.bootstrap.cvm.test()` remain available for compatibility.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' x_sn <- sn::rsn(100, alpha = 3)
#' para.bootstrap.test(x_sn, model = "SN", statistic = "KS", B = 100)
#'
#' x_st <- sn::rst(100, alpha = 3, nu = 8)
#' para.bootstrap.test(x_st, model = "ST", statistic = "CvM", B = 100)
#' }
#'
#' @export
para.bootstrap.test <- function(data = NULL,
                                model = c("SN", "ST"),
                                statistic = c("KS", "CvM"),
                                B = 1000L, seed = 103L,
                                verbose = FALSE) {
  model <- match.arg(toupper(model), c("SN", "ST"))
  statistic <- match.arg(toupper(statistic), c("KS", "CVM"))

  if (identical(model, "SN") && identical(statistic, "KS")) {
    return(sn.para.bootstrap.ks.test(data, B, seed, verbose))
  }
  if (identical(model, "SN") && identical(statistic, "CVM")) {
    return(sn.para.bootstrap.cvm.test(data, B, seed, verbose))
  }
  if (identical(model, "ST") && identical(statistic, "KS")) {
    return(st.para.bootstrap.ks.test(data, B, seed, verbose))
  }
  st.para.bootstrap.cvm.test(data, B, seed, verbose)
}
