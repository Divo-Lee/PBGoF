#' Robust fitting of a univariate skew-t distribution
#'
#' Estimate skew-t parameters using a sequential MLE/MPLE fallback strategy.
#' The returned direct parameterization is `(xi, omega, alpha, nu)`; the
#' centered parameterization is `(mean, sd, gamma1, gamma2)`.
#' @param data A numeric vector with at least 10 finite observations.
#' @param para_form Either `"DP"` or `"CP"`.
#' @return A named numeric vector of estimates and standard errors. For DP,
#'   the components are `xi`, `omega`, `alpha`, `nu` and their standard errors.
#'   CP additionally requires the corresponding skew-t moments to exist.
#' @details Ordinary MLE is attempted first. If it is numerically invalid,
#'   Q-penalized MPLE is attempted with `nlminb` and then `BFGS`. Here robust
#'   means robustness to numerical fitting failure, not robustness to outliers.
#' @examples
#' set.seed(123)
#' x <- sn::rst(100, xi = 0, omega = 1, alpha = 3, nu = 8)
#' st.fit.robust(x, "DP")
#' st.fit.robust(x, "CP")
#' @export
st.fit.robust <- function(data = NULL, para_form = c("DP", "CP")) {
  para_form <- match.arg(para_form)
  if (is.null(data) || !is.numeric(data) || !is.null(dim(data)) || length(data) < 10L || any(!is.finite(data)) || diff(range(data)) == 0) stop("`data` must be a non-constant numeric vector with at least 10 finite observations.", call.=FALSE)
  extract <- function(fit) {
    p <- tryCatch(methods::slot(fit,"param"), error=function(e) NULL)
    v <- tryCatch(methods::slot(fit,"param.var"), error=function(e) NULL)
    if (is.null(p)||is.null(v)||is.null(p$dp)||is.null(v$dp)||length(p$dp)<4L||any(dim(v$dp)<4L)) return(NULL)
    est <- unname(p$dp[seq_len(4L)]); d <- unname(diag(v$dp)[seq_len(4L)])
    tol <- 100*.Machine$double.eps*max(1,max(abs(d),na.rm=TRUE)); if(any(!is.finite(d))||any(d < -tol)) return(NULL); d[d<0] <- 0; se <- unname(sqrt(d))
    if(!all(is.finite(c(est,se)))||est[2L]<=0||est[4L]<=0) return(NULL)
    if(para_form=="DP") return(c(xi=est[1],se.xi=se[1],omega=est[2],se.omega=se[2],alpha=est[3],se.alpha=se[3],nu=est[4],se.nu=se[4]))
    cp <- p$cp; vc <- v$cp; if(is.null(cp)||is.null(vc)||length(cp)<4L||any(dim(vc)<4L)) return(NULL); ce <- unname(cp[seq_len(4L)]); cs <- sqrt(pmax(0,diag(vc)[seq_len(4L)])); if(!all(is.finite(c(ce,cs)))) return(NULL); c(mean=ce[1],se.mean=cs[1],sd=ce[2],se.sd=cs[2],gamma1=ce[3],se.gamma1=cs[3],gamma2=ce[4],se.gamma2=cs[4])
  }
  attempt <- function(...) {
    fit <- tryCatch(sn::selm(data ~ 1, family = "ST", ...),
                    error = function(e) NULL)
    if (is.null(fit)) NULL else extract(fit)
  }
  attempts <- list(
    list(opt.method = "nlminb"),
    list(method = "MPLE", penalty = "Qpenalty", opt.method = "nlminb"),
    list(method = "MPLE", penalty = "Qpenalty", opt.method = "BFGS")
  )
  for (args in attempts) {
    ans <- suppressWarnings(do.call(attempt, args))
    if (!is.null(ans)) return(ans)
  }
  warning("All skew-t estimation procedures failed.",call.=FALSE)
  if(para_form=="DP") return(stats::setNames(rep(NA_real_,8),c("xi","se.xi","omega","se.omega","alpha","se.alpha","nu","se.nu")))
  stats::setNames(rep(NA_real_,8),c("mean","se.mean","sd","se.sd","gamma1","se.gamma1","gamma2","se.gamma2"))
}
