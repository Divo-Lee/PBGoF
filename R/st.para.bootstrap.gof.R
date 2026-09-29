.validate_st <- function(data, B, seed, verbose) {
  if (!is.numeric(data) || !is.null(dim(data)) || length(data) < 10L || any(!is.finite(data)) || diff(range(data)) == 0) stop("`data` must be a non-constant numeric vector with at least 10 finite observations.", call.=FALSE)
  if (length(B)!=1L || !is.numeric(B) || !is.finite(B) || B<1 || B!=floor(B)) stop("`B` must be a positive integer.",call.=FALSE)
  if (!is.null(seed) && (length(seed)!=1L || !is.numeric(seed) || !is.finite(seed) || seed<0 || seed!=floor(seed))) stop("`seed` must be NULL or a non-negative integer.",call.=FALSE)
  if (length(verbose)!=1L || !is.logical(verbose) || is.na(verbose)) stop("`verbose` must be TRUE or FALSE.",call.=FALSE)
  list(data=as.double(data),B=as.integer(B),seed=if(is.null(seed)) NULL else as.integer(seed),verbose=verbose)
}
.st_fit <- function(x) { f<-suppressWarnings(tryCatch(st.fit.robust(x, para_form="DP"),error=function(e) NULL)); if(is.null(f)||any(!is.finite(f[c("xi","omega","alpha","nu")]))||f[["omega"]]<=0||f[["nu"]]<=0) NULL else f[c("xi","omega","alpha","nu")] }
.st_ks <- function(x,d) {u<-sort(sn::pst(x,dp=d)); n<-length(x); max(seq_len(n)/n-u,u-(seq_len(n)-1)/n)*sqrt(n)}
.st_cvm <- function(x,d) {u<-sort(sn::pst(x,dp=d)); n<-length(x); 1/(12*n)+sum((u-(2*seq_len(n)-1)/(2*n))^2)}
.st_boot <- function(data, B, seed, verbose, statistic, name) {
  a <- .validate_st(data, B, seed, verbose)
  fit <- .st_fit(a$data)
  if (is.null(fit)) stop("Skew-t fitting failed for observed data.", call.=FALSE)
  obs <- statistic(a$data, fit)
  if (!is.finite(obs)) stop("Observed goodness-of-fit statistic is not finite.", call.=FALSE)
  if (!is.null(a$seed)) {
    existed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    if (existed) old <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    on.exit({
      if (existed) assign(".Random.seed", old, envir=.GlobalEnv)
      else if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv)
    }, add=TRUE)
    set.seed(a$seed)
  }
  z <- rep.int(NA_real_, a$B)
  for (i in seq_len(a$B)) {
    sample <- sn::rst(length(a$data), dp=fit)
    sample_fit <- .st_fit(sample)
    if (!is.null(sample_fit)) z[i] <- statistic(sample, sample_fit)
  }
  z <- z[is.finite(z)]
  valid <- length(z)
  if (!valid) stop("All bootstrap fits failed.", call.=FALSE)
  if (valid < a$B / 2) warning(sprintf("Only %.1f%% of bootstrap fits were valid.", 100*valid/a$B), call.=FALSE)
  p <- (1 + sum(z >= obs)) / (valid + 1)
  attr(p,"statistic") <- obs; attr(p,"B") <- a$B
  attr(p,"valid") <- valid; attr(p,"failed") <- a$B-valid
  attr(p,"fit") <- fit
  if (a$verbose) cat(sprintf("Observed %s statistic = %.4f\nValid bootstraps = %d/%d\np-value = %.6g\n", name, obs, valid, a$B, p))
  p
}

#' Parametric bootstrap goodness-of-fit tests for a skew-t distribution
#'
#' Tests a fitted univariate skew-t distribution, re-estimating all four
#' parameters in every bootstrap sample with [st.fit.robust()].
#' @param data A numeric vector with at least 10 finite observations.
#' @param B Number of bootstrap samples.
#' @param seed Non-negative integer or `NULL`.
#' @param verbose Logical; print a short summary.
#' @return A numeric bootstrap p-value with attributes `statistic`, `B`,
#'   `valid`, `failed`, and `fit`.
#' @details The KS test uses `sqrt(n) D`; the CvM test uses
#'   `1/(12n) + sum((U[i] - (2i-1)/(2n))^2)`. The finite-simulation correction
#'   `(1 + sum(T.boot >= T.obs))/(1 + B.valid)` is used.
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- sn::rst(100, alpha=3, nu=8)
#' st.para.bootstrap.ks.test(x, B=100)
#' st.para.bootstrap.cvm.test(x, B=100)
#' }
#' @name st-parametric-bootstrap
NULL

#' @rdname st-parametric-bootstrap
#' @export
st.para.bootstrap.ks.test <- function(data=NULL,B=1000L,seed=103L,verbose=FALSE) .st_boot(data,B,seed,verbose,.st_ks,"KS")
#' @rdname st-parametric-bootstrap
#' @export
st.para.bootstrap.cvm.test <- function(data=NULL,B=1000L,seed=103L,verbose=FALSE) .st_boot(data,B,seed,verbose,.st_cvm,"CvM")
