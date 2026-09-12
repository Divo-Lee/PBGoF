#' Parametric bootstrap goodness-of-fit tests for a skew-normal distribution
#'
#' @description
#' Perform a parametric bootstrap goodness-of-fit test for a univariate
#' skew-normal distribution. Model parameters are re-estimated in every
#' bootstrap sample with [sn.fit.robust()].
#'
#' `sn.para.bootstrap.ks.test()` uses the scaled Kolmogorov--Smirnov statistic
#' `sqrt(n) * D`, while `sn.para.bootstrap.cvm.test()` uses the
#' Cramer--von Mises statistic.
#'
#' @param data A numeric vector with at least 10 finite observations.
#' @param B A positive integer giving the number of bootstrap samples.
#' @param seed A single integer used to initialize the bootstrap random-number
#'   stream, or `NULL` to use the current stream. The caller's random-number
#'   state is restored on exit.
#' @param verbose Logical; if `TRUE`, print a short summary.
#'
#' @return A single numeric bootstrap p-value. The value has attributes
#'   `statistic`, `B`, `valid`, and `failed`, containing the observed statistic
#'   and bootstrap diagnostics.
#'
#' @details
#' Failed bootstrap fits are excluded from the p-value calculation. A warning
#' is produced when fewer than half of the requested bootstrap samples are
#' valid. The finite-simulation correction
#' `(1 + sum(T_boot >= T_obs)) / (1 + B_valid)` is used.
#'
#' @references
#' Babu, G. J. and Rao, C. R. (2004). Goodness-of-fit tests when parameters
#' are estimated. *Sankhya*, 66, 63--74.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- sn::rsn(200, xi = 0, omega = 1, alpha = 5)
#' sn.para.bootstrap.ks.test(x, B = 1000)
#' sn.para.bootstrap.cvm.test(x, B = 1000)
#' }
#'
#' @name sn-parametric-bootstrap
NULL


.validate_sn_bootstrap_args <- function(data, B, seed, verbose) {
  if (is.null(data) || !typeof(data) %in% c("integer", "double") ||
      !is.null(dim(data))) {
    stop("`data` must be a numeric vector.", call. = FALSE)
  }
  if (length(data) < 10L) {
    stop("`data` must contain at least 10 observations.", call. = FALSE)
  }
  if (any(!is.finite(data))) {
    stop("`data` must contain only finite values (no NA, NaN, or Inf).",
         call. = FALSE)
  }
  if (diff(range(data)) == 0) {
    stop("`data` must contain at least two distinct values.", call. = FALSE)
  }

  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B < 1 ||
      B > .Machine$integer.max || B != floor(B)) {
    stop("`B` must be a positive integer.", call. = FALSE)
  }
  B <- as.integer(B)

  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
       seed < 0 || seed > .Machine$integer.max || seed != floor(seed))) {
    stop("`seed` must be a single integer or NULL.", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- as.integer(seed)
  }

  if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  list(data = as.double(data), B = B, seed = seed, verbose = verbose)
}


.sn_fit_dp <- function(x) {
  fit <- suppressWarnings(
    tryCatch(
      sn.fit.robust(x, para_form = "DP"),
      error = function(e) NULL
    )
  )

  required <- c("xi", "omega", "alpha")
  if (is.null(fit) || !all(required %in% names(fit))) {
    return(NULL)
  }

  parameters <- unname(fit[required])
  if (!all(is.finite(parameters)) || parameters[2L] <= 0) {
    return(NULL)
  }

  stats::setNames(parameters, required)
}


.sn_ks_statistic <- function(x, parameters) {
  n <- length(x)
  probabilities <- sort(sn::psn(
    x,
    xi = parameters[["xi"]],
    omega = parameters[["omega"]],
    alpha = parameters[["alpha"]]
  ))

  if (any(!is.finite(probabilities))) {
    return(NA_real_)
  }

  index <- seq_len(n)
  d <- max(index / n - probabilities,
           probabilities - (index - 1L) / n)
  sqrt(n) * d
}


.sn_cvm_statistic <- function(x, parameters) {
  n <- length(x)
  probabilities <- sort(sn::psn(
    x,
    xi = parameters[["xi"]],
    omega = parameters[["omega"]],
    alpha = parameters[["alpha"]]
  ))

  if (any(!is.finite(probabilities))) {
    return(NA_real_)
  }

  expected <- (2 * seq_len(n) - 1) / (2 * n)
  1 / (12 * n) + sum((probabilities - expected)^2)
}


.sn_parametric_bootstrap_test <- function(data, B, seed, verbose,
                                          statistic, statistic_name) {
  args <- .validate_sn_bootstrap_args(data, B, seed, verbose)
  data <- args$data
  B <- args$B

  observed_fit <- .sn_fit_dp(data)
  if (is.null(observed_fit)) {
    stop("Skew-normal fitting failed for the observed data.", call. = FALSE)
  }

  observed_statistic <- statistic(data, observed_fit)
  if (!is.finite(observed_statistic)) {
    stop("The goodness-of-fit statistic could not be computed for the observed data.",
         call. = FALSE)
  }

  if (!is.null(args$seed)) {
    seed_existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (seed_existed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (seed_existed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(args$seed)
  }

  n <- length(data)
  bootstrap_statistics <- rep.int(NA_real_, B)
  for (b in seq_len(B)) {
    bootstrap_sample <- sn::rsn(
      n,
      xi = observed_fit[["xi"]],
      omega = observed_fit[["omega"]],
      alpha = observed_fit[["alpha"]]
    )

    bootstrap_fit <- .sn_fit_dp(bootstrap_sample)
    if (!is.null(bootstrap_fit)) {
      bootstrap_statistics[b] <- statistic(bootstrap_sample, bootstrap_fit)
    }
  }

  valid_statistics <- bootstrap_statistics[is.finite(bootstrap_statistics)]
  n_valid <- length(valid_statistics)
  n_failed <- B - n_valid
  if (n_valid == 0L) {
    stop("All bootstrap fits failed; the p-value cannot be computed.",
         call. = FALSE)
  }
  if (n_valid < B / 2) {
    warning(
      sprintf(
        "Only %.1f%% of bootstrap samples were valid; results may be unreliable.",
        100 * n_valid / B
      ),
      call. = FALSE
    )
  }

  p_value <- (sum(valid_statistics >= observed_statistic) + 1) /
    (n_valid + 1)
  attr(p_value, "statistic") <- observed_statistic
  attr(p_value, "B") <- B
  attr(p_value, "valid") <- n_valid
  attr(p_value, "failed") <- n_failed

  if (args$verbose) {
    cat(sprintf("Observed %s statistic = %.4f\n",
                statistic_name, observed_statistic))
    cat(sprintf("Valid bootstraps = %d/%d\n", n_valid, B))
    cat(sprintf("p-value = %.6g\n", p_value))
  }

  p_value
}


#' @rdname sn-parametric-bootstrap
#' @export
sn.para.bootstrap.ks.test <- function(data = NULL, B = 1000L, seed = 103L,
                                      verbose = FALSE) {
  .sn_parametric_bootstrap_test(
    data = data,
    B = B,
    seed = seed,
    verbose = verbose,
    statistic = .sn_ks_statistic,
    statistic_name = "KS"
  )
}


#' @rdname sn-parametric-bootstrap
#' @export
sn.para.bootstrap.cvm.test <- function(data = NULL, B = 1000L, seed = 103L,
                                       verbose = FALSE) {
  .sn_parametric_bootstrap_test(
    data = data,
    B = B,
    seed = seed,
    verbose = verbose,
    statistic = .sn_cvm_statistic,
    statistic_name = "CvM"
  )
}

