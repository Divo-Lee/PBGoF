#' Fast skew-normal goodness-of-fit tests using precomputed quantiles
#'
#' @description
#' Compute a skew-normal goodness-of-fit statistic and obtain an approximate
#' p-value from the package's precomputed 100,000-replicate quantile tables.
#' `PBGoF_ks_test()` uses the Kolmogorov--Smirnov statistic and
#' `PBGoF_cvm_test()` uses the Cramér--von Mises statistic.
#'
#' @param data A numeric vector. Its length must be represented in the
#'   selected quantile table (30 through 500 in the bundled tables).
#' @param ks_table,cvm_table An optional custom quantile table with columns
#'   `n`, `gamma1`, and `q_0.01` through `q_0.99`. If `NULL`, the corresponding
#'   table bundled with PBGoF is used.
#'
#' @return A list with components `statistic`, `n`, `n_used`, `gamma1_hat`,
#'   `gamma1_used`, and `p.value`. `n` is the actual sample size and `n_used`
#'   is the sample size used for statistic scaling and table lookup.
#'
#' @details
#' The observed distribution is fitted in DP form for the test statistic and
#' in CP form for matching the absolute skewness to the simulation table.
#' Both fits use [sn.fit.robust()]. A negative fitted gamma1 is deliberately
#' matched to the table using its absolute value. For the skew-normal family,
#' reflecting the variable changes the signs of the shape parameter alpha and
#' centered skewness gamma1, but not the null distribution of the EDF
#' goodness-of-fit statistics. Mateu-Figueras et al. (2007) therefore use
#' positive shape values when the fitted shape is negative. The bundled PBGoF
#' tables express the same symmetry in the CP scale and are indexed by absolute
#' gamma1. The absolute fitted skewness is rounded to
#' two decimal places and truncated to `[0.01, 0.99]`.
#'
#' For samples larger than 500, all observations are retained for fitting and
#' for constructing the empirical distribution function, but `n_used` is set
#' to 500. Consequently, the external statistic multiplier and quantile-table
#' lookup both use 500. This convention matches the scaling used to construct
#' the bundled tables. Mateu-Figueras et al. (2007) found that, for sample sizes
#' above 500, the quantiles of the EDF statistics were almost identical to
#' those at 500 and recommended using the `n = 500` critical values.
#'
#' Because the tables store percentiles in one-percentage-point increments,
#' the returned p-value is a conservative step-function approximation between
#' 0.01 and 0.99. No interpolation across sample size or skewness is performed.
#'
#' @references
#' Mateu-Figueras, G., Puig, P., and Pewsey, A. (2007). Goodness-of-fit tests
#' for the skew-normal distribution when the parameters are estimated from the
#' data. *Communications in Statistics---Theory and Methods*, 36(9),
#' 1735--1755. \doi{10.1080/03610920601126217}
#'
#' @examples
#' \donttest{
#' set.seed(12345)
#' x <- sn::rsn(200, xi = 0, omega = 1, alpha = 5)
#' PBGoF_ks_test(x)
#' PBGoF_cvm_test(x)
#' }
#'
#' @name PBGoF-precomputed-tests
NULL


.PBGoF_load_quantile_table <- function(table, statistic) {
  if (is.null(table)) {
    table <- switch(
      statistic,
      KS = ks_bootstrap_quantile_table_full,
      CvM = cvm_bootstrap_quantile_table_full,
      stop("Unknown statistic: ", statistic, call. = FALSE)
    )
  }

  if (!is.data.frame(table)) {
    stop("The quantile table must be a data frame.", call. = FALSE)
  }

  q_columns <- grep("^q_[0-9]+(?:\\.[0-9]+)?$", names(table), value = TRUE)
  required <- c("n", "gamma1")
  if (!all(required %in% names(table)) || length(q_columns) == 0L) {
    stop(
      "The quantile table must contain `n`, `gamma1`, and `q_*` columns.",
      call. = FALSE
    )
  }

  probabilities <- suppressWarnings(as.numeric(sub("^q_", "", q_columns)))
  if (any(!is.finite(probabilities)) || any(probabilities <= 0) ||
      any(probabilities >= 1) || anyDuplicated(probabilities)) {
    stop("The quantile table has invalid probability column names.",
         call. = FALSE)
  }

  order_index <- order(probabilities)
  list(
    table = table,
    q_columns = q_columns[order_index],
    probabilities = probabilities[order_index]
  )
}


.PBGoF_table_test <- function(data, table, statistic, statistic_function) {
  # sn.fit.robust() performs the complete numeric-data validation.
  fit_dp <- sn.fit.robust(data, para_form = "DP")
  dp <- fit_dp[c("xi", "omega", "alpha")]
  if (any(!is.finite(dp))) {
    stop("Skew-normal DP fitting failed for the observed data.", call. = FALSE)
  }

  fit_cp <- sn.fit.robust(data, para_form = "CP")
  gamma_hat <- unname(fit_cp[["gamma1"]])
  if (!is.finite(gamma_hat)) {
    stop("Skew-normal CP fitting failed for the observed data.", call. = FALSE)
  }

  n <- length(data)
  table_info <- .PBGoF_load_quantile_table(table, statistic)
  quantile_table <- table_info$table

  if (!is.numeric(quantile_table$n) || !is.numeric(quantile_table$gamma1)) {
    stop("The `n` and `gamma1` columns must be numeric.", call. = FALSE)
  }

  maximum_table_n <- max(quantile_table$n)
  n_used <- min(n, maximum_table_n)
  observed_statistic <- statistic_function(data, dp, n_used)
  if (!is.finite(observed_statistic)) {
    stop("The observed test statistic could not be computed.", call. = FALSE)
  }

  # The SN EDF statistics are invariant to reflection: changing the sign of
  # alpha, and therefore gamma1, does not change their null distributions
  # (Mateu-Figueras et al., 2007, Sections 2 and 4.1). The simulation tables
  # therefore need only the non-negative half of the gamma1 parameter space.
  gamma_used <- min(max(round(abs(gamma_hat), 2L), 0.01), 0.99)

  row_index <- which(
    quantile_table$n == n_used &
      abs(quantile_table$gamma1 - gamma_used) < 1e-8
  )
  if (length(row_index) == 0L) {
    available_n <- sort(unique(quantile_table$n))
    if (!n_used %in% available_n) {
      stop(
        sprintf(
          "No %s quantiles are available for n=%d; available n values range from %s to %s.",
          statistic, n_used, min(available_n), max(available_n)
        ),
        call. = FALSE
      )
    }
    stop(
      sprintf(
        "No %s quantiles are available for n=%d and |gamma1|=%.2f.",
        statistic, n_used, gamma_used
      ),
      call. = FALSE
    )
  }
  if (length(row_index) > 1L) {
    stop("The quantile table contains duplicate `n`/`gamma1` rows.",
         call. = FALSE)
  }

  quantiles <- as.numeric(
    quantile_table[row_index, table_info$q_columns, drop = TRUE]
  )
  if (any(!is.finite(quantiles)) || is.unsorted(quantiles)) {
    stop("The selected row contains invalid or unsorted quantiles.",
         call. = FALSE)
  }

  first_at_or_above <- which(quantiles >= observed_statistic)[1L]
  if (is.na(first_at_or_above)) {
    p_value <- 1 - table_info$probabilities[length(table_info$probabilities)]
  } else {
    p_value <- 1 - table_info$probabilities[first_at_or_above]
  }

  list(
    statistic = unname(observed_statistic),
    n = n,
    n_used = n_used,
    gamma1_hat = gamma_hat,
    gamma1_used = gamma_used,
    p.value = unname(p_value)
  )
}


#' @rdname PBGoF-precomputed-tests
#' @export
PBGoF_ks_test <- function(data, ks_table = NULL) {
  .PBGoF_table_test(
    data = data,
    table = ks_table,
    statistic = "KS",
    statistic_function = function(x, parameters, n_used) {
      # .sn_ks_statistic() uses sqrt(length(x)); only replace this external
      # multiplier. The EDF itself must continue to use every observation.
      .sn_ks_statistic(x, parameters) * sqrt(n_used / length(x))
    }
  )
}


#' @rdname PBGoF-precomputed-tests
#' @export
PBGoF_cvm_test <- function(data, cvm_table = NULL) {
  .PBGoF_table_test(
    data = data,
    table = cvm_table,
    statistic = "CvM",
    statistic_function = function(x, parameters, n_used) {
      sqrt(n_used) * .sn_cvm_statistic(x, parameters)
    }
  )
}
