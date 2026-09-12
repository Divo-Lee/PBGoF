\name{PBGoF-precomputed-tests}
\alias{PBGoF_ks_test}
\alias{PBGoF_cvm_test}
\title{Fast Skew-Normal Goodness-of-Fit Tests Using Precomputed Quantiles}
\usage{
PBGoF_ks_test(data, ks_table = NULL)

PBGoF_cvm_test(data, cvm_table = NULL)
}
\arguments{
\item{data}{A numeric vector whose length is represented in the table.}
\item{ks_table, cvm_table}{Optional custom quantile tables. If \code{NULL},
the corresponding table bundled with PBGoF is used.}
}
\value{
A list containing \code{statistic}, actual sample size \code{n}, lookup and
scaling sample size \code{n_used}, \code{gamma1_hat}, \code{gamma1_used}, and
\code{p.value}.
}
\description{
Computes a skew-normal goodness-of-fit statistic and obtains an approximate
p-value from a precomputed 100,000-replicate quantile table.
}
\details{
The fitted absolute CP skewness is rounded to two decimal places and truncated
to \code{[0.01, 0.99]}. P-values are conservative step-function approximations
between 0.01 and 0.99. No interpolation is performed.

If the fitted \code{gamma1} is negative, table lookup uses
\code{abs(gamma1)}. Reflection of a skew-normal variable changes the signs of
the DP shape parameter and CP skewness, but the null distributions of the EDF
statistics are invariant to this sign change (Mateu-Figueras et al., 2007).
The bundled tables therefore require only non-negative skewness values.

For samples larger than 500, all observations remain in the fit and empirical
distribution function, while the external statistic multiplier and table row
use \code{n_used = 500}.
}
\references{
Mateu-Figueras, G., Puig, P., and Pewsey, A. (2007). Goodness-of-fit tests for
the skew-normal distribution when the parameters are estimated from the data.
\emph{Communications in Statistics---Theory and Methods}, 36(9), 1735--1755.
\doi{10.1080/03610920601126217}
}
\examples{
\donttest{
set.seed(123)
x <- sn::rsn(50, xi = 0, omega = 1, alpha = 3)
PBGoF_ks_test(x)
PBGoF_cvm_test(x)
}
}
