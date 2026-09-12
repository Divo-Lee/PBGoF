# PBGoF

### PBGoF: Parametric Bootstrap Goodness-of-Fit Tests for the Skew-normal Distribution with Estimated Parameters

Hongxiang Li and Tsung Fei Khang

PBGoF is an R package for assessing whether a numeric sample is compatible with a univariate skew-normal distribution when the model parameters are estimated from the same data. It provides Kolmogorov-Smirnov (KS) and Cramer-von Mises (CvM) tests using either a parametric bootstrap or precomputed simulation quantiles, together with robust parameter estimation procedures.

## Installation

Install PBGoF from a local source archive with:

```r
install.packages("PBGoF_0.1.0.tar.gz", repos = NULL, type = "source")
```

Install PBGoF from GitHub with `devtools`:

```r
if (!"devtools" %in% rownames(installed.packages())) {
  install.packages("devtools")
}
devtools::install_github("Divo-Lee/PBGoF")
```

or with `pak`:

```r
if (!"pak" %in% rownames(installed.packages())) {
  install.packages("pak")
}
pak::pkg_install("Divo-Lee/PBGoF")
```

## Dependencies

PBGoF depends on the R packages `sn` and `methods`.

## Methods

### Skew-normal model and composite null hypothesis

Let $X_1,\ldots,X_n$ be an independent sample. Under the null hypothesis,

$$
H_0:\quad X_i \overset{\mathrm{iid}}{\sim} \mathrm{SN}(\xi,\omega,\alpha),
\qquad \xi\in\mathbb{R},\quad \omega>0,\quad \alpha\in\mathbb{R}.
$$

Writing $z=(x-\xi)/\omega$, the skew-normal density is

$$
f_{\mathrm{SN}}(x;\xi,\omega,\alpha)=\frac{2}{\omega} \phi(z) \Phi(\alpha z),
$$

where $\phi$ and $\Phi$ are the standard normal density and distribution functions. Because $\boldsymbol\theta=(\xi,\omega,\alpha)^{\mathsf T}$ is estimated from the sample, this is a composite goodness-of-fit problem; the usual KS null distribution for a fully specified model is not applicable.

PBGoF estimates the model through `sn.fit.robust()`: MLE is attempted first, followed by MPLE and then MPLE with the matching-prior penalty. A fit is accepted only if the estimates and standard errors are finite and the scale is positive.

### EDF statistics with estimated parameters

For order statistics $X_{(1)}\leq\cdots\leq X_{(n)}$, define

$$
U_{(i)}=F_{\mathrm{SN}} \left(X_{(i)};\widehat{\boldsymbol\theta}\right),
\qquad i=1,\ldots,n.
$$

The two-sided Kolmogorov-Smirnov discrepancy and PBGoF scaling are

$$
D_n=\max_{1\leq i\leq n}\left(\frac{i}{n}-U_{(i)},\;U_{(i)}-\frac{i-1}{n}\right),
\qquad T_{\mathrm{KS}}=\sqrt{n_{\mathrm{eff}}}\,D_n.
$$

The Cramer-von Mises statistic is

$$
W_n^2=\frac{1}{12n}+\sum_{i=1}^{n}\left[U_{(i)}-\frac{2i-1}{2n}\right]^2.
$$

The parametric-bootstrap CvM test uses $W_n^2$, whereas the precomputed-quantile CvM test uses

$$
T_{\mathrm{CvM}}=\sqrt{n_{\mathrm{eff}}}\,W_n^2.
$$

Within the table range, $n_{\mathrm{eff}}=n$.

### Parametric-bootstrap calibration

The functions `sn.para.bootstrap.ks.test()` and `sn.para.bootstrap.cvm.test()` reproduce the full estimation procedure:

1. Fit the observed sample and compute $T_{\mathrm{obs}}$.
2. Generate $B$ samples of size $n$ from $\mathrm{SN}(\widehat\xi,\widehat\omega,\widehat\alpha)$.
3. Re-estimate all parameters independently in each bootstrap sample.
4. Compute the same statistic with that sample's fitted parameters.

For $B_{\mathrm{valid}}$ finite bootstrap statistics, PBGoF uses

$$
\widehat p_{\mathrm{boot}}=\frac{1+\displaystyle\sum_{b=1}^{B_{\mathrm{valid}}}\mathbf{1}  \left(T_b^*\geq T_{\mathrm{obs}}\right)}{B_{\mathrm{valid}}+1}.
$$

Failed fits are excluded and reported. Re-estimation in every replicate calibrates the statistic for the composite null rather than incorrectly treating the fitted distribution as fixed.

### Precomputed-quantile calibration

`PBGoF_ks_test()` and `PBGoF_cvm_test()` use tables generated from 100,000 Monte Carlo replicates per available sample-size and centered-skewness combination. The data are fitted in two parameterizations:

- DP $(\xi,\omega,\alpha)$ evaluates the fitted CDF and EDF statistic.
- CP $(\mu,\sigma,\gamma_1)$ matches skewness to the simulation table.

The lookup value is

$$
\gamma_{1,\mathrm{used}}=\min \left(0.99,\max  \left[0.01,\mathrm{round}  \left(\left|\widehat\gamma_1\right|,2\right)\right]\right).
$$

The absolute value follows reflection symmetry: if $X\sim\mathrm{SN}(\xi,\omega,\alpha)$, then $-X$ has shape $-\alpha$. The signs of $\alpha$ and $\gamma_1$ reverse, but the null distributions of the reflection-invariant EDF statistics do not. Therefore, estimated skewness values with the same absolute magnitude but opposite signs use the same reference-table row. PBGoF retains the signed `gamma1_hat` and returns the non-negative lookup value as `gamma1_used`.

The bundled tables cover sample sizes through 500. For $n>500$, all observations remain in the fit and EDF, but

$$
n_{\mathrm{eff}}=\min(n,500)=500
$$

is used for the external statistic multiplier and table lookup, following the finding that fitted-skew-normal EDF critical values above 500 are nearly identical to those at 500.

For stored quantiles $q_p$, define

$$
p^{\*}=\min_{q_p\geq T_{\mathrm{obs}}}p,\qquad \widehat p_{\mathrm{table}}=1-p^{\*}.
$$

The probability grid is 0.01 to 0.99 in increments of 0.01, so the result is a conservative step-function approximation. PBGoF does not interpolate across sample size or skewness.

### Interpretation

A small p-value is evidence against the fitted skew-normal model. A large p-value does not prove skew-normality; it means the test did not detect a departure at the available sample size and calibration resolution. Use `sn.plot.check()` alongside the formal tests.

## References

- Azzalini, A. (1985). A class of distributions which includes the normal ones. *Scandinavian Journal of Statistics*, **12**(2), 171-178.
- Azzalini, A., and Capitanio, A. (2014). *The Skew-Normal and Related Families*. Cambridge University Press.
- Babu, G. J., and Rao, C. R. (2004). Goodness-of-fit tests when parameters are estimated. *Sankhya*, **66**, 63-74.
- Mateu-Figueras, G., Puig, P., and Pewsey, A. (2007). Goodness-of-fit tests for the skew-normal distribution when the parameters are estimated from the data. *Communications in Statistics-Theory and Methods*, **36**(9), 1735-1755. https://doi.org/10.1080/03610920601126217
