# PBGoF

### PBGoF: Parametric Bootstrap Goodness-of-Fit Tests for the Skew-normal Distribution with Estimated Parameters

Hongxiang Li and Tsung Fei Khang

PBGoF is an R package for assessing whether a numeric sample is compatible with a univariate skew-normal distribution when the model parameters are estimated from the same data. It provides Kolmogorov-Smirnov (KS) and Cramér-von Mises (CvM) tests using either a parametric bootstrap or precomputed simulation quantiles, together with robust parameter estimation procedures.

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

Let $X_1,\ldots,X_n$ be an independent sample. Under the null hypothesis, the observations follow a univariate skew-normal distribution,

$$
H_0:\quad X_i \overset{\mathrm{iid}}{\sim} \operatorname{SN}(\xi,\omega,\alpha),
\qquad \xi\in\mathbb{R},\quad \omega>0,\quad \alpha\in\mathbb{R}.
$$

Writing \(z=(x-\xi)/\omega\), its density is

$$
f_{\operatorname{SN}}(x;\xi,\omega,\alpha)
=\frac{2}{\omega}\,\phi(z)\,\Phi(\alpha z),
$$

where \(\phi\) and \(\Phi\) denote the standard normal density and distribution functions. Because \(\boldsymbol\theta=(\xi,\omega,\alpha)^{\mathsf T}\) is estimated from the sample, this is a composite goodness-of-fit problem; therefore, the usual KS distribution for a completely specified null model is not applicable.

PBGoF estimates the model through `sn.fit.robust()`. The function first attempts maximum likelihood estimation (MLE), then maximum penalized likelihood estimation (MPLE), and finally MPLE with the matching-prior penalty. A fitted model is accepted only when its parameter estimates and standard errors are finite and its scale estimate is positive.

### EDF statistics with estimated parameters

Let \(X_{(1)}\leq\cdots\leq X_{(n)}\) denote the order statistics and let

$$
U_{(i)}=F_{\operatorname{SN}}\!\left(X_{(i)};\widehat{\boldsymbol\theta}\right),
\qquad i=1,\ldots,n,
$$

where \(F_{\operatorname{SN}}\) is the skew-normal distribution function evaluated at the fitted direct parameters \(\widehat{\boldsymbol\theta}=(\widehat\xi,\widehat\omega,\widehat\alpha)^{\mathsf T}\).

The two-sided Kolmogorov-Smirnov discrepancy is

$$
D_n=\max_{1\leq i\leq n}
\left {
\frac{i}{n}-U_{(i)},\;
U_{(i)}-\frac{i-1}{n}
\right }.
$$

PBGoF reports the scaled KS statistic

$$
T_{\mathrm{KS}}=\sqrt{n_{\mathrm{eff}}} D_n.
$$

The Cramér-von Mises statistic is

$$
W_n^2=\frac{1}{12n}+
\sum_{i=1}^{n}
\left[
U_{(i)}-\frac{2i-1}{2n}
\right]^2.
$$

The parametric-bootstrap CvM test uses \(W_n^2\), whereas the precomputed-quantile CvM test uses the table-compatible scaling

$$
T_{\mathrm{CvM}}=\sqrt{n_{\mathrm{eff}}}\,W_n^2.
$$

For an ordinary sample within the table range, \(n_{\mathrm{eff}}=n\).

### Parametric-bootstrap calibration

The functions `sn.para.bootstrap.ks.test()` and `sn.para.bootstrap.cvm.test()` account for parameter estimation by reproducing the complete fitting procedure in every bootstrap sample:

1. Fit the skew-normal model to the observed data and compute \(T_{\mathrm{obs}}\).
2. Generate \(B\) independent samples of size \(n\) from \(\operatorname{SN}(\widehat\xi,\widehat\omega,\widehat\alpha)\).
3. Re-estimate all skew-normal parameters separately in every bootstrap sample.
4. Compute the same EDF statistic, using that bootstrap sample's fitted parameters.

If \(B_{\mathrm{valid}}\) bootstrap fits produce finite statistics \(T_1^*,\ldots,T_{B_{\mathrm{valid}}}^*\), PBGoF uses the finite-simulation correction

$$
\widehat p_{\mathrm{boot}}
=\frac{1+\displaystyle\sum_{b=1}^{B_{\mathrm{valid}}}
\mathbf{1}\!\left(T_b^*\geq T_{\mathrm{obs}}\right)}
{B_{\mathrm{valid}}+1}.
$$

Failed fits are excluded and reported. Re-estimating the parameters in every replicate is essential: it calibrates the statistic for the composite null hypothesis rather than treating the fitted distribution as fixed.

### Precomputed-quantile calibration

The functions `PBGoF_ks_test()` and `PBGoF_cvm_test()` provide a faster alternative based on tables generated from 100,000 Monte Carlo replicates for each available combination of sample size and centered skewness. The observed data are fitted twice:

- the direct parameterization (DP), \((\xi,\omega,\alpha)\), is used to evaluate the fitted distribution and calculate the EDF statistic;
- the centered parameterization (CP), \((\mu,\sigma,\gamma_1)\), is used to match the estimated skewness to the simulation table.

The lookup skewness is

$$
\gamma_{1,\mathrm{used}}
=\min\!\left\{0.99,
\max\!\left[0.01,
\operatorname{round}\!\left(\left|\widehat\gamma_1\right|,2\right)
\right]\right\}.
$$

Taking the absolute value is justified by reflection symmetry. If \(X\sim\operatorname{SN}(\xi,\omega,\alpha)\), then \(-X\) has shape \(-\alpha\); the signs of \(\alpha\) and \(\gamma_1\) reverse, but the null distributions of the reflection-invariant EDF statistics do not change. Consequently, fitted skewness values \(-g\) and \(+g\) use the same reference row. PBGoF retains the signed estimate as `gamma1_hat` and reports the non-negative lookup value as `gamma1_used`.

The bundled tables cover sample sizes through 500. For \(n>500\), all observations remain in the parameter fit and empirical distribution function, but PBGoF sets

$$
n_{\mathrm{eff}}=\min(n,500)=500
$$

for the external statistic multiplier and the table lookup. This follows the finding that EDF critical values above \(n=500\) are nearly identical to those at \(n=500\) for the fitted skew-normal model.

For the selected \((n_{\mathrm{eff}},\gamma_{1,\mathrm{used}})\) row, let \(q_p\) be the stored \(p\)-quantile and define

$$
p^{*}=\min\left\{p:q_p\geq T_{\mathrm{obs}}\right\}.
$$

The lookup test returns the upper-tail approximation

$$
\widehat p_{\mathrm{table}}=1-p^{*}.
$$

Because the stored probability grid runs from 0.01 to 0.99 in increments of 0.01, this p-value is a conservative step-function approximation and has no more precision than the table grid. PBGoF does not interpolate across sample size or skewness.

### Interpretation

A small p-value indicates that the observed EDF discrepancy is unusually large under the fitted skew-normal model. It is evidence against that model, not a measure of the practical size of the discrepancy. A large p-value does not prove skew-normality; it means that the selected test did not detect a departure at the available sample size and calibration resolution. Graphical diagnostics from `sn.plot.check()` should be used alongside the formal tests.

## References

- Azzalini, A. (1985). A class of distributions which includes the normal ones. *Scandinavian Journal of Statistics*, **12**(2), 171-178.
- Azzalini, A., and Capitanio, A. (2014). *The Skew-Normal and Related Families*. Cambridge University Press.
- Babu, G. J., and Rao, C. R. (2004). Goodness-of-fit tests when parameters are estimated. *Sankhya*, **66**, 63-74.
- Mateu-Figueras, G., Puig, P., and Pewsey, A. (2007). Goodness-of-fit tests for the skew-normal distribution when the parameters are estimated from the data. *Communications in Statistics-Theory and Methods*, **36**(9), 1735-1755. https://doi.org/10.1080/03610920601126217
