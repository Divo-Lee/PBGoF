# PBGoF

### PBGoF (version 0.2.0): Parametric Bootstrap Goodness-of-Fit Tests for Skew-Normal and Skew-t Distributions

Hongxiang Li, Chenglin Xu, and Tsung Fei Khang

PBGoF provides goodness-of-fit procedures, robust parameter fitting, and
graphical diagnostics for univariate skew-normal and skew-t distributions.
Version 0.1.0 focused on skew-normal models, including robust fitting,
precomputed-quantile tests, and parametric-bootstrap KS and CvM tests. Version
0.2.0 extends the package to skew-t models through `st.fit.robust()`,
skew-t parametric-bootstrap KS and CvM tests, and the unified
`para.bootstrap.test()` interface, which allows users to select either the SN
or ST model. It also adds `st.plot.check()` and the model-selectable
`cdf.plot.check()` and `QQ.plot()` diagnostics for SN and ST models.

## Installation

Install PBGoF from CRAN with:

```r
install.packages("PBGoF")
```

Install PBGoF 0.2.0 from a local source archive with:

```r
install.packages("PBGoF_0.2.0.tar.gz", repos = NULL, type = "source")
```

The original PBGoF GitHub installation instructions remain applicable to PBGoF itself:

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

For skew-t models, PBGoF provides `st.fit.robust()`,
`para.bootstrap.test()`, and `st.plot.check()`.

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

PBGoF estimates the model through `sn.fit.robust()`: MLE is attempted first,
followed by increasingly stabilized penalized fits. A fit is accepted only if
the estimates and standard errors are finite and the scale is positive.

### Skew-t model

The skew-t family extends the skew-normal model by adding a positive
degrees-of-freedom parameter $\nu$. Its direct parameterization is
$(\xi,\omega,\alpha,\nu)$. The extra parameter allows asymmetry and heavy tails
to be represented together; as $\nu$ becomes large, the skew-t model approaches
the skew-normal model.

`st.fit.robust()` first attempts ordinary MLE and then uses Q-penalized MPLE as
a fallback. Its DP output contains `xi`, `omega`, `alpha`, `nu` and their
standard errors. CP output is available when the corresponding moments exist.

### EDF statistics with estimated parameters

For order statistics $X_{(1)}\leq\cdots\leq X_{(n)}$, define

$$
U_{(i)}=F_{\mathrm{SN}} \left(X_{(i)};\widehat{\boldsymbol\theta}\right),
\qquad i=1,\ldots,n.
$$

The two-sided Kolmogorov-Smirnov discrepancy and the statistic used by PBGoF are

$$
D_n=\max_{1\leq i\leq n}\left(\frac{i}{n}-U_{(i)},\;U_{(i)}-\frac{i-1}{n}\right),
\qquad T_{\mathrm{KS}}=\sqrt{n} D_n.
$$

The Cramer-von Mises statistic is

$$
W_n^2=\frac{1}{12n}+\sum_{i=1}^{n}\left[U_{(i)}-\frac{2i-1}{2n}\right]^2.
$$

The parametric-bootstrap CvM test uses $W_n^2$. Within the tabulated sample-size range, the precomputed-quantile CvM test uses

$$
T_{\mathrm{CvM}}=\sqrt{n} W_n^2.
$$

### Parametric-bootstrap calibration

The unified function `para.bootstrap.test()` selects the model through
`model = "SN"` or `model = "ST"`, and selects the statistic through
`statistic = "KS"` or `statistic = "CvM"`. For example,
`para.bootstrap.test(x, model = "SN", statistic = "KS")` and
`para.bootstrap.test(x, model = "ST", statistic = "CvM")` calibrate the
corresponding tests for an estimated model. The previous model- and
statistic-specific functions remain available for compatibility. Both model
families use the same full-estimation bootstrap principle:

1. Fit the observed sample and compute $T_{\mathrm{obs}}$.
2. Generate $B$ samples of size $n$ from the fitted SN or skew-t model.
3. Re-estimate all parameters independently in each bootstrap sample.
4. Compute the same statistic with that sample's fitted parameters.

For $B_{\mathrm{valid}}$ finite bootstrap statistics, PBGoF uses

$$
\widehat p_{\mathrm{boot}}=\frac{1+\displaystyle\sum_{b=1}^{B_{\mathrm{valid}}}\mathbf{1}  \left(T_b^*\geq T_{\mathrm{obs}}\right)}{B_{\mathrm{valid}}+1}.
$$

Failed fits are excluded and reported. Re-estimation in every replicate
calibrates the statistic for the composite null rather than incorrectly
treating the fitted distribution as fixed. For skew-t, all four parameters,
including $\nu$, are re-estimated in every valid replicate.

The precomputed skew-normal tests remain available through
`PBGoF_ks_test()` and `PBGoF_cvm_test()`.

### Precomputed-quantile calibration for skew-normal data

`PBGoF_ks_test()` and `PBGoF_cvm_test()` use tables generated from 100,000 Monte Carlo replicates per available sample-size and centered-skewness combination. The data are fitted in two parameterizations:

- DP $(\xi,\omega,\alpha)$ evaluates the fitted CDF and EDF statistic.
- CP $(\mu,\sigma,\gamma_1)$ matches skewness to the simulation table.

The lookup value is

$$
\gamma_{1,\mathrm{used}}=\min \left(0.99,\max  \left[0.01,\mathrm{round}  \left(\left|\widehat\gamma_1\right|,2\right)\right]\right).
$$

The absolute value follows reflection symmetry: if $X\sim\mathrm{SN}(\xi,\omega,\alpha)$, then $-X$ has shape $-\alpha$. The signs of $\alpha$ and $\gamma_1$ reverse, but the null distributions of the reflection-invariant EDF statistics do not. Therefore, estimated skewness values with the same absolute magnitude but opposite signs use the same reference-table row. PBGoF retains the signed `gamma1_hat` and returns the non-negative lookup value as `gamma1_used`.

The bundled tables cover sample sizes through 500. For $n>500$, all observations remain in the fit and EDF, while `n_used = 500` is used only to select the table row. The KS test compares $\sqrt{n}D_n$ with the stored $n=500$ quantiles. For CvM, the stored $n=500$ quantiles are divided by $\sqrt{500}$ and compared with the observed $W_n^2$. This follows the approximation of using the $n=500$ critical values for larger samples.

Let $q_p^{\dagger}$ denote the selected quantiles expressed on the same scale as the observed statistic; for CvM with $n>500$, $q_p^{\dagger}=q_p/\sqrt{500}$. The table-based p-value is

$$
p^{\*}=\min_{q_p^{\dagger}\geq T_{\mathrm{obs}}}p,\qquad \widehat p_{\mathrm{table}}=1-p^{\*}.
$$

The probability grid is 0.01 to 0.99 in increments of 0.01, so the result is a conservative step-function approximation. PBGoF does not interpolate across sample size or skewness.

### Graphical checks and interpretation

A small p-value is evidence against the fitted model. A large p-value does not
prove skew-normality or skew-t adequacy; it means that the test did not detect
a departure at the available sample size and calibration resolution. Use
`sn.plot.check()` for SN fits and `st.plot.check()` for skew-t fits alongside
the corresponding formal tests. The precomputed-quantile functions remain
specific to the skew-normal family.

## References

- Azzalini, A. (1985). A class of distributions which includes the normal ones. *Scandinavian Journal of Statistics*, **12**(2), 171-178.
- Azzalini, A., and Capitanio, A. (2014). *The Skew-Normal and Related Families*. Cambridge University Press.
- Babu, G. J., and Rao, C. R. (2004). Goodness-of-fit tests when parameters are estimated. *Sankhya*, **66**, 63-74.
- Mateu-Figueras, G., Puig, P., and Pewsey, A. (2007). Goodness-of-fit tests for the skew-normal distribution when the parameters are estimated from the data. *Communications in Statistics-Theory and Methods*, **36**(9), 1735-1755. https://doi.org/10.1080/03610920601126217
