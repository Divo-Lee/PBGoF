# PBGoF

### PBGoF: Parametric Bootstrap Goodness-of-Fit Tests for the Skew-normal Distribution with Estimated Parameters

Hongxiang Li and Tsung Fei Khang

PBGoF: R package provides goodness-of-fit tests for the skew-normal distribution with estimated parameters. Implements Kolmogorov-Smirnov and Cramer-von Mises tests using parametric bootstrap or precomputed simulation quantiles, together with robust parameter estimation procedures.

### Installation:
Install PBGoF from local source with

`install.packages("PBGoF_0.1.0.tar.gz", repos=NULL, type="source")`

Install PBGoF from GitHub with

 `if (!"devtools" %in% installed.packages()) {
  install.packages("devtools")}`
  
 `devtools::install_github("Divo-Lee/PBGoF")`

 or

 `if (!"pak" %in% installed.packages()) {
  install.packages("pak")}`
  
 `pak::pkg_install("Divo-Lee/PBGoF")`
 
### Dependencies:
 `PBGoF` `R` package depends on the following packages: `sn`, `methods`.
