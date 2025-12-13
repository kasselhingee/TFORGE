Author: Kassel Hingee

This directory contains an R package, `TFORGE`, that implements the hypothesis tests for eigenvalues in Hingee, Scealy and Wood (2026, Nonparametric bootstrap inference for the eigenvalues of geophysical tensors, accepted by the *Journal of American Statistical Association*).
A version of this package is available from CRAN at `https://cran.r-project.org/package=TFORGE`.

## Installation
Install this package from github by running:
```
devtools::install_github("https://github.com/kasselhingee/TFORGE")
```
or install from CRAN by running:
```
install.packages("TFORGE")
```

## About
All tests in this package are for hypotheses about the eigenvalues of the (extrinsic) mean of the sampled population(s).
Functions for conducting hypothesis tests start with `test_`.
When the matrices are constrained to a submanifold of the space of symmetric matrices we use the extrinsic mean, which projects the usual linear/Euclidean mean so that it satisfies the same constraints as the data.
Tests can be calibrated using either a chi-squared distribution or bootstrapping; in simulations bootstrapping was more reliable but slower.


The following functions conduct single sample and k-sample tests of the eigenvalues of the (extrinsic) population mean:

+ `test_unconstrained()`
+ `test_fixedtrace()` when matrices have fixed trace
+ `test_ss1()` when the squared eigenvalues of each matrix sums to 1.
+ `test_ss1fixedtrace()` when the squared eigenvalues of each matrix sums to 1 and the trace is zero.

The single sample tests conducted by the above functions test the null hypothesis of a user-provided set of eigenvalues for the extrinsic population mean against the alternative hypothesis that the extrinsic population mean has different eigenvalues.
The k-sample tests conducted by the above functions test the null hypothesis that the extrinsic population means have the same eigenvalues.

Additionally `test_unconstrained_aGOE()` can perform 2-sample tests with calibration by a Gaussian Orthogonal Ensemble (GOE) approximation (Schwartzman et al., 2010, Group comparison of eigenvalues and eigenvectors of diffusion tensors. *Journal of the American Statistical Association*).
 A bootstrapped calibration for this test is also available where the GOE approximation is used to stabilise the scale of the statistic.

There are two functions `conf_fixedtrace()` and `conf_ss1fixedtrace()` for estimating confidence regions.

The above tests all require that the eigenvalues of the population mean are distinct (with degraded performance when eigenvalues are very close to each other).
Eigenvalues are assumed to be in descending order.

Use `test_multiplicity()` to test the eigenvalue-multiplicity of the population mean of a single sample.
A test of the same hypothesis that requires that matrix elements follow a multivariate Gaussian distribution with orthogonally-invariant covariance is also available through `test_multiplicity_OI()` (Schwartzman et al., 2008, Inference for eigenvalues and eigenvectors of Gaussian symmetric matrices, *The Annals of Statistics*).
`test_multiplicity()` can also be applied to matrices with a constrained trace, but use `test_multiplicity_nonnegative()` for matrices constrained to have non-negative eigenvalues.


## Usage
To use the functions starting with `test_`, you must have symmetric matrix data formatted to be suitable for `as_fsm()` or `as_kfsm()` (please see the help for these functions for more detail).
The package includes a vignette demonstrating an application to anisotropy of magnetic susceptibility data.
For example applications of all the hypothesis tests in this package, please see the reproducibility document associated with (Hingee, Scealy and Wood, 2026, Nonparametric bootstrap inference for the eigenvalues of geophysical tensors, accepted by the *Journal of American Statistical Association*).
Below is a toy example using symmetric matrices simulated from a multivariate Normal (aka Gaussian) distribution.

```
# simulate data
flattened_matrices <- rsymm_norm(15, diag(c(3, 2, 1)))
# test multiplicity of eigenvalues of population mean
test_multiplicity(flattened_matrices, mult = c(1, 2))
test_multiplicity(flattened_matrices, mult = c(2, 1))
# test eigenvalues of population mean
test_unconstrained(flattened_matrices, evals = c(3, 2, 1))
```

## Contributing
Contributions welcome.
Please feel free to contact the author of this package and/or make pull requests.

The package includes substantial unit testing. 
When the `FAST_CHECK` environmental variable is `"true"`, then fast unit testing is conducted.
Package checks set `FAST_CHECK` based on the `NOT_CRAN` environmental variable (which is set by the `devtools` package).
These fast unit tests can also be interactively by running `Sys.setenv(FAST_CHECK = "true")` before running `devtools::test()` or similar.
More thorough and slower unit tests can be conducted by running `devtools::test()` or `testthat::test_local()` on the package without setting `FAST_CHECK`.
Note that the bespoke `FAST_CHECK` variable is used to enable fast interactive testing because `testthat::skip_on_cran()` currently never skips interactively, regardless of the value of `NOT_CRAN`.



