#' @details
#' Functions for conducting hypothesis tests start with `test_`.
#' The following can conduct single sample and \eqn{k}-sample tests of eigenvalues using pivotal bootstap methods:
#' + [`test_unconstrained()`]
#' + [`test_fixedtrace()`] when matrices have fixed trace
#' + [`test_ss1()`] when the squared eigenvalues of each matrix sums to \eqn{1}.
#' + [`test_ss1fixedtrace()`] when the squared eigenvalues of each matrix sums to \eqn{1} and the trace is zero.
#' For single sample tests the null hypothesis is a user-provided set of eigenvalues for the population mean.
#' For \eqn{k}-sample tests the null hypothesis is that the eigenvalues of the population means are equal.
#'
#' Additionally [`test_unconstrained_aGOE()`] can perform \eqn{2}-sample tests using an approximation of the Gaussian Orthogonal Ensemble \insertCite{schwartzman2010gr}{TFORGE}.
#'
#' The above tests all require that the eigenvalues of the population mean are distinct (with decreasing performance as eigenvalues get closer to each other).
#' Use [`test_multiplicity()`] to test the eigenvalue-multiplicity of the population mean of a single sample (for unconstrained or fixed trace matrices).
#' A test of the same hypothesis that requires orthogonally-invariant covariance is also available through [`test_multiplicity_OI()`] \insertCite{schwartzman2008in}{TFORGE}.
#' 
#' In this package matrices within the same sample are considered independently and identically distributed.
#' Matrices are stored in their flattened form according to [`vech()`]. See [`fsm`] for details.
#' Samples may be provided as lists of matrices or in their flattened from so long as the column order matches that of [`vech()`].
#' Eigenvalues, when distinct, are assumed to be in descending order.

#' # Acknowledgements
#' Colleagues Andrew T. A. Wood and Janice Scealy played crucial roles in developing the statistical concepts and theory.
#'
#' The package includes `scel.R` for empirical likelihood by  \insertCite{owen:2013;textual}{TFORGE}, which is used for bootstrap resampling weights when the sum of squared eigenvalues is constrained.
#' 
#' This package on Ngunnawal and Ngambri Country. I thank the Country for its influence.
#' @references
#' \insertAllCited{}
"_PACKAGE"

