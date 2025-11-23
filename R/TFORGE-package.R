#' @importFrom Rdpack reprompt
#' @details
#' All tests in this package are for hypotheses about the eigenvalues of the (extrinsic) mean of the sampled population(s).
#' Functions for conducting hypothesis tests start with `test_`.
#' When the matrices are constrained to a submanifold of the space of symmetric matrices we use the extrinsic mean, which projects the usual linear/Euclidean mean so that it satisfies the same constraints as the data.
#' Tests can be calibrated using either a chi-squared distribution or bootstrapping; in simulations bootstrapping was more reliable but slower.
#'
#' The following functions conduct single sample and \eqn{k}-sample tests of the eigenvalues of the (extrinsic) population mean:
#' + [`test_unconstrained()`]
#' + [`test_fixedtrace()`] when matrices have fixed trace
#' + [`test_ss1()`] when the squared eigenvalues of each matrix sums to \eqn{1}.
#' + [`test_ss1fixedtrace()`] when the squared eigenvalues of each matrix sums to \eqn{1} and the trace is zero.
#'
#' The single sample tests conducted by the above functions test the null hypothesis of a user-provided set of eigenvalues for the extrinsic population mean against the alternative hypothesis that the extrinsic population mean has different eigenvalues.
#' The \eqn{k}-sample tests conducted by the above functions test the null hypothesis that the extrinsic population means have the same eigenvalues.
#'
#' Additionally [`test_unconstrained_aGOE()`] can perform \eqn{2}-sample tests with calibration by a Gaussian Orthogonal Ensemble (GOE) approximation \insertCite{schwartzman2010gr}{TFORGE}. A bootstrapped calibration for this test is also available where the GOE approximation is used to stabilise the scale of the statistic.
#'
#' There are two functions [`conf_fixedtrace()`] and [`conf_ss1fixedtrace()`] for estimating confidence regions.
#'
#' The above tests all require that the eigenvalues of the population mean are distinct (with decreasing performance when eigenvalues are very close to each other).
#' Eigenvalues are assumed to be in descending order.
#' 
#' Use [`test_multiplicity()`] to test the eigenvalue-multiplicity of the population mean of a single sample (for unconstrained or fixed trace matrices).
#' A test of the same hypothesis that requires that matrix elements follow a multivariate Gaussian distribution with orthogonally-invariant covariance is also available through [`test_multiplicity_OI()`] \insertCite{schwartzman2008in}{TFORGE}.
#' Eigenvalues are assumed to be in descending order for these multiplicity tests.
#' 
#' In this package, matrices within the same sample are considered independently and identically distributed.
#' Matrices are stored in a flattened form as row-vectors according to [`vech()`] - see [`fsm`] for details.
#' Samples may be provided as lists of matrices or in their flattened from so long as the column order matches that of [`vech()`].

#' # Acknowledgements
#' Colleagues Andrew T. A. Wood and Janice Scealy played crucial roles in developing the statistical concepts and theory.
#' This package was written on Ngunnawal and Ngambri Country.
#'
#' The package includes `scel.R` for empirical likelihood by  \insertCite{owen:2013;textual}{TFORGE}, which is used to estimate optimal weights for weighted bootstrapping of samples of constrained matrices.
#' The `scel.R` file was released in 2014-2015 under the under BSD-3-Clause with 
#' copyright by Board of Trustees, Leland Stanford Junior University
#' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#' "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#' LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#' FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#' HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#' SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#' LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
#' USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
#' AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#' OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#' OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#' 
#' @references
#' \insertAllCited{}
"_PACKAGE"

