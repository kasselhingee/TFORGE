#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @noRd
//' @title Helper for cov_evals().
//' @description Cpp version of `cov_evals_inside()`
//' @return Single numeric value
// [[Rcpp::export]]
double cov_evals_inside_cpp(
  const arma::colvec & vj,
  const arma::colvec & vk,
  const arma::mat & dupmat,
  const arma::mat & mcov
  ) {  
  arma::mat tmp = vj * vk.t();
  arma::mat tmp2 = dupmat.t() * arma::kron(tmp, tmp) * dupmat * mcov;
  double out = arma::trace(tmp2);
  return out;
  }
