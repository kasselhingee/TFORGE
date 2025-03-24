#' @title Simulate Symmetric Matrices with Multivariate Normal Elements
#' @description
#' Simulate symmetric matrices with elements from a multivariate Normal distribution.
#' @param n Number of matrices to generate
#' @param mean A symmetric matrix specifying the mean of the distribution.
#' @param sigma 
#' A covariance matrix for the vectorized lower triangular elements (arranged by [`vech()`]) of the symmetric matrix. 
#' It is passed to [`mvtnorm::rmvnorm()`] without any transformation. 
#' @return A `TFORGE_fsm` object. See [`as_fsm()`].
#' @details 
#' The mean matrix is vectorized using the [`vech()`] function 
#' and then used as the mean vector in the [`mvtnorm::rmvnorm()`] function. The covariance matrix `sigma` is passed unchanged to [`mvtnorm::rmvnorm()`].
#' Symmetric matrices are obtained by applying [`invvech()`] to each simulated vector.
#' @examples 
#' rsymm_norm(100, diag(c(3,2,1)))
#' @export
rsymm_norm <- function(n, mean, sigma = diag(length(vech(mean)))){
  stopifnot(isSymmetric(mean))
  tmp <- mvtnorm::rmvnorm(n, mean = vech(mean), sigma = sigma)
  class(tmp) <- c("TFORGE_fsm", "array")
  return(tmp)
}
#' @rdname rsymm_norm
#' @export
rsymm <- rsymm_norm
