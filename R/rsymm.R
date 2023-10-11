#' @title Simulate Symmetric Matrices with Multivariate Normal Elements
#' @description
#' Simulate symmetric matrices with elements from a multivariate Normal distribution.
#' @param n Number of matrices
#' @param mean Centre of distribution. Must be symmetric.
#' @param sigma The covariance of the matrix elements arranged by [`vech()`]. `sigma` is passed to [`mvtnorm::rmvnorm()`] without any transformation.
#' @return A list of matrices.
#' @details The `mean` is vectorised (with [`vech()`]) and then [`mvtnorm::rmvnorm()`] is called.
#' @examples 
#' rsymm(100, diag(c(3,2,1)))
#' @export
rsymm_norm <- function(n, mean, sigma = diag(length(vech(mean)))){
  stopifnot(isSymmetric(mean))
  tmp <- mvtnorm::rmvnorm(n, mean = vech(mean), sigma = sigma)
  class(tmp) <- c("sst", class(tmp))
  return(tmp)
}
#' @export
rsymm <- rsymm_norm
