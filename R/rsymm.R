#' @title Simulate Symmetric Matrices
#' @description
#' Simulate symmetric matrices with elements from a multivariate Normal distribution.
#' @param n Number of matrices
#' @param mean Centre of disribution. Must be symmetric.
#' @param sigma The standard deviation(s) of the elements.
#' @return A list of matrices.
#' @export
rsymm <- function(n, mean, sigma = diag(length(vech(mean)))){
  stopifnot(isSymmetric(mean))
  tmp <- mvtnorm::rmvnorm(n, mean = vech(mean), sigma = sigma)
  return(apply(tmp, MARGIN = 1, invvech, simplify = FALSE))
}
