#' @title Simulate Symmetric Matrices with Multivariate t Elements
#' @description
#' Simulate symmetric matrices with elements from a multivariate t distribution.
#' @param n Number of matrices
#' @param delta Centre of distribution. Must be symmetric.
#' @param sigma The scale paramater matrix for the elements arranged by [`vech()`]. `sigma` is passed to [`mvtnorm::rmvt()`] without any transformation.
#' @param df Degrees of freedom
#' @return A list of matrices.
#' @details The `delta` is vectorised (with [`vech()`]) and then [`mvtnorm::rmvt()`] is called with `type = "shifted"`
#' @examples 
#' Ys <- rsymm_t(100, delta = matrix(1, nrow = 3, ncol = 3), df = 10, sigma = diag(c(3,2,1,1,1,1)))
#' @export
rsymm_t <- function(n, delta, df = 1, sigma = diag(length(vech(delta)))){
  stopifnot(isSymmetric(delta))
  tmp <- mvtnorm::rmvt(n, delta = vech(delta), df = df, sigma = sigma)
  return(apply(tmp, MARGIN = 1, invvech, simplify = FALSE))
}