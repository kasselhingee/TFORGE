#' @title Simulate Symmetric Matrices
#' @description
#' Simulate symmetric matrices with elements iid from a uniform distributions centred on zero.
#' @param n Number of matrices
#' @param p Number of rows and columns of the matrices
#' @param s Width of the uniform distribution.
#' @return A list of matrices.
rsymm <- function(n, p, s = 1){
  replicate(n, invvech(runif(p*(p+1)/2, min = -s, max = s)), simplify = FALSE)
}
