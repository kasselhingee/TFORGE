#' @title Simulate Symmetric Matrices
#' @description
#' Simulate symmetric matrices with elements iid from a Normal distribution centred on zero.
#' @param n Number of matrices
#' @param p Number of rows and columns of the matrices
#' @param s The standard deviation(s) of the elements.
#' @return A list of matrices.
#' @export
rsymm <- function(n, p, sd = 1){
  replicate(n, invvech(rnorm(p*(p+1)/2, sd = sd)), simplify = FALSE)
}
