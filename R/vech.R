#' @title Stack elements of a symmetric matrix into a vector.
#' @description
#' The `vech` operator as defined by \insertCite{@Section 3.8 @magnus2019ma}{PivotalBootstrapMatrixData}.
#' @details
#' Columns are stacked in order from left to right, entries above the diagonal are ignored.
#' Matrix is not checked for symmetry.
#' `vech` uses [`lower.tri()`].
#' `invvech` is the opposite of `vech`. It could be much more efficient with memory and operations.
#' @param m A symmetric matrix
#' @param name If TRUE vector elements are named `eij` where `i` is the row ang `j` is the column.
#' @export
vech <- function(m, name = FALSE){
  out <- m[lower.tri(m, diag = TRUE)]
  if (name){
    indx <- which(lower.tri(m, diag = TRUE), arr.ind = TRUE)
    nam <- paste0("e", indx[, "row"], indx[, "col"])
    names(out) <- nam
  }
  return(out)
}

#' @describeIn vech The inverse of vech.
#' @export
invvech <- function(x,...){
  # n is s.t. l = 0.5 n (n+1) --> n = 0.5(-1 + sqrt(4l+1))
  n <- (-1 + sqrt(8*length(x) + 1))/2
  m <- matrix(NA, nrow = n, ncol = n)
  m[lower.tri(m, diag = TRUE)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

