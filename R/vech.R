#' @title Stack elements of a symmetric matrix into a vector.
#' @description
#' The `vech` operator as defined by \insertCite{@Section 3.8 @magnus2019ma}{TFORGE}.
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
  n <- dimfromvech(x)
  m <- matrix(NA, nrow = n, ncol = n)
  m[lower.tri(m, diag = TRUE)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

#' @describeIn vech Dimension of matrix represented as a vector by `vech()`
#' @export
dimfromvech <- function(x){
  # n is s.t. l = 0.5 n (n+1) --> n = 0.5(-1 + sqrt(4l+1))
  n <- (-1 + sqrt(8*length(x) + 1))/2
  stopifnot(round(n) == n)
  return(n)
}

#for a vector created by vech, get the elements that come from the diagonal
# if vec is a single number, treat it as the length of vec
isondiag_vech <- function(vec){
  if (length(vec) == 1){vec <- rep(1, vec)}
  mat <- invvech(vec)
  mat[] <- FALSE;
  diag(mat) <- TRUE
  diagels <- vech(mat) > 0.5
  return(diagels)
}
