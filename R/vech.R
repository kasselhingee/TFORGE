#' @title Flatten a symmetric matrix into a vector.
#' @description
#' The `vech` operator as defined by \insertCite{@Section 3.8 @magnus2019ma}{TFORGE} flattens symmetric matrices into vectors. Columns are extracted from left to right with entries above the diagonal ignored.
#' The function `invvech()` is the inverse of `vech()`.
#' The dimension of the matrix can be obtained from its flattened form by `dimfromvech()`.
#' @details
#' The extraction is conveniently performed by `m[lower.tri(m), diag = TRUE]`.
#' The matrix `m` is not checked for symmetry.
#' @param m A symmetric matrix
#' @param name If TRUE vector elements are named `eij` where `i` is the row and `j` is the column.
#' @examples
#' m <- invvech(1:6)
#' dimfromvech(1:6)
#' vech(m)
#' @references \insertAllCited{}
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

#' @rdname vech
#' @param x A flattened symmetric matrix (as a vector).
#' @export
invvech <- function(x){
  n <- dimfromvech(x)
  m <- matrix(NA, nrow = n, ncol = n)
  m[lower.tri(m, diag = TRUE)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

#' @rdname vech
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
