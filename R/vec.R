#' @title Stack elements of a matrix into a vector.
#' @description
#' The `vec` operator as defined by \insertCite{@Section 2.4 @magnus2019ma}{TFORGE}.
#' Columns are stacked in order from left to right.
#' `vec` calls the [`as.vector()`] function, which does exactly as required.
#' `invvec` is the opposite of `vec` and is a wrapper of [`matrix()`] with `byrow = FALSE`.
vec <- function(m){
  as.vector(m)
}

invvec <- function(x,...){
  matrix(x, byrow = FALSE, ...)
}

