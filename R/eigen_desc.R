#' @noRd
#' @title A wrapper around eigen() that sorts eigenvalues according to decreasing order
#' @description `base` `eigen()` orders eigenvalues by absolute size. This wrapper sorts the results in decreasing order, positive to negative.
#' @param ... Passed to [`base::eigen()`].
#' @return Same as [`base::eigen()`], but eigenvalues sorted in decreasing.
eigen_desc <- function(...){
  raw <- base::eigen(...)
  ord <- order(raw$values, decreasing = TRUE)
  raw$values <- raw$values[ord]
  if (!is.null(raw$vectors)){
    raw$vectors <- raw$vectors[, ord]
  }
  return(raw)
}

eigen <- function(...){
  warning("Using base::eigen() without sorting eigenvalues into descending order. Use eigen_desc() instead.")
  base::eigen(...)
}

