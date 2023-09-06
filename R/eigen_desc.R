# a wrapper around eigen that sorts eigenvalues according to descending order, positive to negative
eigen_desc <- function(x, ...){
  raw <- base::eigen(x, ...)
  ord <- order(raw$values, decreasing = TRUE)
  raw$values <- raw$values[ord]
  if (!is.null(raw$vectors)){
    raw$vectors <- raw$vectors[, ord]
  }
  return(raw)
}

eigen <- function(...){
  warning("Using base::eigen() without sorting eigenvalues into descending order")
  base::eigen(...)
}

