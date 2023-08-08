# functions for testing eigenvalues when the trace of the matrices are fixed

#' Test whether the supplied sample(s) have fixed trace
#' @param x Either a list of samples, each sample being a list of matrices, or a single sample as a list of matrices.
#' @param tolerance Tolerance on the relative difference, passed to `all.equal()`
#' @details Credit to Hadley for much of this function: [https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-numeric-vector]
#' @export
#' @return `TRUE` or `FALSE`
hasfixedtrace <- function(x, tolerance = sqrt(.Machine$double.eps)){
  x <- as.mstorsst(x)
  if (inherits(x, "mst")){x <- unlist(x, recursive = FALSE)}
  traces <- vapply(x, function(y){sum(diag(y))}, FUN.VALUE = 1.64)
  tracerange <- range(traces) / mean(traces)
  isTRUE(all.equal(tracerange[1], tracerange[2], tolerance = tolerance))
}

