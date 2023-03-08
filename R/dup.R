#' @title Duplication matrix
#' @description Duplication matrix as defined by \insertCite{@Section 3.8 @magnus2019ma}{PivotalBootstrapMatrixData}
#' @details The entries of the matrix are inferred from the indexes using [`which()`], so creation of the matrix can be computationally intensive for large `n`.

dup <- function(n){
  x <- 1:(n*(n+1)/2)
  m <- invvech(x)
  y <- vec(m)
  rowmap <- match(y, x)
  out <- matrix(0, nrow = length(y), ncol = length(x))
  out[cbind(1:length(rowmap), rowmap)] <- 1
  return(out)
}

