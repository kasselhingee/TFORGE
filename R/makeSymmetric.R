# remove tiny floating point differences to make a matrix exactly symmetric
makeSymmetric <- function(m, tolerance = sqrt(.Machine$double.eps) * 1E2){
  offdiag1 <- m[lower.tri(m)]
  offdiag2 <- t(m)[lower.tri(m)]
  if (!isTRUE(all.equal(offdiag1, offdiag2, tolerance = tolerance))){
    stop("Matrix m is not close to symmetric.")
  }
  m[lower.tri(m)] <- (offdiag1 + offdiag2)/2
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}
