# use matrix exponential to get observations with non-negative eigenvalues
# it isn't quite like entries are lognormally distributed
rsymm_lognorm <- function(n, meanlog, sigmalog = diag(length(vech(meanlog)))){
  stopifnot(isSymmetric(meanlog))
  logX <- rsymm_norm(n, mean = meanlog, sigma = sigmalog)
  X <- apply(logX, 1, function(v){vech(matrixexp(invvech(v)))}) %>% t()
  X <- as_fsm(X)
  return(X)
}

# matrix exponential
matrixexp <- function(m){
  es <- eigen_desc(m)
  es$vectors %*% diag(exp(es$values)) %*% t(es$vectors)
}