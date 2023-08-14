#' @title Methods for testing eigenvalues with sum of squares = 1
#' @param x Multiple samples of matrices, all with the same trace. Or a single sample of matrices. See [`as.mstorsst()`] for required structure.
stat_ss1 <- function(x, evals = NULL, evecs = NULL, NAonerror = FALSE){
  x <- as.mstorsst(x)
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  
  # means and eigenspaces used in multiple parts, so calculate first here:
  mns <- lapply(x, mmean)
  ess <- lapply(mns, eigen)
  evalsav <- lapply(ess, "[[", "values")
  
  # first get each Omega2j:
  covars_unconstrained <- mapply(cov_evals, ms = x, evecs = lapply(ess, "[[", "vectors"), av = mns, SIMPLIFY = FALSE)
  Deltas <- lapply(evalsav, function(d2){amaral2007Lemma1(d2/sqrt(sum(d2^2)))})
  Omega2s <- mapply(function(d2, covar, Delta){
    Delta %*% covar %*% t(Delta) / sum(d2^2)
  }, d2 <- evalsav, covar = covars_unconstrained, Delta = Deltas)
  
  
}

# the construction of matrix A in Lemma1 Amaral, G. J. A., Dryden, I. L., & Wood, A. T. A. (2007). Pivotal Bootstrap Methods for k-Sample Problems in Directional Statistics and Shape Analysis. Journal of the American Statistical Association, 102(478), 695–707. http://www.jstor.org/stable/27639898
 
amaral2007Lemma1 <- function(m){
  d <- length(m)
  finalelement <- m[d]
  # the following is a rearrangement to make result less computationally sensitive to size of final element.
  A1 <- diag(1, d-1) - m[-d] %*% Conj(t(m[-d]))/(1+abs(finalelement))
  if (finalelement < 0){
    A1 <- -A1
  }
  A <- cbind(A1, -m[-d])
  return(A)
}