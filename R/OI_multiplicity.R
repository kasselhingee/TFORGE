#' @title Test of eigenvalue multiplicity assuming orthogonally invariant covariance
#' @description 
#' Given a sample from a population of symmetric matrices with Gaussian-distributed elements and orthogonally-invariant covariance, corollary 4.3 by \insertCite{schwartzman2008in;textual}{TFORGE} provides a method to test the eigenvalue multiplicity of the mean matrix.
#' Orthogonally-invariant covariance is a strong assumption and may not be valid; consider using [`test_multiplicity()`] if you are unsure.
#' @inheritParams test_multiplicity
#' @param refbasis Converts matrices to new basis.
#' @details 
#' The orthogonally invariant covariance matrix is estimated by [`estimate_OIcov()`]. The maximum-likelihood estimate of the population mean under the null hypothesis is computed according to \insertCite{@Theorem 4.2, @schwartzman2008in}{TFORGE}. 
#' @export
test_multiplicity_OI <- function(x, mult, B = "chisq", refbasis = diag(1, sum(mult))){
  # convert x to new basis
  x <- as_fsm(x)
  mats <- apply(x, 1, invvech, simplify = FALSE)
  mats <- lapply(mats, function(m){t(refbasis) %*% m})
  x <- as_fsm(mats)
  
  if (B == "chisq"){
    out <- chisq_calib(x = x, stat_multiplicity_OI, mult = mult, df = 0.5 * sum(mult * (mult + 1)) - length(mult))
    return(out)
  }
  
  x_std <- standardise_multiplicity(x, mult)
  res <- bootresampling(x, x_std, 
                        stat = stat_multiplicity_OI,
                        B = B,
                        mult = mult)
  return(res)
}

stat_multiplicity_OI <- function(x, mult){
  mn <- colMeans(x)
  es_mn <- eigen_desc(invvech(mn))
  Mhat <- es_mn$vectors %*% diag(blk(es_mn$values, mult)) %*% t(es_mn$vectors)
  OIparams <- estimate_OIcov(x, Mhat)
  stat <- (nrow(x)/OIparams$scalesq) * sum((es_mn$values - blk(es_mn$values, mult))^2)
  attr(stat, "scalesq") <- OIparams$scalesq
  attr(stat, "tau") <- OIparams$tau
  return(stat)
}

# block into multiplicities a set of ordered eigenvalues, following Definition 4.1 Schwartzman (2008) for blk
# @param mult A vector giving the multiplicity of eigenvalues in descending order of eigenvalue size.
# @returns A vector of new eigenvalues
blk <- function(evals, mult){
  #indices of multiplicity - same as stat_multiplicity()
  cmult <- cumsum(mult)
  stopifnot(sum(mult) == length(evals))
  esvalstart <- c(0, cmult[-length(cmult)]) + 1
  idxs <- lapply(1:length(mult), function(i){
    esvalstart[i] : cmult[i]
  })
  multmeans <- lapply(idxs, function(idx) {mean(evals[idx])})
  blkvals <- unlist(mapply(rep, multmeans, mult, SIMPLIFY = FALSE))
  blkvals
}
