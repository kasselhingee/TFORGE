#' @title Test of eigenvalue multiplicity assuming orthogonally invariant covariance
#' @description 
#' Given a sample from a population of symmetric matrices with Gaussian-distributed elements and orthogonally-invariant covariance, corollary 4.3 by \insertCite{schwartzman2008in;textual}{TFORGE} provides a method to test the eigenvalue multiplicity of the mean matrix.
#' Orthogonally-invariant covariance is a strong assumption and may not be valid; consider using [`test_multiplicity()`] if you are unsure.
#' @inheritParams test_multiplicity
#' @details 
#' The orthogonally invariant covariance matrix is estimated [`estimate_OIparams()`]. The maximum-likilihood estimate of the population mean under the null hypothesis is computed according to \insertCite{@Theorem 4.2, @schwartzman2008in}{TFORGE}. 
#' @export
test_multiplicity_OI <- function(Ysample, mult){
  mn <- colMeans(Ysample)
  es_mn <- eigen_desc(invvech(mn))
  Mhat <- es_mn$vectors %*% diag(blk(es_mn$values, mult)) %*% t(es_mn$vectors)
  OIparams <- estimate_OIcov(Ysample, Mhat)
  stat <- (nrow(Ysample)/OIparams$scalesq) * sum((es_mn$values - blk(es_mn$values, mult))^2)
  pval <- 1 - pchisq(stat, df =  0.5 * sum(mult * (mult + 1)) - length(mult))
  return(list(
    pval = pval,
    stat = stat,
    scalesq = OIparams$scalesq,
    tau = OIparams$tau
  ))
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
