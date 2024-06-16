# File for Single Sample Schwartzman 2008 Methods

#' @title Estimate of Orthogonally Invariant Covariance Parameters
#' @description
#' Implements MLE for \eqn{\tau}{tau} and \eqn{\sigma^2}{s^2} given in Lemma 3.3 of Schwartzman et (2008). The MLE for the population mean is required.
#' @param Mhat An MLE for the population mean (typically under a hypothesis)
#' @param ms A single sample in [`sst()`] format.
#' @param tau If supplied only \eqn{\sigma^2}{s^2} will be estimated.
#' @returns A named list of \eqn{\sigma^2}{s^2} and \eqn{\tau}
#' @export
estimateOIparams <- function(ms, Mhat, tau = NULL){
  ms <- as_fsm(ms)
  p <- as.integer((-1 + sqrt(8*ncol(ms) + 1))/2)
  q <- p * (p+1)/2
  Ybar <- colMeans(ms)
  Yerr <- t(t(ms) - Ybar)
  # estimate tau
  if (is.null(tau)){
    numerator <- mean(OIinnerprod_sst(Yerr, Yerr, 1, (p+1)/2)) +
      OIinnerprod_sst(matrix(Ybar - vech(Mhat), nrow = 1), 
                                 matrix(Ybar - vech(Mhat), nrow = 1), 1, (p+1)/2)
    ondiag <- isondiag_vech(ncol(ms))
    trYi2Ybar <- rowSums(Yerr[, ondiag])
    denominator <- (ncol(ms) - 1) * (mean(trYi2Ybar^2) + sum((Ybar - vech(Mhat))[ondiag])^2)
    tau <- -numerator/denominator
  }
  
  # now to estimate scale
  normYbarMhat <- OIinnerprod_sst(matrix(Ybar - vech(Mhat), nrow = 1),
                                  matrix(Ybar - vech(Mhat), nrow = 1),
                                  1, tau)
  sYerr2 <- mean(OIinnerprod_sst(Yerr, Yerr, 1, tau))
  q <- p * (p+1)/2
  scalesq <- sYerr2/q + normYbarMhat/q
  return(list(scalesq = scalesq, tau = tau))
}

#' @title Test OI Covariance
#' @description
#' Proposition 3.1 of Schwartzman et al (2008) provides a test of OI covariance with \eqn{\tau > 0} against unrestricted covariance structure.
#' The current method is not giving correct p-values. The statistic appears to be missing an offset of `q=p*(p+1)/2` to have the correct asymptotic distribution.
#' @details
#' The parameters \eqn{\tau}{tau} and \eqn{\sigma^2}{s^2} are estimated using [`estimateOIparams()`] using the sample average as the estimate of population mean.
#'
#' Should \eqn{\tau < 0}, Schwartzman et al (2008) after the proof of proposition 3.1 says the asymptotic distribution is "guaranteed only if \eqn{\tau < 0}", 
#' which is opposite to the statement in proposition 3.1 itself.
#' @returns A list of the p value, statistic, and estimated \eqn{\tau}{tau} and \eqn{\sigma^2}{s^2}.
#' @export
testOIcov <- function(ms){
  p <- as.integer((-1 + sqrt(8*ncol(ms) + 1))/2)
  q <- ncol(ms)
  if (nrow(ms) <= q * (q+3)/2){warning(sprintf("Only %i samples, but more than %i=q(q+3)/2 is usually required.", nrow(ms), q * (q+3)/2))}
  Ybar <- colMeans(ms)
  OIparams <- estimateOIparams(ms, Ybar)
  if (OIparams$tau < 0){warning(sprintf("Estimated tau=%f is smaller than 0. The distribution of the test statistic is only valid when tau > 0.", OIparams$tau))}
  covhat <- S_mcovar(t(t(ms)-Ybar))
  stat <- nrow(ms) * (q * log(OIparams$scalesq) - log(1-p*OIparams$tau) - determinant(covhat, logarithm = TRUE)$modulus)
  stat <- as.numeric(stat)
  pval <- 1-pchisq(stat, df = q*(q+1)/2 - 2)
  return(list(
    pval = pval,
    stat = stat,
    scalesq = OIparams$scalesq,
    tau = OIparams$tau
  ))
}

#' @title Create an Orthogonally Invariant Covariance Matrix from Parameters
#' @details According to Schwartzman et al (2008), orthogonally invariant covariance is such that:
#' 1. the covariance between diagonal elements of the random matrix is \eqn{\sigma^2(I_p + c)}{s^2(I_p + c)}, where \eqn{I_p} is the identity matrix of size \eqn{p}, the number of rows of the random matrix,
#' and \eqn{c} is related to \eqn{\tau}{tau} by \eqn{\tau = c/(1+pc)}.
#' 2. the covariance between the offdiagonal elements is \eqn{\sigma^2/2 I_{q-p}}{s^2/2 I_{q-p}} where \eqn{q = p(p+1)/2} is the number of different elements allowed in a symmetric matrix with \eqn{p} rows.
#' 3. the diagonal elements are independent of the off-diagonal elements.
#' 
#' For positive definiteness, \eqn{\tau} must be smaller than \eqn{1/p}
#' @param vectorisor Either 'vech' or 'vecd'. The covariance matrix in Schwartzman et al (2008) represents a random (data) matrix as a vector using [`vecd()`], which puts a \eqn{\sqrt{2}}{`sqrt(2)`} weight on the off diagonals, and has a different ordering to the [`vech()`] used by Hingee, Scealy and Wood (see [`vech2vecd_mat()`]).
#' @param s The scale \eqn{\sigma}{s} for the covariance matrix.
#' @param tau The offset \eqn{\tau}{tau} for the covariance related to the diagonal elements of the random matrix.
#' @export
covOI <- function(p, s, tau, vectorisor = "vech"){
  if (tau >= 1/p){warning("tau larger than 1/p and 'covariance' will not be positive definite")}
  offsetc <- tau/(1-tau*p)
  diagcov <- diag(1, p) + offsetc
  offdiagcov <- diag(1, p*(p+1)/2 - p) #also p(p-1)/2
  fullcov <- s^2 * blockdiag(diagcov, offdiagcov)
  switch(vectorisor,
         "vecd" = return(fullcov),
         "vech" = {
           inv_vech2vecd <- solve(vech2vecd_mat(nrow(fullcov)))
           return(inv_vech2vecd %*% fullcov %*% t(inv_vech2vecd))
         }
  )
}

# Schwartzman 2008 inner product in equations (10) and (12). This is the fast looking equation 12.
OIinnerprod <- function(A, B, s, tau){
  stopifnot(isSymmetric(A))
  stopifnot(isSymmetric(B))
  (sum(A * B) - tau * sum(diag(A)) * sum(diag(B)))/s^2
}

# Avecs and Bvecs in vech representation, rowwise between each row of A to each row of B
OIinnerprod_sst <- function(Avecs, Bvecs, s, tau){
  stopifnot(ncol(Avecs) == ncol(Bvecs))
  stopifnot(nrow(Avecs) == nrow(Bvecs))
  isdiag <- isondiag_vech(ncol(Avecs))
  elprods <- Avecs * Bvecs
  tmp <- rowSums(elprods) + rowSums(elprods[,!isdiag, drop = FALSE]) - 
    tau * rowSums(Avecs[, isdiag, drop = FALSE]) * rowSums(Bvecs[, isdiag, drop = FALSE])
  tmp/s^2
}


# block into multiplicities a set of ordered eigenvalues, following Definition 4.1 Schwartzman (2008) for blk
#' @param mult A vector giving the multiplicity of eigenvalues in descending order of eigenvalue size.
#' @returns A vector of new eigenvalues
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

#' @title Test of multiplicity given in Corllary 4.3 Schwartzman et al 2008
#' @details An estimate of the scale of the covariance is required, which uses [`estimateOIparams()`] and the MLE under the null hypothesis given by Theorem 4.2 of Schwartzman 2008
#' @export
test_multiplicity_OI <- function(Ysample, mult){
  mn <- colMeans(Ysample)
  es_mn <- eigen_desc(invvech(mn))
  Mhat <- es_mn$vectors %*% diag(blk(es_mn$values, mult)) %*% t(es_mn$vectors)
  OIparams <- estimateOIparams(Ysample, Mhat)
  stat <- (nrow(Ysample)/OIparams$scalesq) * sum((es_mn$values - blk(es_mn$values, mult))^2)
  pval <- 1 - pchisq(stat, df =  0.5 * sum(mult * (mult + 1)) - length(mult))
  return(list(
    pval = pval,
    stat = stat,
    scalesq = OIparams$scalesq,
    tau = OIparams$tau
  ))
}
