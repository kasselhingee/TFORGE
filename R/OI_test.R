#' @title Test OI Covariance
#' @description
#' Proposition 3.1 of Schwartzman et al (2008) provides a test of OI covariance with \eqn{\tau > 0} against unrestricted covariance structure.
#' The current method is not giving correct p-values. The statistic appears to be missing an offset of `q=p*(p+1)/2` to have the correct asymptotic distribution.
#' @details
#' The parameters \eqn{\tau}{tau} and \eqn{\sigma^2}{s^2} are estimated using [`estimate_OIcov()`] using the sample average as the estimate of population mean.
#'
#' Should \eqn{\tau < 0}, Schwartzman et al (2008) after the proof of proposition 3.1 says the asymptotic distribution is "guaranteed only if \eqn{\tau < 0}", 
#' which is opposite to the statement in proposition 3.1 itself.
#' @return A list of the p value, statistic, and estimated \eqn{\tau}{tau} and \eqn{\sigma^2}{s^2}.
#' @export
test_OIcov <- function(x){
  p <- as.integer((-1 + sqrt(8*ncol(x) + 1))/2)
  q <- ncol(x)
  if (nrow(x) <= q * (q+3)/2){warning(sprintf("Only %i samples, but more than %i=q(q+3)/2 is usually required.", nrow(x), q * (q+3)/2))}
  Ybar <- colMeans(x)
  OIparams <- estimate_OIcov(x, Ybar)
  if (OIparams$tau < 0){warning(sprintf("Estimated tau=%f is smaller than 0. The distribution of the test statistic is only valid when tau > 0.", OIparams$tau))}
  covhat <- S_mcovar(t(t(x)-Ybar))
  stat <- nrow(x) * (q * log(OIparams$scalesq) - log(1-p*OIparams$tau) - determinant(covhat, logarithm = TRUE)$modulus)
  stat <- as.numeric(stat)
  pval <- 1-pchisq(stat, df = q*(q+1)/2 - 2)
  return(list(
    pval = pval,
    stat = stat,
    scalesq = OIparams$scalesq,
    tau = OIparams$tau
  ))
}


