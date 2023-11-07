# File for Single Sample Schwartzman 2008 Methods

#' Estimate of Orthogonally Invariant Covariance Parameters
#' @description
#' Implements MLE for \eqn{\tao}{tau} and \eqn{\sigma^2}{s^2} given in Lemma 3.3 of Schwartzman et (2008). The MLE for the population mean is required.
NULL

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
covOI <- function(p, s, tau, vectorisor = "vech"){
  if (tau >= 1/p){warning("tau larger than 1/p and 'covariance' will not be positive definite")}
  offsetc <- tau/(1-tau*p)
  diagcov <- diag(1, p) + offsetc
  offdiagcov <- diag(1, p*(p+1)/2 - p) #also p(p-1)/2
  fullcov <- s^2 * blockdiag(diagcov, offdiagcov)
  switch(vectorisor,
         "vecd" = return(fullcov),
         "vech" = return(t(vech2vecd_mat(nrow(fullcov))) %*% fullcov %*% vech2vecd_mat(nrow(fullcov)))
  )
}

