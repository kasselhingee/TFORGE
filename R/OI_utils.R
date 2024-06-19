# File for Single Sample Schwartzman 2008 Methods

#' @title Estimate parameters of orthogonally invariant covariance
#' @description
#' Orthogonally-invariant covariance is a restrictive structure, but if it holds then a suite of tools is available \insertCite{schwartzman2008in}{TFORGE}.
#' Any orthogonally-invariant covariance can be specified by just two parameters \eqn{\tau}{tau} and \eqn{\sigma^2}{s^2}.
#' For a Gaussian-distributed elements, the parameters \eqn{\tau}{tau} and \eqn{\sigma^2}{s^2} can be estimated by maximum-likelihood if provided the data and the maximum-likelihood estimate of the population mean \insertCite{@Lemma 3.3, @schwartzman2008in}{TFORGE}.
#' @param Mhat A maximum-likelihood estimate of the population mean
#' @param x A single sample of symmetric matrices (passed to [`as_fsm()`]).
#' @param tau If supplied only \eqn{\sigma^2}{s^2} will be estimated.
#' @returns A named list of \eqn{\sigma^2}{s^2} and \eqn{\tau}
#' @details
#' # Orthoganally-Invariant Covariance
#' A symmetric random matrix \eqn{Y} with a Gaussian distribution has orthogonally-invariant covariance if and only if \eqn{Q Y Q^T} has the same distribution as \eqn{Y} for any orthogonal matrix \eqn{Q}.
#' 
#' Using the parameterisation of \eqn{\tau}{tau} and \eqn{\sigma^2}{s^2} by \insertCite{schwartzman2008in;textual}{TFORGE}:
#'  + the covariance of the off-diagonal elements of \eqn{Y} is \eqn{I\sigma^2/2}{Is^2/2} where \eqn{I} is the identity matrix of the correct size.
#'  + the covariance of the diagonal elements of \eqn{Y} is \eqn{\sigma^2 (I + 1 1^T \tau/(1-\tau p) )}{s^2 (I + 1 1^T \tau/(1-\tau p) )} where \eqn{p} is the number of columns of \eqn{Y} and \eqn{1} is the vector of ones.
#'  + the covariance between diagonal elements and non-diagonal elements is zero (i.e. they are independent).
#' @references
#' \insertAllCited{}
#' @export
estimate_OIcov <- function(x, Mhat, tau = NULL){
  x <- as_fsm(x)
  p <- as.integer((-1 + sqrt(8*ncol(x) + 1))/2)
  q <- p * (p+1)/2
  Ybar <- colMeans(x)
  Yerr <- t(t(x) - Ybar)
  # estimate tau
  if (is.null(tau)){
    numerator <- mean(OIinnerprod_fsm(Yerr, Yerr, 1, (p+1)/2)) +
      OIinnerprod_fsm(matrix(Ybar - vech(Mhat), nrow = 1), 
                                 matrix(Ybar - vech(Mhat), nrow = 1), 1, (p+1)/2)
    ondiag <- isondiag_vech(ncol(x))
    trYi2Ybar <- rowSums(Yerr[, ondiag])
    denominator <- (ncol(x) - 1) * (mean(trYi2Ybar^2) + sum((Ybar - vech(Mhat))[ondiag])^2)
    tau <- -numerator/denominator
  }
  
  # now to estimate scale
  normYbarMhat <- OIinnerprod_fsm(matrix(Ybar - vech(Mhat), nrow = 1),
                                  matrix(Ybar - vech(Mhat), nrow = 1),
                                  1, tau)
  sYerr2 <- mean(OIinnerprod_fsm(Yerr, Yerr, 1, tau))
  q <- p * (p+1)/2
  scalesq <- sYerr2/q + normYbarMhat/q
  return(list(scalesq = scalesq, tau = tau))
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
OIcov <- function(p, s, tau, vectorisor = "vech"){
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
OIinnerprod_fsm <- function(Avecs, Bvecs, s, tau){
  stopifnot(ncol(Avecs) == ncol(Bvecs))
  stopifnot(nrow(Avecs) == nrow(Bvecs))
  isdiag <- isondiag_vech(ncol(Avecs))
  elprods <- Avecs * Bvecs
  tmp <- rowSums(elprods) + rowSums(elprods[,!isdiag, drop = FALSE]) - 
    tau * rowSums(Avecs[, isdiag, drop = FALSE]) * rowSums(Bvecs[, isdiag, drop = FALSE])
  tmp/s^2
}



