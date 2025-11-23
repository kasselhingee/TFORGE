#' @title Simulate Symmetric Matrices with Multivariate t Elements
#' @description
#' Simulate symmetric matrices with elements from a multivariate t distribution.
#' @param sigma The scale parameter matrix for the elements arranged by [`vech()`]. `sigma` is passed to [`mvtnorm::rmvt()`] without any transformation.
#' @return A `TFORGE_fsm` object. See [`as_fsm()`].
#' @param n Number of matrices to generate.
#' @param mean A symmetric matrix specifying the mean of the distribution.
#' @param df Degrees of freedom for the t distribution.
#' 
#' @details 
#' The function uses Representation A in \insertCite{lin1972ch}{TFORGE} to simulate multivariate-t vectors.
#' The mean matrix is vectorised using the [`vech()`] function 
#' and then used as the mean vector in the [`mvtnorm::rmvt()`] function.
#' The scale parameter matrix `sigma` is passed unchanged to [`mvtnorm::rmvt()`].
#' The covariance of the resulting vectors is `sigma * df / (df - 2)`.
#' Symmetric matrices are obtained by applying [`invvech()`] to each simulated vector.
#' @examples 
#' rsymm_t(100, mean = matrix(1, nrow = 3, ncol = 3), df = 10, sigma = diag(c(3,2,1,1,1,1)))
#' 
#' @references \insertAllCited{}
#' @export
rsymm_t <- function(n, mean, df = 1, sigma = diag(length(vech(mean)))){
  stopifnot(inherits(mean, "matrix"))
  stopifnot(isSymmetric(mean))
  mean <- vech(mean, name = TRUE)
  zeros <- rep(0, length(mean))
  nuS2 <- stats::rchisq(n, df)
  S2 <- nuS2/df
  sim_l <- lapply(S2, function(x){
    mvtnorm::rmvnorm(n = 1, mean = zeros, sigma = sigma / x)
  })
  out <- do.call(rbind, sim_l)
  out <- t(t(out) + mean)
  colnames(out) <- names(mean)
  class(out) <- c("TFORGE_fsm", "array")
  return(out)
}
