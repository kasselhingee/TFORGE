#' @title Simulate Symmetric Matrices with Multivariate t Elements
#' @description
#' Simulate symmetric matrices with elements from a multivariate t distribution.
#' @param n Number of matrices
#' @param delta Centre of distribution. Must be symmetric.
#' @param sigma The scale paramater matrix for the elements arranged by [`vech()`]. `sigma` is passed to [`mvtnorm::rmvt()`] without any transformation.
#' @param df Degrees of freedom
#' @return A list of matrices.
#' @details Using the use Lin 1972 representation A to simulate about mean of zero. Then shift all simulated values by delta. The covariance of the result is sigma *  df / (df - 2).
#' @examples 
#' Ys <- rsymm_t(100, delta = matrix(1, nrow = 3, ncol = 3), df = 10, sigma = diag(c(3,2,1,1,1,1)))
#' @export
rsymm_t <- function(n, delta, df = 1, sigma = diag(length(vech(delta)))){
  stopifnot(isSymmetric(delta))
  delta <- vech(delta, name = TRUE)
  zeros <- rep(0, length(delta))
  nuS2 <- rchisq(n, df)
  S2 <- nuS2/df
  sim_l <- lapply(S2, function(x){
    mvtnorm::rmvnorm(n = 1, mean = zeros, sigma = sigma / x)
  })
  out <- do.call(rbind, sim_l)
  out <- t(t(out) + delta)
  colnames(out) <- names(mean)
  class(out) <- c("TFORGE_fsm", class(out))
  return(out)
}