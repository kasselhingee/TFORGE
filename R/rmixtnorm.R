#' @title Simulate a mixture of T and Normal
#' @param mean The mean of the distribution
#' @param tdelta The mean of the T component
#' @param n The number of (matrix samples)
#' @param cov The covariance for both the Normal and T components
#' @param df The degrees of freedom of the T component
#' @param w The chance of drawing a matrix from the T component (as opposed to the Normal component)
#' @export
rmixtnorm <- function(n, mean, w = 0.2, cov = diag(length(vech(mean))), tdelta = mean + 1, df = 10){
  distA <- runif(n) < w
  stopifnot(df > 2)
  if (sum(distA) > 0){YA <- rsymm_t(sum(distA), delta = tdelta, df = df, sigma = cov * (df - 2)/df)
  } else { YA <- NULL }
  # from law of total expectation the mean of the Gaussian component must be (delta - w * tdelta)/(1-w)
  if (sum(!distA) > 0){
    YB <- rsymm_norm(sum(!distA), mean = (delta - w * tdelta)/(1-w), sigma = cov)
  } else {YB <- NULL}
  c(YA, YB)
}