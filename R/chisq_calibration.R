#' Internal: Chi Squared Calibration for Testing
#' @description
#' Similar to `bootresampling()`, but uses a chi-squared distribution instead of an empirical distribution from bootstrapping.
#' @inheritParams bootresampling
#' @export
chisq_calib <- function(x, stat, df, ...){
  x <- as_flat(x)
  t0 <- stat(x, ...)
  pval <- 1-stats::pchisq(t0, df)
  
  out <- list(
    pval = pval,
    t0 = t0,
    df = df
  )
  class(out) <- c("TFORGE", class(out))
  return(out)
}
