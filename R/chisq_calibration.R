#' Chi Squared Calibration for Testing
#' @description
#' Similar to [`boot_calib()`], but uses chi-squared calibration instead of bootstrapping.
#' @param df Degrees of freedom of the chi-squared distribution
#' @inheritParams boot_calib
#' @return 
#' A list of
#'  + `pval` the `p`-value from the test
#'  + `t0` the statistic for the observations `x`
#'  + `df` The degrees of freedom of the chi-squared distribution
#' 
#' The returned object has class `TFORGE` (same as [`boot_calib()`]) for easy use of `print()`.
#' @export
chisq_calib <- function(x, stat, df, ...){
  x <- as_flat(x)
  t0 <- stat(x, ...)
  pval <- 1-stats::pchisq(t0, df)
  attributes(pval) <- NULL
  
  out <- list(
    pval = pval,
    t0 = t0,
    df = df
  )
  class(out) <- c("TFORGE", class(out))
  return(out)
}
