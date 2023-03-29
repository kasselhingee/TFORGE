#' @title The sample mean, difference to the mean and covariance of elements for a set of matrices.
#' @param ms Set of matrices as a list or a 3-array.
#' @details Uses [`abind::abind()`] which can be very slow for large numbers of matrices as it tries to store them all in memory.
mmean <- function(ms){
  return(msum(ms)/length(ms))
}

merr <- function(ms, mean = mmean(ms)){
  if(is.array(ms)){if (length(dim(ms))==3){
    ms <- lapply(1:dim(ms)[3], function(i) ms[,,i])
  }}
  lapply(ms, function(m){m - mean})
}

mcovar <- function(merr){
  if(is.array(merr)){if (length(dim(merr))==3){
    merr <- lapply(1:dim(merr)[3], function(i) merr[,,i])
  }}
  mcovarls <- lapply(merr, function(m){
    x <- vech(m)
    x %*% t(x)
  })
  mmean(mcovarls)
}

msum <- function(ms){
  if(is.array(ms)){if (length(dim(ms))==3){
    ms <- lapply(1:dim(ms)[3], function(i) ms[,,i])
  }}
  return(purrr::reduce(ms, `+`))
}
