#' @title The sample mean, difference to the mean and covariance of elements for a set of matrices.
#' @param ms Set of matrices as a list or a 3-array.
mmean <- function(ms){
  stopifnot(inherits(ms, "sst"))
  invvech(colMeans(ms))
}

merr <- function(ms, mean = mmean(ms)){
  stopifnot(inherits(ms, "sst"))
  mean = vech(mean)
  out <- t(t(ms) - mean)
  # class(out) <- c("sst", class(out))
  return(out)
}

mcovar <- function(merr){
  stopifnot(inherits(merr, "sst"))
  out <- cov(merr)
  indx <- which(lower.tri(merr[[1]], diag = TRUE), arr.ind = TRUE)
  nam <- paste0("e", indx[, "row"], indx[, "col"])
  colnames(out) <- rownames(out) <- nam
  return(out)
}

msum <- function(ms){
  if (inherits(ms, "sst")){return(invvech(colSums(ms)))}
  else{return(purrr::reduce(ms, `+`))}
}
