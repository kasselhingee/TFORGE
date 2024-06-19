# @title The sample mean, difference to the mean and covariance of elements for a set of matrices.
# @param ms Set of matrices as a list or a 3-array.
mmean <- function(ms){
  stopifnot(inherits(ms, "TFORGE_fsm"))
  invvech(colMeans(ms))
}

merr <- function(ms, mean = mmean(ms)){
  stopifnot(inherits(ms, "TFORGE_fsm"))
  mean = vech(mean)
  out <- t(t(ms) - mean)
  # class(out) <- c("TFORGE_fsm", class(out))
  return(out)
}

mcovar <- function(merr){
  stopifnot(inherits(merr, "TFORGE_fsm"))
  out <- cov(merr)
  indx <- which(lower.tri(invvech(merr[1, ]), diag = TRUE), arr.ind = TRUE)
  nam <- paste0("e", indx[, "row"], indx[, "col"])
  colnames(out) <- rownames(out) <- nam
  return(out)
}

msum <- function(ms){
  if (inherits(ms, "TFORGE_fsm")){return(invvech(colSums(ms)))}
  else{return(purrr::reduce(ms, `+`))}
}
