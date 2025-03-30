# @title The sample mean, difference to the mean and covariance of elements for fsm.
# @param ms Set of matrices as a list or a 3-array.

# mean
mmean <- function(ms){
  stopifnot(inherits(ms, "TFORGE_fsm"))
  invvech(colMeans(ms))
}

# difference to mean
merr <- function(ms, mean = mmean(ms)){
  stopifnot(inherits(ms, "TFORGE_fsm"))
  mean = vech(mean)
  out <- t(t(ms) - mean)
  # class(out) <- c("TFORGE_fsm", class(out))
  return(out)
}

# covariance
mcovar <- function(merr){
  stopifnot(inherits(merr, "TFORGE_fsm"))
  out <- stats::cov(merr)
  indx <- which(lower.tri(invvech(merr[1, ]), diag = TRUE), arr.ind = TRUE)
  nam <- paste0("e", indx[, "row"], indx[, "col"])
  colnames(out) <- rownames(out) <- nam
  return(out)
}

# sum
msum <- function(ms){
  if (inherits(ms, "TFORGE_fsm")){return(invvech(colSums(ms)))}
  else{return(purrr::reduce(ms, `+`))}
}
