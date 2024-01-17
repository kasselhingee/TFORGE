#' @title A function for helping perform single-sample and k-sample bootstrap tests
#' @param x Observations as a list of matices, or list of list of matrices.
#' @param stdx List of matrices standardized to satisfy the null OR sampling weights for each matrix in `x`, in the same structure as `x`.
#' @param stat Function to compute the statistic
#' @param B The number of bootstrap samples to use
#' @param ... Passed to `stat`
#' @param NAonerror Passed to `stat` if `stat` has a formal argument called `NAonerror`. If `FALSE` then bootstrap resamples that lead to matrix inversion errors will return a statistic value of `NA`. `NAonerror` is *not* passed to the call of `stat` applied to the original data `x`.
#' @export
bootresampling <- function(x, stdx, stat, B, NAonerror = TRUE, ...){
  stopifnot(is_single_whole_number(B))
  x <- as.mstorsst(x)
  t0 <- stat(x, ...)
  exargs <- list(...)
  if ("NAonerror" %in% names(formals(stat))){ exargs <- c(exargs, NAonerror = NAonerror) }
  if (inherits(x, "mst")){
    if (inherits(stdx[[1]][[1]], "numeric")){
      #stdx is weights for an mst because first element of first sample is not a matrix/array, but just a numeric
      nullt_l <- replicate(B, catch_do.call(stat, c(list(multisample(x, prob = stdx)), exargs)), simplify = FALSE)
    } else {
      nullt_l <- replicate(B, catch_do.call(stat, c(list(multisample(stdx)), exargs)), simplify = FALSE)
    }
  } else if (inherits(x, "sst")){
    if (inherits(stdx[[1]], "numeric")){
      #stdx is weights for an sst because first element sample is not a matrix/array, but just a numeric
      nullt_l <- replicate(B, catch_do.call(stat, c(list(samplesst(x, prob = stdx, replace = TRUE)), exargs)), simplify = FALSE)
    } else {
      nullt_l <- replicate(B, catch_do.call(stat, c(list(samplesst(stdx, replace = TRUE)), exargs)), simplify = FALSE)
    }
  }
  
  nullt <- simplify2array(nullt_l)
  messages <- replicate(length(nullt), vector(mode = "character", length = 0))
  if (any(is.na(nullt))){
    warning(sprintf("The statistic could not be calculated for %i bootstrap resamples. See the returned nullt_messages for more information.", sum(is.na(nullt))))
    messages[is.na(nullt)] <- vapply(nullt_l[is.na(nullt)], attr, "message", FUN.VALUE = "abcd")
  }
  
  pval <- mean(nullt > t0, na.rm = TRUE)
  out <- list(
    pval = pval,
    t0 = t0,
    nullt = nullt,
    stdx = stdx,
    B = B,
    nullt_messages = messages
  )
  class(out) <- c(class(out), "tensorboot")
  return(out)
}

#equivalent of sample, but for multiple samples
#' @title Resample a multisample `mst` object
#' @param x an `mst`
#' @param prob weights. If present, must have the same structure as `x`
#' @export
multisample <- function(x, prob = NULL){
  if (is.null(prob)){
    out <- lapply(x, samplesst, replace = TRUE)
    class(out) <- c(class(out), "mst")
    return(out)
  }
  else {
    stopifnot(length(prob) == length(x))
    stopifnot(all(vapply(x, length, 2) == vapply(prob, length, 2)))
    out <- mapply(samplesst, x, prob = prob, MoreArgs = list(replace = TRUE), SIMPLIFY = FALSE)
    class(out) <- c(class(out), "mst")
    return(out)
  }
}

# an internal version of sample that automatically marks the result as an sst
samplesst <- function(x, prob = NULL, replace = TRUE){
  stopifnot(inherits(x, "sst"))
  idx <- sample.int(nrow(x), prob = prob, replace = replace)
  class(x) <- "matrix"
  out <- x[idx, ]
  class(out) <- c("sst", class(out))
  return(out)
}

is_single_whole_number <- function(x, tol = .Machine$double.eps^0.5) {
  if (length(as.vector(drop(x))) > 1){return(FALSE)}
  abs(x - round(x)) < tol
}

#' @export
print.tensorboot <- function(x, ...){
  x <- x[c("pval", "t0")]
  class(x) <- "list"
  NextMethod("print")
}

# For catching specific errors in above bootresampling, and any new ones
catch_do.call <- function(stat, args){
  tryCatch(do.call(stat, args),
           est_evals_not_descending = function(e){
             out <- NA_real_
             attr(out, "message") <- e$message
    out
  },
  matrixsingular = function(e){
    out <- NA_real_
    attr(out, "message") <- e$message
    out
  },
  zerocovariance = function(e){
    out <- NA_real_
    attr(out, "message") <- e$message
    out
  },
  error = function(e){
    out <- NA_real_
    attr(out, "message") <- e$message
    out
  })
}
