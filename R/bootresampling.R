#' @title Internal: Bootstrap calibration for single-sample and k-sample tests
#' @description
#' Performs a bootstrap hypothesis test using the supplied statistic.
#' The `stdx` parameter is used to define an empirical distribution that satisfies the null hypothesis. 
#' @param x Symmetric matrix observations. Either a [`fsm`] or a [`kfsm`].
#' @param stdx Either a [`fsm`] or [`kfsm`] of matrices transformed to satisfy the null hypothesis OR sampling weights for each matrix in `x` for weighted bootstrap calibration (sampling weights should be optimised and also generate an empirical distribution that satisfies the null hypothesis).
#' @param stat Function to compute the statistic.
#' @param B The number of bootstrap samples to use.
#' @param ... Passed to `stat`
#' @details
#' The function `stat` is applied to `x` and all resamples, with the result returned in the `t0` and `nullt` elements of the returned object, respectively.
#' 
#' Errors in evaluating `stat` on resamples are recorded in the `nullt_messages` and lead to `NA` values for the statistic in the `nullt` element of the returned object.
#' 
#' The `p`-value is the fraction of non-`NA` resample statistic values that are greater than `stat` applied to `x`.
#' @return
#' A list of
#'  + `pval` the `p`-value from the test
#'  + `t0` the statistic for the observations `x`
#'  + `nullt` The statistic evaluated on the resamples
#'  + `stdx` The `stdx` passed into `bootresampling()`
#'  + `B` The number of resamples requested
#'  + `nullt_messages` Any error messages for the corresponding resample
#'
#' The returned object has a bespoke class `TFORGE` for easy use of `print()`.
#' @export
bootresampling <- function(x, stdx, stat, B,  ...){
  stopifnot(is_single_whole_number(B))
  x <- as_flat(x)
  t0 <- stat(x, ...)
  exargs <- list(...) #extra arguments
  if (inherits(x, "TFORGE_kfsm")){
    if (inherits(stdx[[1]], "numeric")){ #stdx is a vector of weights - sample from x using these weights
      nullt_l <- replicate(B, catch_do.call(stat, c(list(multisample(x, prob = stdx)), exargs)), simplify = FALSE)
    } else { #if stdx isn't a vector of weights, assume it is a standardised version of x - sample from stdx
      nullt_l <- replicate(B, catch_do.call(stat, c(list(multisample(stdx)), exargs)), simplify = FALSE)
    }
  } else if (inherits(x, "TFORGE_fsm")){ #repeat above, but for a single sample data set
    if (inherits(stdx, "numeric")){
      nullt_l <- replicate(B, catch_do.call(stat, c(list(sample_fsm(x, prob = stdx, replace = TRUE)), exargs)), simplify = FALSE)
    } else {
      nullt_l <- replicate(B, catch_do.call(stat, c(list(sample_fsm(stdx, replace = TRUE)), exargs)), simplify = FALSE)
    }
  }
 
  # Extract results and capture any failure messages 
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
  class(out) <- c("TFORGE", class(out))
  return(out)
}

#' @noRd
#' @title Resample a multisample `TFORGE_kfsm` object
#' @description equivalent of sample, but for multiple samples
#' @param x an `TFORGE_kfsm`
#' @param prob weights. If present, must have the same structure as `x`
multisample <- function(x, prob = NULL){
  if (is.null(prob)){
    out <- lapply(x, sample_fsm, replace = TRUE)
    class(out) <- c("TFORGE_kfsm", class(out))
    return(out)
  }
  else {
    stopifnot(length(prob) == length(x))
    stopifnot(all(vapply(x, nrow, 2) == vapply(prob, length, 2)))
    out <- mapply(sample_fsm, x, prob = prob, MoreArgs = list(replace = TRUE), SIMPLIFY = FALSE)
    class(out) <- c("TFORGE_kfsm", class(out))
    return(out)
  }
}

# an internal version of sample that automatically marks the result as an TFORGE_fsm
sample_fsm <- function(x, prob = NULL, replace = TRUE){
  stopifnot(inherits(x, "TFORGE_fsm"))
  idx <- sample.int(nrow(x), prob = prob, replace = replace)
  out <- x[idx, ]
  class(out) <- c("TFORGE_fsm", class(out))
  return(out)
}

is_single_whole_number <- function(x, tol = .Machine$double.eps^0.5) {
  if (length(as.vector(drop(x))) > 1){return(FALSE)}
  abs(x - round(x)) < tol
}

#' @export
print.TFORGE <- function(x, ...){
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
