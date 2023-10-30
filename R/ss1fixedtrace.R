#' @title Statistic and test of 3-matrices with trace 0 and sum of sq eigenvalues equal to 1
#' @details Warning: the null distribution of `stat_ss1fixedtrace()` for multisamples does not appear to be chi-sq.
#' @param x Multiple samples of matrices, all with the same trace. Or a single sample of matrices. See [`as.mstorsst()`] for required structure.
#' @param evals If supplied the eigenvalues of the null hypothesis and `evals` must sum to the trace of the matrices. For the multisample statistic this should be `NULL` and is estimated within the function.
#' @export
stat_ss1fixedtrace <- function(x, evals = NULL){
  x <- as.mstorsst(x)
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  stopifnot(all(dim(x[[1]][[1]]) == c(3,3)))
  if (is.null(evals) && (length(x) == 1)){warning("evals must be supplied for a meaningful statistic since x is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){warning("evals supplied, returned statistic is not a statistic for common eigenvalues between groups")}
  
  # some building blocks
  A0 <- matrix(c( 0, 1,-1,
                 -1, 0, 1,
                  1,-1, 0), ncol = 3, nrow = 3, byrow = TRUE)
  av <- lapply(x, mmean)
  ess <- lapply(av, eigen_desc)
  naveval <- lapply(ess, function(a){a$values/sqrt(sum(a$values^2))})
  
  # The Omegas
  Omegas <- mapply(function(ms, av, ess){
    covar_unconstrained <- cov_evals_est(ms = ms, evecs = ess$vectors, av = av)
    projmat <- (diag(1, nrow(av)) - (ess$values %*% t(ess$values)/sum(ess$values^2)))/sqrt(sum(ess$values^2))
    projmat %*% covar_unconstrained %*% projmat
    },
    ms = x,
    av = av,
    ess = ess,
    SIMPLIFY = FALSE)
  
  # now for the eigenvalue for the null
  if (is.null(evals)){
    #estimate according to equation after (41)
    mats <- mapply(function(aveval, Omega){
      t(A0) %*% aveval %*% t(aveval) %*% A0 / drop(t(aveval) %*% t(A0) %*% Omega %*% A0 %*% aveval)
      },
      aveval = naveval,
      Omega = Omegas,
      SIMPLIFY = FALSE)
    mat <- purrr::reduce(mats, `+`)
    eigenmat <- eigen_desc(mat)
    idx <- max(which(eigenmat$values > 0)) 
    #zero eigenvalue corresponds to the 1,1,1 eigenvector so skip if the scalar product is anything like 1
    # should be easy thershold because all other eigenvectors are perpendicular
    if (abs(eigenmat$vectors[, idx] %*% rep(1, 3)/ sqrt(3)) > 0.8){idx <- idx - 1}
    d0 <- eigenmat$vectors[, idx]
    # make d0 DOT evalsav have as much positive sign as possible
    dotprds <- vapply(naveval, function(v){v %*% d0}, FUN.VALUE = 0.1)
    avsign <- mean(sign(dotprds))
    if (avsign < 0){
      d0 <- -1 * d0
    }
    if (!all(order(d0, decreasing = TRUE) == 1:length(d0))){
      d0 <- descendingordererror(d0)
    }
  } else {
    d0 <- sort(evals / sqrt(sum(evals^2)), decreasing = TRUE)
    if (!isTRUE(all.equal(sum(d0), sum(diag(x[[1]][[1]]))))){stop("Provided evals do not sum to trace of observations.")}
  }
  
  # now the statistic (40) for each sample:
  persamplestat <- mapply(function(d3, Omega, n){
    n * (t(d0) %*% t(A0) %*% d3)^2 / (t(d0) %*% t(A0) %*% Omega %*% A0 %*% d0)
  },
  d3 = naveval,
  Omega = Omegas,
  n = lapply(x, nrow),
  SIMPLIFY = FALSE
  )
  stat <- drop(purrr::reduce(persamplestat, `+`))
  attr(stat, "null_evals") <- drop(d0)
  return(stat)
}

#' @describeIn stat_ss1fixedtrace Bootstrap test using `stat_ss1fixedtrace()`.
#' @param maxit Passed to [`el.test()`]
#' @details The number of iterations in [`el.test()`] can have a big influence on the result, and it seems the mean with the best empirical likelihood can often be on the boundary of the convex hull of the data.
#' `maxit = 25` is too small. Perhaps `maxit = 1000`?
#' @export
test_ss1fixedtrace <- function(x, evals = NULL, B, maxit = 25, sc = TRUE){
  x <- as.mstorsst(x)
  if (inherits(x, "sst")){x <- as.mstorsst(list(x))}
  if (is.null(evals) && (length(x) == 1)){stop("evals must be supplied for a meaningful test since mss is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}
  
  t0 <- stat_ss1fixedtrace(x, evals = evals)
  d0 <- attr(t0, "null_evals")
  
  # means corresponding to NULL and d0
  nullmeans <- lapply(x, elnullmean, d0 = d0, getcbound = TRUE)
  
  # el weights
  wts <- mapply(opt_el.test, ms = x, mu = nullmeans, maxit = maxit, sc = sc, SIMPLIFY = FALSE)
  
  #check the weights
  lapply(wts, function(w){
    if (abs(length(w) - sum(w)) > 0.9){
      # above sees if weight sums to n (otherwise should sum to k < n being number of points in face). 
      warning(sprintf("Empirical likelihood weights sum to %0.1f, which suggests the mean is on a face of the convex hull.", sum(w)))
    }
  })
  if (any(vapply(wts, sum, FUN.VALUE = 1.3) < 0.5)){
    warning("Empirical likelihood finds mean is outside the convex hull of a sample.")
    # above sees if weight sums to n (otherwise should sum to k < n being number of points in face). Assume proposed mean is outside convex hull and with pval of zero, t0 of +infty
    return(list(
      pval = 0,
      t0 = t0,
      nullt = NA,
      stdx = wts,
      B = NA
    ))
  }
  
  res <- bootresampling(x, wts, 
                        stat = stat_ss1fixedtrace,
                        B = B,
                        evals = evals)
  return(res)
}
