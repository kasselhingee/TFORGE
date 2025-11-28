#' @title Test eigenvalues when trace=0 and sum of square eigenvalues = 1
#' @description
#' For a single sample, test eigenvalues of the population mean.
#' For multiple samples, test for equality of the eigenvalues of the population means.
#' This function is for 3x3 symmetric matrices with trace of zero and sum of squared eigenvalues of one.
#' These constraints combine so that the space of possible sets of (ordered) eigenvalues is a circle.
#' The test statistic is calculated by `stat_ss1fixedtrace()`.
#' @details
#' Test hypotheses described below.
#'
#' The sum of squared eigenvalues constraint forces the set of eigenvalues to lie on a sphere and the trace constraint forces eigenvalues onto a plane.
#' Combined the constraints force eigenvalues onto a circle in 3D Euclidean space.
#' The test statistic accounts for these constraints by projecting eigenvalues onto a line tangential to this circle and orthogonal to the null-hypothesis eigenvalues.
#'
#' Weighted bootstrap calibration is used (see 'Weighted Bootstrapping' below).
#'
#' Eigenvalues must be distinct.

#' @inheritParams test_unconstrained
#' @inheritParams test_fixedtrace
#' @inheritSection test_unconstrained Hypotheses
#' @inheritSection test_fixedtrace Weighted Bootstrapping
#' @inherit test_unconstrained return
#' @references \insertAllCited{}
#' @export
test_ss1fixedtrace <- function(x, evals = NULL, B = 1000, maxit = 25){
  x <- as_flat(x)
  stopifnot(has_ss1(x))
  stopifnot(has_fixedtrace(x))
  if (inherits(x, "TFORGE_fsm")){x <- as_flat(list(x))}
  if (is.null(evals) && (length(x) == 1)){stop("evals must be supplied for a meaningful test since mss is a single sample")}
  if (!is.null(evals) && (length(x) > 1)){stop("evals cannot be supplied when testing common eigenvalues between groups")}
  # check that evals satisfy ss1 constraint
  if (!is.null(evals)){if (abs(sqrt(sum(evals^2)) - 1) > sqrt(.Machine$double.eps)){stop("Square of evals do not sum to 1")}}
  if (!is.null(evals)){if (!isTRUE(all.equal(sum(evals), sum(diag(invvech(x[[1]][1, ])))))){stop("Provided evals do not sum to trace of observations.")}}
  
  # chisq calibration quick exit
  if (B == "chisq"){
    df <- (dim_fsm_kfsm(x) - 2) * (length(x) - is.null(evals))
    return(chisq_calib(x, stat_ss1fixedtrace, df = df, evals = evals))
  }

  # statistic and eigenvalues of the null hypothesis  
  t0 <- stat_ss1fixedtrace(x, evals = evals)
  d0 <- attr(t0, "null_evals")
  
  # means corresponding to NULL and d0
  nullmeans <- lapply(x, elnullmean, d0 = d0, getcbound = TRUE)
  
  # find weights that maximise empirical likelihood
  wts <- mapply(opt_el.test, ms = x, mu = nullmeans, maxit = maxit, SIMPLIFY = FALSE)
  
  #check the weights, return early if not good
  if (!wtsokay(wts)){
    out <- list(
      pval = 0,
      t0 = t0,
      nullt = NA,
      stdx = wts,
      B = NA
    )
    class(out) <- c("TFORGE", class(out))
    return(out)
  }
  
  res <- boot_calib(x, wts, 
                        stat = stat_ss1fixedtrace,
                        B = B,
                        evals = evals)
  return(res)
}

#' @rdname test_ss1fixedtrace
#' @export
stat_ss1fixedtrace <- function(x, evals = NULL){
  x <- as_flat(x)
  if (inherits(x, "TFORGE_fsm")){x <- as_flat(list(x))}
  stopifnot(ncol(x[[1]]) == 6) #corresponds to 3x3 matrix
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
    covar_unconstrained <- cov_evals_est(x = ms, evecs = ess$vectors, av = av)
    projmat <- (diag(1, nrow(av)) - (ess$values %*% t(ess$values)/sum(ess$values^2)))/sqrt(sum(ess$values^2))
    projmat %*% covar_unconstrained %*% projmat
    },
    ms = x,
    av = av,
    ess = ess,
    SIMPLIFY = FALSE)
  
  # now for the eigenvalue for the null
  if (is.null(evals)){
    mats <- mapply(function(n, aveval, Omega){
      n * t(A0) %*% aveval %*% t(aveval) %*% A0 / drop(t(aveval) %*% t(A0) %*% Omega %*% A0 %*% aveval)
      },
      n = vapply(x, nrow, FUN.VALUE = 1),
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
    if (abs(sqrt(sum(evals^2)) - 1) > sqrt(.Machine$double.eps)){warning("Normalising eigenvalues supplied to stat so that the squared eigenvalues sum to 1.")}
    d0 <- sort(evals / sqrt(sum(evals^2)), decreasing = TRUE)
    if (!isTRUE(all.equal(sum(d0), sum(diag(invvech(x[[1]][1, ])))))){stop("Provided evals do not sum to trace of observations.")}
  }
  
  # now the statistic for each sample:
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

