#' @title Test eigenvalue multiplicity
#' @description
#' Given a sample from a population, uses bootstrap resampling to test multiplicity hypotheses of the population mean's eigenvalues.
#' The test statistic is computed by `stat_multiplicity()`, which includes a uniformly random rotation of eigenvectors associated with each eigenvalue of multiplicity greater than 1.
#' The null hypothesis is the multiplicity of eigenvalues of the population mean.
#' Bootstrap resampling is conducted from the null hypothesis, which uses the original sample converted to satisfy the null hypothesis by `standardise_multiplicity()`.
#' @details
#' This hypothesis test works on unconstrained symmetric matrices or matrices constrained to have fixed trace. It may work poorly on matrices with other constraints, on samples smaller than 15, or samples of multimodal populations.
#' 
#' Due to the random rotation of the eigenvectors, use [`set.seed()`] before `stat_multiplicity()` or `test_multiplicity()` if you want the answer to be repeatable.
#' @param x A single sample of matrices (passed to [`as_fsm()`]).
#' @param mult A vector specifying the eigenvalue multiplicity under the null hypothesis in descending order of eigenvalue size.
#' @param B The number of bootstrap samples
#' @examples
#' x <- rsymm_norm(15, mean = diag(c(2, 1, 1, 0)))
#' test_multiplicity(x, mult = c(1, 2, 1))
#' @return
#' + `test_multiplicity()` returns a `TFORGE` object including the p-value of the test (slot `pval`) and the statistic for `x` (slot `t0`). See [`bootresampling()`].
#' @export
test_multiplicity <- function(x, mult, B = 1000){
  x <- as_flat(x)
  ms_std <- standardise_multiplicity(x, mult)
  res <- bootresampling(x, ms_std, 
    stat = stat_multiplicity,
    B = B,
    mult = mult)
  return(res)
}

#' @rdname test_multiplicity
#' @param evecs For debugging only. Supply eigenvectors of population mean.
#' @return
#' + `stat_multiplicity()` returns a single value.
#' @export
stat_multiplicity <- function(x, mult, evecs = NULL){
  av <- mmean(x)
  if (sum(mult) != ncol(av)){
    stop(paste("Sum of mult = ", mult, "is not equal to ", ncol(av), collapse = " "))
  }
  stopifnot(all(mult > 0))
  stopifnot(any(mult != 1))
  if (!is.null(evecs)){
    warning("evecs should only be supplied for debugging")
    es <- list()
    es$vectors <- evecs
  } else {
    es <- eigen_desc(av)
  }

  #indices of multiplicity
  cmult <- cumsum(mult)
  esvalstart <- c(0, cmult[-length(cmult)]) + 1
  idxs <- lapply(1:length(mult), function(i){
    esvalstart[i] : cmult[i]
  })
  # Get evals after random rotations of the eigenvectors for each eigenvalue
  # To uniformly randomly rotate them in each block, need uniform rotation matrices
  # these can be found by uniform simulation within the basis of the given eigenvectors
  # I can represent the existing basis as e1,e2,e3 etc, and create a new basis by using the random matrix to project e1,e2,e3 onto new directions of e1,e2,e3.
  es$vectors <- rotevecs(es$vectors, idxs)
  es$values <- diag(t(es$vectors) %*% av %*% es$vectors)
  
  # the random variables xi in sets per multiplicity because the weight matrix is different in each one
  xi <- xiget(es$values, mult, idxs)

  C0 <- mcovar(merr(x, mean = av)) # the covariance between elements of xi
  covar <- xicovar(mult, idxs, es$vectors, C0/nrow(x))

  return(drop(t(xi) %*% solve_error(covar) %*% xi))
}

#this function turns the qYbarq into vectors perpendicular to 1,1,1,1,1,1...
evalsvs111 <- function(evals){
    if (length(evals) < 2){return(NULL)}
    wmat <- helmertsub(length(evals))
    return(wmat %*% evals)
}

# the random variables xi in sets per multiplicity because the weight matrix is different in each one
# @param idxs gives the locations in the full eigenvalue vector of the eigenvalues corresponding to each multiplicity
xiget <- function(evals, mult, idxs){
  xi <- lapply(1:length(mult), function(j){
    evals <- evals[ idxs[[j]] ]
    return(evalsvs111(evals))
  })
  xi <- unlist(xi)
  return(xi)
}

# compute/estimate covariance of xi, Cav is the covariance of av = C0/n
xicovar <- function(mult, idxs, evecs, Cav){
  # for each pair of distinct eigenvalues with multiplicity>1, create the matrix of covariance between xi, these are entries in the block-matrix form of V3 in my notes 
  valpairs <- expand.grid(1:length(idxs), 1:length(idxs))
  blocks <- apply(valpairs, MARGIN = 1, function(jk){
    if (mult[jk[1]] < 1.5){return(NULL)}
    if (mult[jk[2]] < 1.5){return(NULL)}
    helmertsub(mult[jk[1]]) %*% covarbetweenevals(jk[1], jk[2], idxs, evecs, Cav) %*% t(helmertsub(mult[jk[2]]))
  }, simplify = FALSE)

  #bind the blocks together appropriately
  # bind by row, columns still seperate
  cblocks <- lapply(1:length(mult), function(k){
     purrr::reduce(blocks[valpairs[, 2] == k], `rbind`)
  })
  covar <- purrr::reduce(cblocks, `cbind`)
  return(covar)
}

# for each pair of distinct eigenvalues with multiplicity>1, create the matrix of covariance between eigenvalues (I've called this matrix A_jk in my notes)
# j and k refer to index of *distinct* eigenvalues, i.e. entries of mult
# idxs is a list of vectors, each vector contains the columns of evecs that have columns corresponding to the vectors eigenvalue
covarbetweenevals <- function(j, k, idxs, evecs, Cav){
  Dp <- dup(nrow(evecs))
  idxj = idxs[[j]]
  idxk = idxs[[k]]
  if ((length(idxj) == 1 ) || (length(idxk) == 1)){return(NULL)}
  vecpairs <- expand.grid(idxj, idxk)
  At <- apply(vecpairs, MARGIN = 1, function(uv){
    # qjuqkv <- evecs[, ] %*% t(evecs[, uv[2]])
    qjuqju <- evecs[, uv[1]] %*% t(evecs[, uv[1]])
    qkvqkv <- evecs[, uv[2]] %*% t(evecs[, uv[2]])
    t(vec(qjuqju)) %*% Dp %*% Cav %*% t(Dp) %*% vec(t(qkvqkv))
  }, simplify = TRUE)
  dim(At) <- c(length(idxj), length(idxk)) # converts to matrix taking advantage that the first index is filled fastest
  return(At) # this is actually the transpose of whats in my notes (where rows correspond to idxk)
}

#' @rdname test_multiplicity
#' @return
#' + `standardise_multiplicity()` returns a sample of matrices stored in flattened form (a `TFORGE_fsm`).
#' @export
standardise_multiplicity <- function(x, mult){
  x <- as_fsm(x)
  av <- mmean(x)
  stopifnot(sum(mult) == ncol(av))
  stopifnot(all(mult > 0))
  stopifnot(any(mult != 1))
  es <- eigen_desc(av, symmetric = TRUE)

  #indices
  cmult <- cumsum(mult)
  esvalstart <- c(0, cmult[-length(cmult)]) + 1


  # get eigenvalues for each multiplicity by averaging 
  evals <- vapply(1:length(mult), function(i){
    mean(es$values[esvalstart[i] : cmult[i]])
  }, FUN.VALUE = 1.1)
  
  #make matices with the correspond eigenvectors (with eval of 1 for now)
  mats <- lapply(1:length(mult), function(i){
    idx <- esvalstart[i] : cmult[i]
    msum(lapply(idx, function(k) es$vectors[, k] %*% t(es$vectors[, k])))
  })

  #combine evals and mats to create a matrix with the desired multiplicity of eigenvalues and the same eigenvectors of av
  newM <- msum(mapply('*', evals, mats, SIMPLIFY = FALSE))

  # make new sample out of newM and errors
  newM <- vech(newM)
  av <- vech(av)
  out <- t(t(x) - av + newM)
  return(out)
}

# random matrices in the uniform distribution on the Stiefel manifold can be obtained as the L in X = TL decomposition of a random Normal matrix (elements iid normal), where
# T is lower triangular (and square) and L is unitary (not necessarily square).
#[ P67-68 Gupta and Nagar 1999]
# This is not quite the QR decomposition in LINPACK which is X = QR such that Q is orthonormal (square matrix) and R is upper triangular.
# t(X) = QR => X = R^T Q^T.
# R^T is lower triangular ==> T
# Q^T is orthonormal ==> L.
# So when X if square the Q of the QR decomposition of t(X) is the transpose of a uniformly random matrix on the Stiefel manifold.
# @param p the nrows and columns
runifortho <- function(p){
  # rstiefel::rustiefel(p,p)
  m <- matrix(stats::rnorm(p * p),
         nrow = p,
         ncol = p)
  out <- qr.Q(qr(t(m), LAPACK = TRUE))
  return(t(out))
}

# randomly rotate eigenvectors for each eigenspace
rotevecs <- function(evecs, idxs){
  rotevecs <- lapply(idxs, function(idxforeval){
    if (length(idxforeval) == 1){return(evecs[, idxforeval, drop = FALSE])}
    rot <- runifortho(length(idxforeval))
    t(rot %*% t(evecs[, idxforeval]))
  })
  do.call(cbind, rotevecs)
}
