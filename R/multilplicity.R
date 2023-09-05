#' @name multiplicity
#' @title Tools for testing a specified multiplicity for a single sample
#' @param ms Sample of matrices
#' @param mult A vector giving the multiplicity of eigenvalues in descending order of eigenvalue size.
NULL

#' @describeIn multiplicity Bootstrap test
#' @param B The number of bootstrap samples
#' @export
test_multiplicity <- function(ms, mult, B){
  ms <- as.mstorsst(ms)
  ms_std <- standardise_multiplicity(ms, mult)
  res <- bootresampling(ms, ms_std, 
    stat = stat_multiplicity,
    B = B,
    mult = mult)
  return(res)
}

#' @describeIn multiplicity Test statistic
#' @export
stat_multiplicity <- function(ms, mult, NAonerror = FALSE){
  av <- mmean(ms)
  if (sum(mult) != ncol(av)){
    stop(paste("Sum of mult = ", mult, "is not equal to ", ncol(av), collapse = " "))
  }
  stopifnot(all(mult > 0))
  stopifnot(any(mult != 1))
  es <- eigen(av)

  #indices
  cmult <- cumsum(mult)
  esvalstart <- c(0, cmult[-length(cmult)]) + 1
  idxs <- lapply(1:length(mult), function(i){
    esvalstart[i] : cmult[i]
  })

  #prod w evecs - should be equal to evals since vectors calculated from av
  evals <- es$values
    warning("need to sort values in case of negative evals")
  #vapply(1:ncol(es$vectors), function(i){
  #  t(es$vectors[, i]) %*% av %*% es$vectors[, i]
  #}, FUN.VALUE = 1.1)

  # the random variables xi in sets per multiplicity because the weight matrix is different in each one
  xi <- xiget(evals, mult, idxs)

  C0 <- mcovar(merr(ms, mean = av)) # the covariance between elements of xi
  covar <- xicovar(mult, idxs, es$vectors, C0/length(ms))

  return(drop(t(xi) %*% solve_NAonerror(covar, NAonerror) %*% xi))
}

#this function turns the qYbarq into vectors perpendicular to 1,1,1,1,1,1...
evalsvs111 <- function(evals){
    if (length(evals) < 2){return(NULL)}
    wmat <- helmertsub(length(evals))
    return(wmat %*% evals)
}

# the random variables xi in sets per multiplicity because the weight matrix is different in each one
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

#' @describeIn multiplicity Standardise a sample to satisfy the (null) hypothesis of the given eigenvalue multiplicity in `mult`.
#' @export
standardise_multiplicity <- function(ms, mult){
  av <- mmean(ms)
  stopifnot(sum(mult) == ncol(av))
  stopifnot(all(mult > 0))
  stopifnot(any(mult != 1))
  es <- eigen(av, symmetric = TRUE)

  #indices
  cmult <- cumsum(mult)
  esvalstart <- c(0, cmult[-length(cmult)]) + 1


  # get eigenvalues for each multiplicity by averaging 
  evals <- vapply(1:length(mult), function(i){
    mean(es$values[esvalstart[i] : cmult[i]])
    warning("need to sort values in case of negative evals")
  }, FUN.VALUE = 1.1)
  
  #make matices with the correspond eigenvectors (with eval of 1 for now)
  mats <- lapply(1:length(mult), function(i){
    idx <- esvalstart[i] : cmult[i]
    msum(lapply(idx, function(k) es$vectors[, k] %*% t(es$vectors[, k])))
  })

  #combine evals and mats to create a matrix with the desired multiplicity of eigenvalues and the same eigenvectors of av
  newM <- msum(mapply('*', evals, mats, SIMPLIFY = FALSE))

  # make new sample out of newM and errors
  out <- lapply(merr(ms, mean = av), '+', newM)
  class(out) <- c(class(out), "sst")
  return(out)
}

