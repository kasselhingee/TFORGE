#' @name specifiedmultiplicity
#' @title Tools for testing a specified multiplicity for a single sample
#' @param ms Sample of matrices
#' @param mult A vector giving the multiplicity of eigenvalues in descending order of eigenvalue size.
NULL

#' @describeIn specifiedmultiplicity Standardise a sample to satisfy the (null) hypothesis of the given eigenvalue multiplicity
standardise_specifiedmultiplicity <- function(ms, mult){
  av <- mmean(ms)
  stopifnot(sum(mult) == ncol(av))
  stopifnot(all(mult > 0))
  es <- eigen(av, symmetric = TRUE)

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
  return(lapply(merr(ms, mean = av), '+', newM))
}

