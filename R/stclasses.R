#' The ms and ss classes

#' @param x An object. 
#' @param ... Passed to `isSymmetric()` for testing whether matrices are symmetric.
#' @details 
#' + If `x` is a list of list of equal-sized matrices then it returns
#' `x` with class `ms` added.
#' @return An object with class `mst` or `sst`.
#' @export
as.mstorsst <- function(x, ...){
  if (inherits(x, "mst")){return(x)}  #isa() requires a match on all elements of the class attribute, so inherits() more suitable
  if (inherits(x, "sst")){return(x)}
  if (inherits(x, "matrix")){return(as.sst(x))}
  if (inherits(x, "list")){
    x <- lapply(x, as.sst, ...) #does nothing if elements already preprocessed
    dims <- vapply(x, ncol, FUN.VALUE = 3)
    if (length(unique(dims)) != 1){stop("Matrices in samples have different sizes.")}
      class(x) <- c("mst", class(x))
      return(x)
  }
  stop("Could not convert to mst or sst.")
}

as.sst <- function(x, ...){
  if (inherits(x, "sst")){return(x)}
  if (inherits(x, "matrix")){
    # check that columns numbers make sense
    tryCatch(dimfromvech(x[1, ]),
             error = function(e){
      if (grepl("round", e$message)){stop("Number of columns of x don't correspond to a possible result of vech()")}
      })
    class(x) <- c("sst", class(x))
    return(x)}

  # 

  # now if its a list of matrices
  if (is.list(x)) {if (!all(vapply(x, inherits, "matrix", FUN.VALUE = FALSE))){stop("Some elements are not matrices")}}
  dims <- do.call(rbind, lapply(x, dim))
  if (dims[1,2] != dims[1,2]){stop("Matrices are not square.")}
  if (length(unique(dims[,1])) != 1){stop("Some matrices are different sizes.")}
  if (length(unique(dims[,2])) != 1){stop("Some matrices are different sizes.")}
  if (!all(vapply(x, isSymmetric, FUN.VALUE = FALSE, ...))){stop("Some matrices are not symmetric according to default limits in isSymmetric().")}
  xvec <- lapply(x, vech)
  xmat <- do.call(rbind, xvec)
  class(xmat) <- c("sst", class(xmat))
  return(xmat)
}

