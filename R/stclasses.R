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
  if (inherits(x, "list")){
    if (all(vapply(x, inherits, "sst", FUN.VALUE = FALSE))){ #sst objects already :)
      class(x) <- c(class(x), "mst")
      return(x)
      }
    if (all(vapply(x, inherits, "list", FUN.VALUE = FALSE))){ #a list of lists: possible multisample
      x <- lapply(x, as.sst, ...) #does nothing if elements already preprocessed
      dims <- do.call(rbind, lapply(x, function(y){dim(y[[1]])}))
      if (length(unique(dims[,2])) != 1){stop("Matrices in samples have different sizes.")}
      if (length(unique(dims[,2])) != 1){stop("Matrices in samples have different sizes.")}
      class(x) <- c(class(x), "mst")
      return(x)
    }
    # if not a list of lists, then assume a list of matrices (i.e. a single sample)
    if (all(vapply(x, inherits, "matrix", FUN.VALUE = FALSE))){return(as.sst(x, ...))}
  }
  stop("Could not convert to mst or sst.")
}

as.sst <- function(x, ...){
  if (inherits(x, "sst")){return(x)}
  if (inherits(x, "matrix")){stop("x is a single matrix")}
  if (!all(vapply(x, inherits, "matrix", FUN.VALUE = FALSE))){stop("Some elements are not matrices")}
  dims <- do.call(rbind, lapply(x, dim))
  if (dims[1,2] != dims[1,2]){stop("Matrices are not square.")}
  if (length(unique(dims[,2])) != 1){stop("Some matrices are different sizes.")}
  if (length(unique(dims[,2])) != 1){stop("Some matrices are different sizes.")}
  if (!all(vapply(x, isSymmetric, FUN.VALUE = FALSE, ...))){stop("Some matrices are not symmetric according to default limits in isSymmetric().")}
  xvec <- lapply(x, vech)
  xmat <- do.call(rbind, xvec)
  class(xmat) <- c("sst", class(xmat))
  return(xmat)
}

#' @export
`[.sst` <- function(x, i, j, ...){
  class(x) <- "matrix" #so it uses the default array indexing
  apply(x[i, , drop = FALSE], 1, invvech, simplify = FALSE)
}

#' @export
`[[.sst` <- function(x, i, j, ...){
  stopifnot(is.vector(i))
  class(x) <- "matrix" #so it uses the default array indexing
  if (length(i) == 2){
    return(invvech(x[i[[1]], ])[ i[[2]] ])
  }
  if (length(i) == 1){
    return(invvech(x[i, ]))
  }
}


