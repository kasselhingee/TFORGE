#' @title Flat storage of symmetric matrices
#' @name fsm
#' @description
#' The `TFORGE_fsm` class, short for 'Flat Symmetric Matrices', is for storing a collection of symmetric matrices with each matrix stored as a row according to [`vech()`].
#' The `TFORGE_fsm` class is itself a thin wrapper of the array class.
#' So, for example, `x[1, ]` will return the flattened-form of the first matrix in the collection, and `invvech(x[1,])` will be the first matrix in non-flat form.
#' The `TFORGE_kfsm` class is for a collection of multiple `TFORGE_fsm`.
#' The function `as_flat()` automatically converts data to either `TFORGE_kfsm` or `TFORGE_fsm`.

#' @examples
#' x <- list(list(matrix(c(1,2,3,2,4,5,3,5,6), 3), 
#'                matrix(c(2,3,4,3,5,6,4,6,7), 3)),
#'           list(matrix(c(0.1,0.2,0.3,0.2,0.4,0.5,0.3,0.5,0.6), 3), 
#'                matrix(c(0.2,0.3,0.4,0.3,0.5,0.6,0.4,0.6,0.7), 3)))
#' as_kfsm(x)
#' summary(as_flat(x))
NULL

#' @describeIn fsm Automatically convert to either a `TFORGE_kfsm` or a `TFORGE_fsm`.
#' @details
#' The matrices inside `x` must all have the same dimension.
#' 
#' The function `as_flat()` automatically chooses between a `TFORGE_kfsm` or a `TFORGE_fsm`:
#' + If `x` is a list of symmetric matrices then it will return a `TFORGE_fsm`.
#' + If `x` is a list of lists of equal-sized matrices then it returns a `TFORGE_kfsm`, with each element of the larger list a `TFORGE_fsm`.
#' + If `x` is a list of 2D arrays, each satisfying `as_fsm()`, then `as_flat()` will return a `TFORGE_kfsm`.
#' + In the rare case that `x` is a list of 2D arrays of flattened matrices, but the 2D arrays happen to be perfectly symmetric (requires size of collections to perfectly relate to the dimension of the matrix observations) then `as_flat()` will mistakenly treat each element of `x` as a symmetric matrix and return a `TFORGE_fsm`.
#' @param x For `as_fsm()` a list of symmetric matrices, or a 2D array of flattened matrices. For `as_kfsm()` a list of objects that can be passed to `as_fsm()` (i.e. a list of a lists of matrices, or a list of 2D arrays). For `as_flat()`, `x` can be suitable for either `as_fsm()` or `as_kfsm()`.
#' @param ... Passed to [`isSymmetric()`] for testing whether matrices are symmetric.
#' @return An object with class `TFORGE_kfsm` or `TFORGE_fsm`.
#' @export
as_flat <- function(x, ...){
  if (inherits(x, "TFORGE_kfsm")){return(x)}  #isa() requires a match on all elements of the class attribute, so inherits() more suitable
  if (inherits(x, "TFORGE_fsm")){return(x)}
  if (inherits(x, "array")){return(as_fsm(x))}

  return(as_kfsm(x, ...))
}

#' @describeIn fsm Convert multiple collections of matrices into a `kfsm`. `x` must be a list, with each entry of `x` a separate collection of matrices.
#' @export
as_kfsm <- function(x, ...){
  if (inherits(x, "TFORGE_kfsm")){return(x)}  #isa() requires a match on all elements of the class attribute, so inherits() more suitable
  if (!inherits(x, "list")){stop("x must be a list")}

  # try converting each entry to an fsm
  val <- try(x <- as_fsm(x, ...), silent = TRUE)
  if (!inherits(val, "try-error")){return(x)}
  else { #if the TFORGE_fsm fails then...
    x <- lapply(x, as_fsm, ...) #does nothing if elements already preprocessed
    dims <- vapply(x, ncol, FUN.VALUE = 3)
    if (length(unique(dims)) != 1){stop("Matrices in samples have different sizes.")}
    class(x) <- c("TFORGE_kfsm", class(x))
    return(x)
  }
}

#' @describeIn fsm For `x` a list of symmetric matrices of the same size, flattens `x` into a 2D array the `i`th row is a flattened version `vech(x[[i]])` of the `i`th matrix of  `x`. If `x` is already flattened then `as_fsm()` will check that the number of columns are consistent with a flattened symmetric matrix.
#' @export
as_fsm <- function(x, ...){
  if (inherits(x, "TFORGE_fsm")){return(x)}
  if (inherits(x, "array")){
    # check that columns numbers make sense
    tryCatch(dimfromvech(x[1, ]),
             error = function(e){
      if (grepl("round", e$message)){stop("Number of columns of x don't correspond to a possible result of vech()")}
      })
    class(x) <- c("TFORGE_fsm", class(x))
    return(x)
  }

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
  class(xmat) <- c("TFORGE_fsm", setdiff(class(xmat), "matrix"))
  return(xmat)
}

#' @export
summary.TFORGE_fsm <- function(object, ...){
  p <- dim_fsm_kfsm(object)
  list(Size = c("Number of matrices" = nrow(object), "Unflattened matrix ncol=nrow" = p, "Flattened ncol" = ncol(object)),
       #sprintf("%i flattened symmetric %i x %i matrices", nrow(object), p, p),
       "Element Summary" = NextMethod())
}

#' @export
summary.TFORGE_kfsm <- function(object, ...){
  p <- dim_fsm_kfsm(object)
  list(c("Number of samples" = length(object),
    "Unflattened matrix ncol=nrow" = p,
    "Flattened ncol" = ncol(object[[1]])),
    "Sample sizes" = vapply(object, nrow, FUN.VALUE = 1))
}

dim_fsm_kfsm <- function(object){
  if (inherits(object, "TFORGE_fsm")){return(dimfromvech(object[1, , drop = TRUE]))}
  if (inherits(object, "TFORGE_kfsm")){return(dimfromvech(object[[1]][1, , drop = TRUE]))}
  else {stop("object needs to have class TFORGE_fsm or TFORGE_kfsm")}
}
