#' @rdname fsm
#' @title Flat symmetric matrix classes for storing and checking symmetric matrices
#' @description
#' The `fsm` class, short for 'Flat Symmetric Matrices' is for a collection symmetric matrices with each matrix stored as a row according to [`vech()`].
#' The `fsm` class is itself a thin wrapper of the matrix class.
#' So, for example, `x[1, ]` will return the flattened-form of the first matrix in the collection, and `invvech(x[1,])` will be the first matrix in non-flat form.
#'
#' The `kfsm` class is for a collection of multiple `fsm`.
#' @examples
#' x <- list(list(matrix(c(1,2,3,2,4,5,3,5,6), 3), 
#'                matrix(c(2,3,4,3,5,6,4,6,7), 3)),
#'           list(matrix(c(0.1,0.2,0.3,0.2,0.4,0.5,0.3,0.5,0.6), 3), 
#'                matrix(c(0.2,0.3,0.4,0.3,0.5,0.6,0.4,0.6,0.7), 3)))
NULL

#' @param x An object. 
#' @param ... Passed to `isSymmetric()` for testing whether matrices are symmetric.
#' @details 
#' + If `x` is a list of list of equal-sized matrices then it returns
#' `x` with class `ms` added.
#' If x is a list of symmetric matrices then it will become an `TFORGE_fsm`.
#' In the rare case that `x` is a list, and each element is a matrix *vectorised* matrices such that each element of `x` is symmetric then `as.mstorsst()` will mistakenly treat each each element of `x` as a symmetric tensor and return an `TFORGE_fsm` object.
#' @return An object with class `TFORGE_kfsm` or `TFORGE_fsm`.
#' @export
as.mstorsst <- function(x, ...){
  if (inherits(x, "TFORGE_kfsm")){return(x)}  #isa() requires a match on all elements of the class attribute, so inherits() more suitable
  if (inherits(x, "TFORGE_fsm")){return(x)}
  if (inherits(x, "matrix")){return(as_fsm(x))}

  return(as_kfsm(x, ...))
}

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
#' @param ... is passed to [`isSymmetric()`] for checking symmetric matrices.
#' @export
as_fsm <- function(x, ...){
  if (inherits(x, "TFORGE_fsm")){return(x)}
  if (inherits(x, "matrix")){
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
  class(xmat) <- c("TFORGE_fsm", class(xmat))
  return(xmat)
}

#' @export
summary.TFORGE_fsm <- function(object, ...){
  p <- dimfromvech(object[1, , drop = TRUE])
  list(Size = c("Number of matrices" = nrow(object), "Unflattened matrix ncol=nrow" = p, "Flattened ncol" = ncol(object)),
       #sprintf("%i flattened symmetric %i x %i matrices", nrow(object), p, p),
       "Element Summary" = NextMethod())
}

