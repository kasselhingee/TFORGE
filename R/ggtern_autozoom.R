#' @title Function for computing ternary plot limits that zoom onto the data
#' @description Extracts data clumsily from built ggplot2 object. Would be smoother to use the data and aes inheritance, but that is beyond me at the moment.
#' @examples
#' Y <- matrix(runif(20, 0.3, 0.37), ncol = 2)
#' Y <- cbind(Y, 1-rowSums(Y))
#' a <- ggtern::ggtern(data.frame(Y)) +
#'  ggplot2::geom_point(ggplot2::aes(x = X1, y = X2, z = X3))
#' autozoom(a)
#' @export
autozoom <- function(plotobj){
  arr <- ggplot2::layer_data(plotobj, 1)[, c("x", "y", "z")]
  zoomlimit <- getternzoomlimits(arr)
  plotobj + 
    ggtern::tern_limit(L = zoomlimit[["L"]], T = zoomlimit[["T"]], R = zoomlimit[["R"]])
}

getternzoomlimits <- function(arr){
  colmaxs <- apply(arr, 2, max)
  names(colmaxs) <- c("L", "T", "R")
  colmaxs <- round(colmaxs, 2) #need supplied maxes to be roundish numbers for tern_limit to work
  colmins <- apply(arr, 2, min)
  names(colmins) <- c("L", "T", "R")
  x <- colmaxs
  A <- diag(-1, 3, 3) + 1
  Ainv <- solve(A)
  # tern_limit() solves the equation Ax = 1-colmaxs to get corresponding minimum limits
  # below keeps expanding maximums until minimums are sufficiently small
  while (any(x >= colmins)){
    incdir <- ((which(x >= colmins)[[1]]) %% 3) + 1
    colmaxs[incdir] <- colmaxs[incdir]+1E-2
    x <- round(Ainv %*% (1-colmaxs), 3)
  }
  return(colmaxs)
}

