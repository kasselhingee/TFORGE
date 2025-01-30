#' @noRd
#' @title Duplication matrix
#' @description Duplication matrix as defined by \insertCite{@Section 3.8 @magnus2019ma}{TFORGE}
#' @details The entries of the matrix are inferred from the indexes using [`which()`], so creation of the matrix can be computationally intensive for large `n`.
#' 
#' The function is memoised so calling dup on the same n does not take any more time.


dup_direct <- function(n){
  x <- 1:(n*(n+1)/2)
  m <- invvech(x)
  y <- vec(m)
  rowmap <- match(y, x)
  out <- matrix(0, nrow = length(y), ncol = length(x))
  out[cbind(1:length(rowmap), rowmap)] <- 1
  return(out)
}

# function f must have atomic inputs only
memoise_simple_function <- function(f) {
  cache <- new.env(parent = emptyenv())  # Create an isolated environment
  
  function(...) {
    key <- paste(..., collapse = "_")  # Simple key (works for atomic inputs)
    if (exists(key, envir = cache)) {
      return(get(key, envir = cache))  # Return cached result
    }
    
    result <- f(...)  # Compute function
    assign(key, result, envir = cache)  # Store result in cache
    return(result)
  }
}

dup <- memoise_simple_function(dup_direct)