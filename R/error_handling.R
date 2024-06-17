# wrapper around solve that returns a matrix of NA if couldn't solve
solve_error <- function(A){
  erroraction <- function(e){
    if (!grepl("singular", e$message)){stop(e)}
    stop(structure(
      class = c("matrixsingular", "error", "condition"),
      list(message = e$message,
           call = e$call)
    ))
  }
  out <- tryCatch(solve(A), error = erroraction)
  out
}

