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

descendingordererror <- function(d0){
  # good help on withRestarts and related here: http://adv-r.had.co.nz/beyond-exception-handling.html
  withRestarts(
    stop(structure(
      class = c("est_evals_not_descending", "error", "condition"),
      list(message = paste("Estimated common eigenvalues are not in descending order:", paste(d0, collapse = " ")),
           call = sys.call(-1))
    )),
    ignore = function() d0,
    sort = function() sort(d0, decreasing = TRUE, na.last = TRUE),
    use_NA = function() NA * d0,
    use_value = function(xx) xx
  )
}
