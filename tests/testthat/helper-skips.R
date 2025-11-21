# function skips if FAST_CHECK is "true"
fast_check_on <- function() {
  if (identical(Sys.getenv("FAST_CHECK"), "true")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
skip_if_fast_check <- function() {
  if (fast_check_on()) {
    testthat::skip("Slow test skipped.")
  } else {
    invisible()
  }
}


