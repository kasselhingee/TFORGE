# function skips if FAST_CHECK is "true"
skip_if_fast_check <- function() {
  if (identical(Sys.getenv("FAST_CHECK"), "true")) {
    testthat::skip("Slow test skipped.")
  } else {
    invisible()
  }
}


