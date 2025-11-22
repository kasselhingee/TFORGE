# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(TFORGE)

# from part of testthat:::on_cran (ignoring the interactive part)
# do a fast check unless NOT_CRAN is explicitly true
# see tests/testthat/helper-skips.R
if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))){
  Sys.setenv("FAST_CHECK" = "true")
}

test_check("TFORGE")
