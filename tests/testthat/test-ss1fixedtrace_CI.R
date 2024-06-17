test_that("conf_ss1fixedtrace() contains population mean about 95% of the time", {
  set.seed(345)
  
  popmeanincr <- replicate(100, {
    Y <- rsymm_norm(50, diag(c(1/sqrt(2), 0, -1/sqrt(2))), sigma = diag(6) * 0.01)
    Y <- projtrace_fsm(Y) #this shifts the distribution if the trace from rsymm_norm isn't symmertic about zero
    Y <- normL2evals_sst(Y) #replace eigenvalues with normalised ones. This changes the distribution, but I think it is symmetric about the mean normalised eigenvalues - just like averages of directions.
    stopifnot(has_fixedtrace(Y, tolerance = 1E10 * sqrt(.Machine$double.eps)))
    stopifnot(has_ss1(Y))
    cr <- suppressWarnings(conf_ss1fixedtrace(Y, 0.05, B = 100))
    conf_ss1fixedtrace_inregion(c(1, 0, -1), cr)
  })
  expect_lt(abs(mean(popmeanincr) - 0.95),
            2*sd(popmeanincr)/sqrt(length(popmeanincr)))
})
