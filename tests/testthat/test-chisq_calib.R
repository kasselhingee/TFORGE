test_that("chisq_calib() is really just calling pchisq()", {
  set.seed(354)
  x <- rchisq(1, df = 4)
  res <- chisq_calib(rsymm_norm(3, diag(1, 3)), stat = function(y){x}, df = 4) #stat here overrides whatever x is for this test
  expect_equal(res$t0, x)
  expect_equal(res$pval, 1 - stats::pchisq(x, 4))
})