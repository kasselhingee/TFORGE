test_that("stat_specifiedevals() doesn't reject for simulation of single sample from null and rejects otherwise", {
  set.seed(13131)
  Ysample <- rsymm(50, 3)
  Ysample <- lapply(Ysample, `+`, diag(c(3,2,1)))
  res <- test_specifiedevals(Ysample, c(3,2,1), 100)
  expect_gt(res$pval, 0.2)
})

test_that("stat_specifiedevals() is zero for standarised sample", {
  set.seed(13131)
  Ysample <- rsymm(50, 3)
  Ysample <- lapply(Ysample, `+`, diag(c(3,2,1)))
  Ystdsample <- standardise_specifiedevals(Ysample, c(3,2,1))
  expect_equal(stat_specifiedevals(Ystdsample, c(3,2,1)), 0)
})

