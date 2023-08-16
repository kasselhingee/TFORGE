test_that("stat_specifiedevals() has correct null distribution", {
  set.seed(1365)
  vals <- replicate(100, {
    Ysample <- rsymm(50, diag(c(3,2,1)))
    stat_specifiedevals(Ysample, c(3,2,1))
  })

  qqplot(vals, y = rchisq(1000, df = 3))
  res <- ks.test(vals, "pchisq", df = 3)
  expect_gt(res$p.value, 0.2)
})

test_that("stat_specifiedevals() doesn't reject for simulation of single sample from null", {
  set.seed(13131)
  Ysample <- rsymm(50, diag(c(3,2,1)))
  res <- test_specifiedevals(Ysample, c(3,2,1), 100)
  expect_gt(res$pval, 0.2)
})

test_that("stat_specifiedevals() is zero for standarised sample", {
  set.seed(13131)
  Ysample <- rsymm(50, diag(c(3,2,1)))
  Ystdsample <- standardise_specifiedevals(Ysample, c(3,2,1))
  expect_equal(stat_specifiedevals(Ystdsample, c(3,2,1)), 0)
})

