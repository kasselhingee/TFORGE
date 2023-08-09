test_that("test_commonevals() doesn't reject for simulation of multi sample from null", {
  set.seed(13)
  Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE)

  res <- test_commonevals(Ysamples, 20)
  expect_gt(res$pval, 0.2)
})

test_that("test_commonevals() reject for simulation of multi sample not from null", {
  set.seed(13)
  Ysamples <- list(
    rsymm(50, diag(c(3,2,1))),
    rsymm(50, diag(c(4,3,2)))
  )
  res <- test_commonevals(Ysamples, 20)
  expect_lt(res$pval, 0.05)
})

test_that("test_commonevals() repeats under set.seed()", {
  set.seed(34)
  Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE) 
  set.seed(134)
  res1 <- test_commonevals(Ysamples, 10)
  set.seed(134)
  res2 <- test_commonevals(Ysamples, 10)
  expect_equal(res1, res2)
})
