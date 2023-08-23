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

test_that("stat_commonevals_ksample() has correct null distribution", {
  set.seed(13131)
  vals <- replicate(100, {
    Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE)
    stat_commonevals_ksample(Ysamples)
  })
  
  # qqplot(vals, y = rchisq(1000, df = (5-1)*3))
  res <- ks.test(vals, "pchisq", df = (5-1)*3)
  expect_gt(res$p.value, 0.2)
})

test_that("test_specifiedevals() from NULL has uniform p-values", {
  pvals <- vapply(13 + (1:100), function(seed){
    set.seed(seed)
    Ysample <- rsymm(50, diag(c(3,2,1)))
    res <- test_specifiedevals(Ysample, c(3,2,1), 100)
    res$pval}, FUN.VALUE = 1.1)
  # qqplot(pvals, y = runif(1000))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})


test_that("test_commonevals() from NULL mst has uniform p-values", {
  set.seed(13)
  pvals <- pbreplicate(100, {
    Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE)
    res <- test_commonevals(Ysamples, B = 100)
    res$pval
  }, cl = 2)
  qqplot(pvals, y = runif(1000))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
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
  set.seed(134)
  Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE) 
  set.seed(1345)
  res1 <- test_commonevals(Ysamples, 10)
  set.seed(1345)
  res2 <- test_commonevals(Ysamples, 10)
  expect_equal(res1, res2)
})


test_that("stat_specifiedevals() is zero for standarised sample", {
  set.seed(13131)
  Ysample <- rsymm(50, diag(c(3,2,1)))
  Ystdsample <- standardise_specifiedevals(Ysample, c(3,2,1))
  expect_equal(stat_specifiedevals(Ystdsample, c(3,2,1)), 0)
})
