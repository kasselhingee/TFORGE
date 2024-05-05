test_that("stat for sst has correct null distribution", {
  vals <- vapply(1:100, function(seed){
    set.seed(seed)
    Ysample <- rsymm(50, diag(c(3,2,1)))
    stat <- stat_unconstrained(Ysample, c(3,2,1))
    expect_equal(as.numeric(stat), stat_specifiedevals(Ysample, c(3,2,1)))
    stat
  }, FUN.VALUE = 1.3)

  # qqplot(vals, y = rchisq(1000, df = 3))
  res <- ks.test(vals, "pchisq", df = 3)
  expect_gt(res$p.value, 0.2)
})

test_that("stat for sst, specified evecs, has correct null distribution", {
  vals <- vapply(1:100, function(seed){
    set.seed(seed)
    Ysample <- rsymm(50, diag(c(3,2,1)))
    stat <- stat_unconstrained(Ysample, evals = c(3,2,1), evecs = diag(1, 3))
    expect_equal(as.numeric(stat), stat_specifiedevals(Ysample, c(3,2,1), evecs = diag(1, 3)))
    stat
  }, FUN.VALUE = 1.3)
  
  # qqplot(vals, y = rchisq(1000, df = 3))
  res <- ks.test(vals, "pchisq", df = 3)
  expect_gt(res$p.value, 0.2)
})

test_that("stat for mst has correct null distribution", {
  set.seed(13131)
  vals <- replicate(100, {
    Ysamples <- lapply(c(2000,2000,1000,1000,1000), function(n) rsymm(n, diag(c(3,2,1))))
    stat <- stat_unconstrained(Ysamples)
    stat
  })
  
  # qqplot(vals, y = rchisq(1000, df = (5-1)*3))
  res <- ks.test(vals, "pchisq", df = (5-1)*3)
  expect_gt(res$p.value, 0.2)
})

test_that("stat for mst w specified evecs has INcorrect null distribution", {
  set.seed(1311)
  vals <- replicate(100, {
    Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE)
    suppressWarnings(stat_unconstrained(Ysamples, evecs = diag(1, 3)))
  })
  
  # qqplot(vals, y = rchisq(1000, df = (5-1)*3))
  res <- ks.test(vals, "pchisq", df = (5-1)*3)
  expect_lt(res$p.value, 0.2)
})

test_that("test sst from NULL has uniform p-values", {
  pvals <- vapply(13 + (1:100), function(seed){
    set.seed(seed)
    Ysample <- rsymm(50, diag(c(3,2,1)))
    set.seed(seed+1)
    res <- test_unconstrained(Ysample, c(3,2,1), B = 100)
    set.seed(seed+1)
    res2 <- test_specifiedevals(Ysample, c(3,2,1), B = 100)
    expect_equal(res[c("pval", "nullt")], res2[c("pval", "nullt")])
    res$pval}, FUN.VALUE = 1.1)
  # qqplot(pvals, y = runif(1000))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})

test_that("test from NULL mst has uniform p-values", {
  set.seed(13)
  pvals <- vapply(13 + (1:100), function(seed){
    set.seed(seed)
    Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE)
    set.seed(seed+1)
    res <- test_unconstrained(Ysamples, B = 100)
    # set.seed(seed+1)
    # res2 <- test_commonevals(Ysamples, B = 100)
    # expect_equal(res[c("pval", "nullt")], res2[c("pval", "nullt")])
    res$pval
  }, FUN.VALUE = 1.3)
  # qqplot(pvals, y = runif(1000))
  res <- suppressWarnings({ks.test(pvals, "punif")})
  expect_gt(res$p.value, 0.05)
})


test_that("test reject for simulation of multi sample not from null", {
  set.seed(13)
  Ysamples <- list(
    rsymm(50, diag(c(3,2,1))),
    rsymm(50, diag(c(4,3,2)))
  )
  res <- test_unconstrained(Ysamples, B = 20)
  expect_lt(res$pval, 0.05)
})

test_that("test repeats under set.seed()", {
  set.seed(134)
  Ysamples <- replicate(5, rsymm(50, diag(c(3,2,1))), simplify = FALSE) 
  set.seed(1345)
  res1 <- test_unconstrained(Ysamples, B = 10)
  set.seed(1345)
  res2 <- test_unconstrained(Ysamples, B = 10)
  expect_equal(res1, res2)
})


test_that("stat is zero for standarised sample", {
  set.seed(13131)
  Ysample <- rsymm(50, diag(c(3,2,1)))
  Ystdsample <- standardise_specifiedevals(Ysample, c(3,2,1))
  expect_equal(as.numeric(stat_unconstrained(Ystdsample, c(3,2,1))), 0)
})
