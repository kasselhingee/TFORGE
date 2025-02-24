# Note that it seems really hard to simulate with fixed trace and non-negative eigenvalues, with mean with specified eigenvalues
# Instead I've used the empirical distribution trick with a very large sample
logevals <- c(rep(3, 3), rep(2, 2), 1, 0.5)-5
mult <- c(3,2,1,1)
# creating an empirical population that satisfies the null
set.seed(5)
bigtmp <- rsymm_norm(10000, diag(logevals), sigma = diag(1, sum(mult) * (sum(mult) + 1) / 2))
nonnegevals <- t(apply(bigtmp, 1, function(v){
  es <- eigen_desc(invvech(v))
  vech(es$vectors %*% diag(exp(es$values)) %*% t(es$vectors))
  }))
nonnegevals <- as_fsm(nonnegevals)
mmean(nonnegevals)
cov_evals_est(nonnegevals)
aves <- eigen_desc(mmean(nonnegevals))
nullmean <- multiplicity_nullmean(mmean(nonnegevals), mult)
# expect_equal(t(aves$vectors) %*% nullmean %*% aves$vectors,
#              diag(eigen_desc(nullmean)$values))

wts <- elwts_fixedtrace(nonnegevals, nullmean, maxit = 25)
# expect_true(wtsokay(wts))
specialsample <- function(size){
  idx <- sample.int(nrow(nonnegevals), size = size, prob = wts, replace = TRUE)
  out <- nonnegevals[idx, ]
  class(out) <- c("TFORGE_fsm", class(out))
  return(out)
}

test_that("stat has correct null distribution", {
  # simulate check distribution of stat:
  set.seed(10)
  vals <- replicate(1000, { #asymptotics work slower with the specialsample() distribution than a normal hence larger sample size here
    Ysample <- specialsample(200)
    stat_multiplicity(Ysample, mult = mult)
  })
  
  # qqplot(vals, y = rchisq(1000, df = sum(mult-1)))
  expect_gt(ks.test(vals, "pchisq", df = sum(mult-1))$p.value, 0.2)
})

test_that("test has uniform distribution", {
  skip_on_cran() #test very slow
  set.seed(5)#set.seed(1331)
  vals <- pbapply::pbreplicate(1000, { #1000 for more thorough
    Ysample <- specialsample(200)
    # B = 100 for more thorough
    test_multiplicity_nonnegative(Ysample, mult = mult, B = 1000)$pval
  })
  
  # qqplot(vals, y = runif(1000))
  expect_gt(suppressWarnings(ks.test(vals, "punif"))$p.value, 0.05) #above 0.2 if above two thoroughness measures taken
})

test_that("test size at small samples: Error not uniform.", {
  set.seed(5)#set.seed(1331)
  vals <- pbapply::pbreplicate(1000, { #1000 for more thorough
    Ysample <- specialsample(100) #at n = 15 and 30: null mean never in convex hull
    res <- test_multiplicity_nonnegative(Ysample, mult = mult, B = 1000)
    res[c("pval", "B")]
  }, cl = 3)
  sum(is.na(unlist(vals[2, ])))
  pvals <- unlist(vals[1, ])
  mean(pvals <= 0.05)
  qqplot(pvals, y = runif(1000))
  ks.test(pvals, "punif")
})

test_that("chisq: test has uniform distribution", {
  set.seed(3) #set.seed(1331)
  vals <- replicate(100, { #1000 for more thorough
    Ysample <- specialsample(100)
    test_multiplicity_nonnegative(Ysample, mult = mult, B = "chisq")$pval
  })
  
  # qqplot(vals, y = runif(1000))
  expect_gt(suppressWarnings(ks.test(vals, "punif"))$p.value, 0.15)
})

test_that("test rejects some incorrect hypotheses simulated from unconstrained case", {
  set.seed(13321)
  Ysample <- rsymm(100, diag(c(rep(3, 3), rep(2, 2), 1, 0.5)),
                   sigma = 0.1 * diag(7 * (7 + 1)/2))
  set.seed(3653)
  res <- test_multiplicity_nonnegative(Ysample, mult = c(3,2,1,1), 100)
  expect_gt(res$pval, 0.1)
  
  res <- test_multiplicity_nonnegative(Ysample, mult = c(3,1,1,1,1), 100)
  expect_gt(res$pval, 0.1)
  res <- test_multiplicity_nonnegative(Ysample, mult = c(2,1,1,1,1,1), 100)
  expect_gt(res$pval, 0.1)
  
  # set.seed(3542)
  # pvals_2311 <- replicate(100, test_multiplicity(Ysample, mult = c(2,3,1,1), 100)$pval, cl = 3)
  set.seed(3542)
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(2,3,1,1), 100)$pval, 0.05)
  set.seed(35423) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(2,2,2,1), 100)$pval, 0.05)
  set.seed(35424) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(4,1,1,1), 100)$pval, 0.05)
  set.seed(35425) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(3,1,3), 100)$pval, 0.05)
  set.seed(35426) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(3,3,1), 100)$pval, 0.05)
  
  set.seed(35427) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(2,3,2), 100)$pval, 0.05)
  set.seed(3541) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(3,2,2), 100)$pval, 0.05)
  # note that at B=100 there is still a lot a randomness in the output pvalue
})

test_that("test rejects some incorrect hypotheses p = 3", {
  set.seed(13321)
  Ysample <- rsymm(30, diag(c(2,2,1)-3), sigma = diag(3 * (3+1) /2))
  nonnegevals <- t(apply(Ysample, 1, function(v){
    es <- eigen_desc(invvech(v))
    vech(es$vectors %*% diag(exp(es$values)) %*% t(es$vectors))
  }))
  nonnegevals <- as_fsm(nonnegevals)
  eigen_desc(mmean(nonnegevals))$values
  set.seed(3654)
  res <- test_multiplicity_nonnegative(nonnegevals, mult = c(2,1), B = 1000)
  expect_gt(res$pval, 0.1)
  
  set.seed(3542)
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(1, 2), B = 1000)$pval, 0.05)
  
  set.seed(3541)
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(3), B = 1000)$pval, 0.05)
})

test_that("test rejects some incorrect hypotheses", {
  set.seed(13321)
  Ysample <- specialsample(30)
  
  set.seed(1)
  expect_gt(test_multiplicity_nonnegative(Ysample, mult = mult, 1000)$pval, 0.05)
  
  set.seed(3542)
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(2,3,1,1), 1000)$pval, 0.05)
  set.seed(35423) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(2,2,2,1), 1000)$pval, 0.05)
  set.seed(35424) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(4,1,1,1), 1000)$pval, 0.05)
  set.seed(35425) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(3,1,3), 1000)$pval, 0.05)
  set.seed(35426) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(3,3,1), 1000)$pval, 0.05)
  
  set.seed(35427) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(2,3,2), 1000)$pval, 0.05)
  set.seed(3541) 
  expect_lt(test_multiplicity_nonnegative(Ysample, mult = c(3,2,2), 1000)$pval, 0.05)
  # note that at B=100 there is still a lot a randomness in the output pvalue
})

test_that("fixed trace from projection preserved by standardisation and ignored by stat", {
  set.seed(134)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  Ysample <- rsymm_norm(1000, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2))
  Ysample_n <- project_trace(Ysample)
  
  # check that only first element of Hevals  is changed by trace fix
  evals <- eigen_desc(mmean(Ysample))$values
  evals_n <- eigen_desc(mmean(Ysample_n))$values
  expect_equal((helmert(sum(mult)) %*% evals)[-1],
    (helmert(sum(mult)) %*% evals_n)[-1])
  
  std <- standardise_multiplicity(Ysample, mult)
  std_n <- standardise_multiplicity(Ysample_n, mult)
  expect_true(all.equal(as_fsm(apply(std, 1, function(m){projtrace_matrix(invvech(m))}, simplify = FALSE)), 
                        std_n, check.attributes = FALSE))
  
  set.seed(123)
  stat_orig <- stat_multiplicity(Ysample, mult = mult)
  set.seed(123)
  stat_n <- stat_multiplicity(Ysample_n, mult = mult)
  expect_equal(stat_orig, stat_n)
})
