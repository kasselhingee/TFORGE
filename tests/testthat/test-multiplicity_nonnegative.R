# Note that it seems really hard to simulate with fixed trace and non-negative eigenvalues, with mean with specified eigenvalues
# Instead I've used the empirical distribution trick with a very large sample
evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
mult <- c(3,2,1,1)
# creating an empirical population that satisfies the null
set.seed(5)
bigtmp <- rsymm_norm(10000, diag(evals), sigma = diag(1, sum(mult) * (sum(mult) + 1) / 2)/10)
nonnegevals <- apply(bigtmp, 1, function(v){min(eigen_desc(invvech(v))$values)}) >= 0 #time consuming truncation
bigtmp <- as_fsm(bigtmp[nonnegevals, ])
big <- normalise_trace(bigtmp)
big <- as_fsm(big)
nullmean <- multiplicity_nullmean(mmean(big), mult)
wts <- elwts_fixedtrace(big, nullmean, maxit = 25)
expect_true(wtsokay(wts))
specialsample <- function(size){
  idx <- sample.int(nrow(big), size = size, prob = wts, replace = TRUE)
  out <- big[idx, ]
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
    test_multiplicity_nonnegative(Ysample, mult = mult, B = 100)$pval
  })
  
  # qqplot(vals, y = runif(1000))
  expect_gt(suppressWarnings(ks.test(vals, "punif"))$p.value, 0.05) #above 0.2 if above two thoroughness measures taken
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

test_that("test rejects some incorrect hypotheses", {
  set.seed(13321)
  Ysample <- specialsample(100)
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
