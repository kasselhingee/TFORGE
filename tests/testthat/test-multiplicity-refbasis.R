test_that("project_basis independent of subspace basis", {
  p <- 5
  # Reference basis:
  # set.seed(1)
  refbasis <- runifortho(p)
  
  # eigenspace size with arbitrary basis
  q <- 2
  # set.seed(2)
  q_abasis <- runifortho(p)[, 1:q]
  bas1 <- project_basis(q_abasis, refbasis)
  # rotation in eigenspace
  # set.seed(3)
  R <-  runifortho(q)
  # apply rotation to q_abasis
  q_abasis2 <- q_abasis %*% R
  bas2 <- project_basis(q_abasis2, refbasis)
  expect_equal(abs(t(bas1) %*% bas2), diag(nrow = q))
})

test_that("project_basis onto full space preserves ref", {
  p <- 4
  # Reference basis:
  set.seed(1)
  refbasis <- runifortho(p)
  q_abasis <- runifortho(p)
  bas1 <- project_basis(q_abasis, refbasis)
  expect_equal(abs(t(bas1) %*% refbasis), diag(nrow = p))
})

test_that("project_basis depends on reference", {
  p <- 5
  q <- 2
  # Reference basis:
  set.seed(1)
  refbasis <- runifortho(p)
  refbasis2 <- runifortho(p)
  q_abasis <- runifortho(p)[, 1:q]
  bas1 <- project_basis(q_abasis, refbasis)
  bas2 <- project_basis(q_abasis, refbasis2)
  expect_error(expect_equal(abs(t(bas1) %*% bas2), diag(nrow = q), tolerance = 0.1))
})

test_that("project_basis returns vector for 1d subspace", {
  p <- 3
  q <- 1
  # Reference basis:
  set.seed(1)
  refbasis <- runifortho(p)
  q_abasis <- runifortho(p)[, 1:q, drop = FALSE]
  bas1 <- project_basis(q_abasis, refbasis)
  expect_equal(drop(t(bas1) %*% q_abasis), 1)
})

test_that("stat_multiplicity() is slightly different with random, cannonical and other refbasis", {
  set.seed(13321)
  Ysample <- rsymm(20, diag(c(rep(3, 3), rep(2, 2), 1, 0.5)),
                   sigma = 0.1 * diag(7 * (7 + 1)/2))
  set.seed(3654)
  statr <- stat_multiplicity(Ysample, mult = c(2,3,1,1), refbasis = "random")
  statc <- stat_multiplicity(Ysample, mult = c(2,3,1,1), refbasis = diag(1, 7))
  set.seed(3)
  stata <- stat_multiplicity(Ysample, mult = c(2,3,1,1), refbasis = runifortho(7))
  expect_error(expect_error(statr, statc))
  expect_error(expect_error(stata, statc))
})


test_that("stat has INcorrect null distribution using sample evecs", {
  skip_on_cran()
  set.seed(1331)
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  vals <- replicate(1000, {
    Ysample <- rsymm_norm(100, diag(evals), sigma = diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    evecs <- eigen_desc(mmean(Ysample))$vectors
    stat_multiplicity(Ysample, mult = mult, refbasis = evecs)
  })
  
  qqplot(vals, y = rchisq(1000, df = sum(mult-1)))
  expect_lt(ks.test(vals, "pchisq", df = sum(mult-1))$p.value, 0.001)
})

test_that("test has uniform distribution with misuse of sample evecs", {
  skip_on_cran() #test very slow
  evals <- c(rep(3, 3), rep(2, 2), 1, 0.5)
  mult <- c(3,2,1,1)
  set.seed(5)#set.seed(1331)
  vals <- pbapply::pbreplicate(100, { #1000 for more thorough
    Ysample <- rsymm_norm(100, diag(evals), sigma = 0.001 * diag(1, sum(mult) * (sum(mult) + 1) / 2) )
    evecs <- eigen_desc(mmean(Ysample))$vectors
    test_multiplicity(Ysample, mult = mult, B = 20, refbasis = evecs)$pval
  })
  
  # qqplot(vals, y = runif(1000))
  expect_lt(suppressWarnings(ks.test(vals, "punif"))$p.value, 0.0001)
})

