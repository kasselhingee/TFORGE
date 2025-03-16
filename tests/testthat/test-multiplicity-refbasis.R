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


