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
