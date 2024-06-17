test_that("multisample() listens to weights and gives correct sizes", {
  set.seed(13)
  Ysamples <- as_flat(list(rsymm_norm(5, diag(c(3,2,1))),
                rsymm_norm(6, diag(c(3,2,1))),
                rsymm_norm(7, diag(c(3,2,1)))))
  
  #sample just the first
  w <- list(c(1, rep(0,4)),
            c(1, rep(0,5)),
            c(1, rep(0,6))
            )
  
  newY <- multisample(Ysamples, prob = w)
  expect_equal(lapply(newY, "[[", 1), lapply(Ysamples, "[[", 1))
  expect_equal(lapply(newY, "[[", 5), lapply(Ysamples, "[[", 1))
  expect_equal(lapply(newY, nrow), lapply(Ysamples, nrow))
  
  set.seed(315)
  newY <- multisample(Ysamples)
  expect_equal(lapply(newY, nrow), lapply(Ysamples, nrow))
  expect_error(expect_equal(lapply(newY, "[[", 1), lapply(Ysamples, "[[", 1)))
})
