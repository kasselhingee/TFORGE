
test_that("multisample() gives correct sizes", {})

test_that("multisample() listens to weights", {
  set.seed(13)
  Ysamples <- list(rsymm_norm(5, diag(c(3,2,1))),
                rsymm_norm(6, diag(c(3,2,1))),
                rsymm_norm(7, diag(c(3,2,1))))
  
  #sample just the first
  w <- list(c(1, rep(0,4)),
            c(1, rep(0,5)),
            c(1, rep(0,6))
            )
  
  newY <- multisample(Ysamples, w = w)
  expect_equal(lapply(Ysamples, "[[", 1), lapply(newY, "[[", 1))
  expect_equal(lapply(Ysamples, "[[", 1), lapply(newY, "[[", 5))
})