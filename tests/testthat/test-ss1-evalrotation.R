# see if it is possible to find a rotation matrix that gives desired eigenvalues
test_that("a rotation of eigenvalues to the null exists", {
  set.seed(1333)
  Y <- rsymm_norm(50, diag(c(3,2,1)), sigma = diag(0.7, 6))
  Y <- normL2evals_sst(Y)
  
  # rotate eigenvalues so that mean of Y has eigenvalues c(1,1,1)/sqrt(3)
  nullevals <- normalise_ss1_vec(c(1,1,1))
  rot_evals <- function(Y, rot){
    Yrot <- apply(Y, 1, FUN = function(v){
      es <- eigen_desc(invvech(v))
      es$values <- drop(rot %*% es$values)
      makeSymmetric(es$vectors %*% diag(es$values) %*% t(es$vectors))
    }, simplify = FALSE)
    return(as_fsm(Yrot))
  }
  tominimise <- function(x, Y, nullevals){ #x is vector of length 3
    rot <- sphm:::cayley(x) #convert x to an orthogonal matrix
    # apply rotation to sample
    Yrot <- rot_evals(Y, rot)
    # get evals of average
    av <- mmean(Yrot)
    avevals <- eigen_desc(av)$values
    avevals <- normalise_ss1_vec(avevals)
    acos(sum(avevals * nullevals))
  }
  tominimise(runif(3, -1, 1), Y, nullevals)
  myoptimise <- function(start){
    optim(start, tominimise,
        Y = Y, nullevals = nullevals,
        control = list(abstol = 1E-2,
                       reltol = 0,
                       maxit = 200))
  }
  #many starting locations
  starts <- expand.grid(seq(-1, 1, by = 0.5), seq(-1, 1, by = 0.5), seq(-1, 1, by = 0.5))
  res <- pbapply::pbapply(starts, 1, myoptimise)
  summary(vapply(res, "[[", "value", FUN.VALUE = .1))
  # the best arent that close to the nullevals
  res[[2]]$value
  bestavnull <- mmean(rot_evals(Y, sphm:::cayley(res[[2]]$par)))
  expect_equal(eigen_desc(bestavnull)$values, nullevals, tolerance = 1E-2) #vs weighted bootstrap which has exactly the null evals
})
