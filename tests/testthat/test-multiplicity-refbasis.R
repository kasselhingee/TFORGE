test_that("project_basis independent of subspace basis", {
  p <- 4
  # Reference basis:
  # set.seed(1)
  refbasis <- mclust::randomOrthogonalMatrix(p, p)
  
  # eigenspace size with arbitrary basis
  q <- 2
  # set.seed(2)
  q_abasis <- mclust::randomOrthogonalMatrix(p, p)[, 1:q]
  bas1 <- project_basis(q_abasis, refbasis)
  # rotation in eigenspace
  # set.seed(3)
  R <-  mclust::randomOrthogonalMatrix(q, q)
  # apply rotation to q_abasis
  q_abasis2 <- q_abasis %*% R
  bas2 <- project_basis(q_abasis2, refbasis)
  expect_equal(bas1, bas2)
})

test_that("project_basis onto full space preserves ref", {
  p <- 3
  # Reference basis:
  set.seed(1)
  refbasis <- mclust::randomOrthogonalMatrix(p, p)
  q_abasis <- mclust::randomOrthogonalMatrix(p, p)
  bas1 <- project_basis(q_abasis, refbasis)
  expect_equal(bas1, refbasis)
})

test_that("project_basis depends on reference", {
  p <- 3
  q <- 2
  # Reference basis:
  set.seed(1)
  refbasis <- mclust::randomOrthogonalMatrix(p, p)
  refbasis2 <- mclust::randomOrthogonalMatrix(p, p)
  q_abasis <- mclust::randomOrthogonalMatrix(p, p)[, 1:q]
  bas1 <- project_basis(q_abasis, refbasis)
  bas2 <- project_basis(q_abasis, refbasis2)
  expect_error(expect_equal(bas1, bas2))
})

# project ref basis to eigenspace
projbasis1 <- q_abasis %*% t(q_abasis) %*% refbasis
svd(projbasis1)

# rotation in eigenspace
set.seed(2)
R <-  mclust::randomOrthogonalMatrix(q, q)
# apply rotation to q_abasis
q_abasis2 <- q_abasis %*% R
projbasis2 <- q_abasis2 %*% t(q_abasis2) %*% refbasis
all.equal(projbasis1, projbasis2)
all.equal(q_abasis2 %*% t(q_abasis2), q_abasis %*% t(q_abasis)) #projection matrices for a space don't care about how subspace is parameterised.
# Therefore using the SVD/eigenvectors of the projection matrix is a fine thing to do in the sense that they are independent of the parameterisation of the space.
# But I think it is only looking stable because of the computation accuracy

projmat <- q_abasis %*% t(q_abasis)
all.equal(svd(projmat), svd(projmat %*% refbasis))
alignedsvd <- svd(projmat %*% diag(1 + (p:1)/10000) %*% refbasis, nu = q, nv = q)
all.equal(alignedsvd$u %*% t(alignedsvd$u), projmat)
# Fails: all.equal(alignedsvd$v %*% t(alignedsvd$v), projmat)

# alignedsvd$u are candidate axes for the subspace.


