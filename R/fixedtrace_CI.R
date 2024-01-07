
# For an ellipse size a and b on the x and y axis
# return the (x,y) location corresponding to the ray at the angle intersecting the unit circle, then scaled to have the correct a or b
# that is gives the (x,y) such that (x/a, y/b) (which is on the unit circle) corresponds to the given angle
# angle may be a vector of angles
regularellipse <- function(angle, a, b){
  x <- a * cos(angle)
  y <- b * sin(angle)
  return(cbind(x = x, y = y))
}

# convert points from regular ellipse to the ft plane with rotation by eigenvectors of Omega
# @param pts from regularellipse()
# @param evecs Eigenvectors of Omega
ellipseinftplane <- function(pts, evecs){
  rotpts <- evecs %*% t(pts)
  Hfull <- helmert(ncol(pts) + 1)
  cbind(0, t(rotpts)) %*% Hfull
}

# combine ellipseinftplane and regular ellipse to get locations around the mean
ellipseftcentre <- function(angle, a, b, evecs, ctrevals){
  locs2d <- regularellipse(angle, a = a, b = b)
  locsplane <- ellipseinftplane(locs2d, evecs = evecs)
  locsplanearoundcenter <- t(ctrevals - t(locsplane))
  return(locsplanearoundcenter)
}

#' @title Confidence region for 3x3 Tensors with Fixed Trace Constraint
#' @param x A sample of 3x3 tensors.
#' @param alpha Significance level
#' @param B Number of bootstrap resamples.
#' @param pts Number of points on the boundary of the region to compute.
#' @return A list:
#' + `est`: the eigenvalues of the mean tensor
#' + `boundary`: npts x 3 matrix giving the boundary of the region (each row corresponds to a point on the boundary and the columns are the first, second and final eigenvalue).
#' + `inregion`: A function that accepts a vector of eigenvalues and returns `TRUE` when the eigenvalues are in the confidence region and `FALSE` otherwise.
#' + `Omega`: The estimated covariance of the (projected) eigenvalues
#' + `threshold`: The threshold on the statistic used.
#' @export
conf_fixedtrace <- function(x, alpha = 0.05, B = 1000, npts = 1000){
  stopifnot(ncol(x) == 6)
  p <- 3
  x <- as.sst(x)

  # sample mean
  av <- mmean(x)
  av_ess <- eigen_desc(av)
  av_eval <- av_ess$values
  size <- nrow(x)

  # resampling and computing stat_fixedtrace()
  res <- bootresampling(x, x, stat = stat_fixedtrace, B = B, evals = av_eval)
  statthreshold <- quantile(res$nullt, probs = 1-alpha, names = FALSE, type = 1)

  # now compute boundary of region
  Omega <- cov_evals_ft(x, evecs = av_ess$vectors, av = av)
  Omega_ess <- eigen_desc(Omega)
  a <- sqrt(statthreshold * Omega_ess$values[1]/size)
  b <- sqrt(statthreshold * Omega_ess$values[2]/size)
  stopifnot(all(Omega_ess$values >= 0))
  bdrypts <- ellipseftcentre(angle = seq(0, 2*pi, length.out = npts),
                         a = a,
                         b = b,
                         evecs = Omega_ess$vectors,
                         ctreval = av_eval
                         )

  # also create a function that tests whether inside the region using the statistic directly
  H <- helmertsub(p)
  inregion <- function(evals){
    statval <- size * t(av_eval - evals) %*% t(H) %*% Omega_ess$vectors %*% diag(1/Omega_ess$values) %*% t(Omega_ess$vectors) %*% H %*% (av_eval - evals)
    return(drop(statval) <= statthreshold)
  }

  # now check if close to the descending order boundary by solving for the intersection of the boundary lines and the ellipse
  # the sqrt(6) comes from the definition of the helmertsub matrix
  # see ConfidenceRegionsFT.pdf EQ1 and EQ2 for obtaining the equations
  # first two evals:
  ft <- sum(av_eval)
  A <- t(Omega_ess$vectors) %*% (H %*% av_eval + (2*ft/sqrt(6)) * rbind(0, 1))
  B <- sqrt(6)  * t(Omega_ess$vectors)  %*%  rbind(0, 1)
  # b^2 - 4ac >= 0
  invLambda <- diag(1/Omega_ess$values)
  e1e2intersect <- drop(
  (t(A) %*% invLambda %*% B + t(B) %*%  invLambda %*% A)^2 -
    4 * (t(B) %*% invLambda %*% B) * ((t(A) %*% invLambda %*% A) - statthreshold/size)
  ) >= 0
  
  # second two evals
  A <- t(Omega_ess$vectors) %*% (H %*% av_eval + (ft/2) * rbind(1/sqrt(2), 1/sqrt(6)))
  B <- (3/2) * t(Omega_ess$vectors) %*% rbind(1/sqrt(2), 1/sqrt(6))
  e2e3intersect <- drop(
    (t(A) %*% invLambda %*% B + t(B) %*%  invLambda %*% A)^2 -
      4 * (t(B) %*% invLambda %*% B) * ((t(A) %*% invLambda %*% A) - statthreshold/size)
  ) >= 0
  
  
  if (e1e2intersect | e2e3intersect){
    desc <- switch(e1e2intersect + 2 * e2e3intersect,
           "largest two eigenvalues",
           "smallest two eigenvalues",
           "largest two eigenvalues or smallest two eigenvalues")
    warning(sprintf("Confidence region includes eigenvalues where the %s are not in descending order.", desc))
  } 

  return(list(
    est = av_eval,
    boundary = bdrypts,
    inregion = inregion,
    Omega = Omega,
    threshold = statthreshold,
    bootstat = res$nullt
  ))
}
