
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
  locsplanearoundcenter <- t(ctrevals - t(locs))
  return(locsplanearoundcenter)
}

#' @title Confidence for 3x3 Tensors with Fixed Trace Constraint
#' @param x A sample of 3x3 tensors.
#' @param alpha Significance level
#' @param B Number of bootstrap resamples.
#' @param pts Number of points on the boundary of the region to compute.
#' @export
conf_fixedtrace <- function(x, alpha = 0.05, B = 1000, npts = 1000){
  stopifnot(ncol(x) == 6)
  p <- 3
  x <- as.sst(x)

  # sample mean
  av <- mmean(x)
  av_ess <- eigen_desc(av)
  av_eval <- av_ess$values

  # resampling and computing stat_fixedtrace()
  res <- bootresampling(x, x, stat = stat_fixedtrace, B = B, evals = av_eval)
  statthreshold <- quantile(res$nullt, probs = 1-alpha)

  # now compute boundary of region
  Omega <- cov_evals_ft(x, evecs = av_ess$vectors, av = av)
  Omega_ess <- eigen_desc(Omega)
  stopifnot(all(Omega_ess$values >= 0))
  bdrypts <- ellipseftcentre(angle = seq(0, 2*pi, length.out = npts),
                         a = sqrt(statthreshold * Omega_ess$values[1]),
                         b = sqrt(statthreshold * Omega_ess$values[2]),
                         evecs = Omega_ess$vectors,
                         ctreval = av_eval
                         )

  # also create a function that tests whether inside the region using the statistic directly
  H <- helmertsub(p)
  inregion <- function(evals){
    statval <- t(av_eval - evals) %*% t(H) %*% Omega_ess$vectors %*% diag(1/Omega_ess$values) %*% t(Omega_ess$vectors) %*% (av_eval - evals)
    return(statval <= statthreshold)
  }

  # now check if close to the descending order boundary empirically
  e1e2 <- optim(par = c(0, pi/2), 
        fn = function(angle){
      bdrypt <- ellipseftcentre(angle = angle,
                         a = sqrt(statthreshold * Omega_ess$values[1]),
                         b = sqrt(statthreshold * Omega_ess$values[2]),
                         evecs = Omega_ess$vectors,
                         ctreval = av_eval
                         )
      (bdrypt[1] - bdrypt[2])^2 
      },
  method = "CG")
  e2e3 <- optim(par = c(0, pi/2), 
        fn = function(angle){
      bdrypt <- ellipseftcentre(angle = angle,
                         a = sqrt(statthreshold * Omega_ess$values[1]),
                         b = sqrt(statthreshold * Omega_ess$values[2]),
                         evecs = Omega_ess$vectors,
                         ctreval = av_eval
                         )
      (bdrypt[2] - bdrypt[3])^2 
      },
  method = "CG")
  if ((e1e2$value < sqrt(.Machine$double.eps)) |
      (e2e3$value < sqrt(.Machine$double.eps)) ){
    warning("Confidence region intersects the change in order boundary, the confidence region should not be trusted in this situation")
  }

  return(list(
    boundary = bdrypts,
    inregion = inregion
  ))
}
