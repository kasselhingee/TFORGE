
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
  t(Hfull) %*% rbind(0, rotpts)
}
