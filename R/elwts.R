# for a single sample compute corresponding weights that lead to emp.lik.
elwts_fixedtrace <- function(x, nullmean, maxit = 25){
  scelres <- emplik(x, vech(nullmean), itermax = maxit)
  wts <- as.vector(scelres$wts) * nrow(x)
  if (!isTRUE(scelres$converged)){
    wts <- -1 * wts
  }
  return(wts)
}


#function finds the best c and weights such that weighted average of data is c.mu and 
#empirical likelihood maximised
# note that the profile likelihood function (result of emplik()) has convex superlevel sets according to Theorem 3.2 (Owen 2001).
# So there is unique minimum value to the problem where the mean lies on a line.
# @returns If the method doesn't converge then the negative of the weights is returned
# @param ms A single sample of symmetric tensors
# @param mu Proposed mean up-to-constant c. It is assumed to have an attribute "c_range" range gives a range of values of c, passed to `optimize()`
opt_el.test <- function(ms, mu, maxit = 25){
  if (attr(mu, "c_range")[["min"]] > attr(mu, "c_range")[["max"]]){
    return(rep(0, nrow(ms))) #no pluasible values of c - avoids error triggered in optimise
  }
  bestmult <- stats::optimise(f = function(x){#optim warns that Nelder-Mead unreliable on 1 dimension so using Brent here instead
    scelres <- emplik(ms, x*vech(mu), itermax = maxit)
    return(-scelres$logelr)
  },
  lower = attr(mu, "c_range")[["min"]], 
  upper = attr(mu, "c_range")[["max"]]) 
  scelres <- emplik(ms, bestmult$minimum*vech(mu), itermax = maxit)
  wts <- as.vector(scelres$wts) * nrow(ms) #the multiple here is to match the weights put out by emplik::el.test()

  # check and warn about results
  if (!isTRUE(scelres$converged)){
    wts <- -1 * wts
  }
  wts
}

# A function for checking that a list of weights is okay
wtsokay <- function(wts){
  if (is.atomic(wts) && is.numeric(wts)){wts <- list(wts)} #for single sample situations a vector of weights can be passed - it will be treated a list with one element.
  if (any(vapply(wts, sum, FUN.VALUE = 1.3) < 0)){ #negative weights
    warning("The proposed null mean appears to be outside the convex hull for at least one sample (emplik() did not converge)" , call. = FALSE)
    return(FALSE)
  }
  if (any(vapply(wts, sum, FUN.VALUE = 1.3) < 0.5)){# weights sum to zero
    warning("The proposed null mean appears to be outside the convex hull for at least one sample (empirical likelihood weights are zero)", call. = FALSE)
    return(FALSE)
  }  
  lapply(wts, function(w){
    if (abs(length(w) - sum(w)) > 0.9){
      # above sees if weight sums to n (otherwise should sum to k < n being number of points in face). 
      warning(sprintf("Empirical likelihood weights sum to %0.1f, which suggests the proposed null mean is on a face of the convex hull.", sum(w)), call. = FALSE)
    }
  })
  return(TRUE)
}