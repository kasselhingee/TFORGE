testthat("stat_commoneigenvals() doesn't reject for simulation of single sample from null", {
  Ysample <- rsymm(50, 3)
  Ysample <- lapply(Ysample, `+`, diag(c(3,2,1)))
  

})
