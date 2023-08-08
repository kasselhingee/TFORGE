There are the following types of functions

+ statistics: usually z^2 values that are asymptotically chi-sq
  + that takes lists of tensors, or lists of lists for multisample
  + returns the value of the stat, __with other useful information as attributes__
  + stats for given values (single sample) accept given parameter values
    + naming convention: `stat_ss_??`
  + multisample stats use the single sample stats after estimating
    + naming convention: `stat_ms_??`
+ estimators: for the common value between multiple samples
+ resampling from NULL preparation methods: 
  + either transform the sample with replacement OR
  + find weights and sample according to weights
+ general bootstrappy resampling thing that takes:
  + a (multi)sample (the original data), a NULLified sample (or empty), weights (or empty), and a stat function, and applies stat to the first sample and all resamples. Possibly in parallel.
+ functions that combine others into user friendly:
  + tests
  + confidence regions
+ utility functions
+ detector of constraints

