# Prepare the Gonjo Basin Data

# Download the data from https://doi.org/10.5281/zenodo.3666760
# This data has a creative commons 4.0 license
download.file(
  paste0(
    "https://zenodo.org/records/3666760/files/",
    "AMS%20data%20of%20the%20Gonjo%20Basin.xlsx?download=1"
  ),
  "gonjo.xlsx"
)

# Read data into R
ams <- readxl::read_xlsx("gonjo.xlsx", skip = 1, .name_repair = "minimal")
names(ams)[14:19] <- paste0(names(ams)[14:19], "I")
names(ams)[20:25] <- paste0(names(ams)[20:25], "T")

# There are two sets of columns starting with 'K', referring to in-situ (I) and tilt-corrected values (T) respectively.
# The AMS tensor eigenvectors and eigenvalues are stored separately.
# Below we extract first the eigenvalues, then the eigenvectors, then combine both into a symmetric tensor.

# Eigenvalues:
Tau3 <- 1 / (ams[, "P"] + ams[, "F"] + 1)
Tau2 <- 1 / (ams[, "L"] + 1 + (1 / ams[, "F"]))
Tau1 <- 1 - Tau2 - Tau3
evals <- data.frame(Tau1, Tau2, Tau3)

# First eigenvector:
inc_new <- ams[, "K1incI", drop = TRUE]
inc_new2 <- 90 - inc_new
theta <- inc_new2 * (pi / 180)
phi <- ams[, "K1decI", drop = TRUE] * (pi / 180)
y1 <- matrix(0, nrow(ams), 3)
y1[, 1] <- cos(theta)
y1[, 2] <- sin(theta) * cos(phi)
y1[, 3] <- sin(theta) * sin(phi)

#Second eigenvector:
inc_new <- ams[, "K2incI", drop = TRUE]
inc_new2 <- 90 - inc_new
theta <- inc_new2 * (pi / 180)
phi <- ams[, "K2decI", drop = TRUE] * (pi / 180)
y2 <- matrix(0, nrow(ams), 3)
y2[, 1] <- cos(theta)
y2[, 2] <- sin(theta) * cos(phi)
y2[, 3] <- sin(theta) * sin(phi)

#Third eigenvector:
inc_new <- ams[, "K3incI", drop = TRUE]
inc_new2 <- 90 - inc_new
theta <- inc_new2 * (pi / 180)
phi <- ams[, "K3decI", drop = TRUE] * (pi / 180)
y3 <- matrix(0, nrow(ams), 3)
y3[, 1] <- cos(theta)
y3[, 2] <- sin(theta) * cos(phi)
y3[, 3] <- sin(theta) * sin(phi)

#Combine into a list of symmetric matrices:
AMSmatrices <- lapply(1:nrow(ams), function(i) {
  evecs <- cbind(y1[i, ], y2[i, ], y3[i, ])
  eval <- unlist(evals[i, ])
  m <- evecs %*% diag(eval) %*% t(evecs)
  # make exactly symmetric - removing floating point asymmetries
  m[lower.tri(m)] <- (m[lower.tri(m)] + t(m)[lower.tri(m)]) / 2
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
})

# Save into the package
Gonjo <- list(
  matrices = AMSmatrices,
  datatable = ams
)
usethis::use_data(Gonjo, overwrite = TRUE)
