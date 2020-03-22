X <- matrix(rnorm(80), ncol = 4)
Z <- matrix(rnorm(80), ncol = 4)
X_loc <- matrix(exp(40), ncol = 2)
Z_loc <- matrix(exp(40), ncol = 2)

form = ~ x_ * z_ + x_ * d
