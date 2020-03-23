X <- matrix(rnorm(80), ncol = 4,
            dimnames = list(NULL, c("x.1", "x.2", "x.3", "x.4")))
Z <- matrix(rnorm(80), ncol = 4,
            dimnames = list(NULL, c("z.1", "z.2", "z.3", "z.4")))
X_loc <- matrix(rexp(40), ncol = 2)
Z_loc <- matrix(rexp(40), ncol = 2)
w = rexp(20)

bla <- OOI(~ x_ * z_ + x_ * D, X = X, Z = Z,
           X.location = X_loc, Z.location = Z_loc, wgt = w)


