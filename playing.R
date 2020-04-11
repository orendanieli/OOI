n = 100
X <- matrix(rnorm(4 * n), ncol = 4,
            dimnames = list(NULL, c("x.1", "x.2", "x.3", "x.4")))
Z <- matrix(rnorm(4 * n), ncol = 4,
            dimnames = list(NULL, c("z.1", "z.2", "z.3", "z.4")))
X_loc <- matrix(runif(2 * n, 40, 42), ncol = 2)
Z_loc <- matrix(runif(2 * n, 40, 42), ncol = 2)
w = rexp(n)

bla <- OOI(~ x_ * z_ + x_ * d, X = X, Z = Z, sim.factor = 4,
           X.location = X_loc, Z.location = Z_loc)

bla <- OOI(~ x_ * d, X = X, sim.factor = 4,
           X.location = X_loc, Z.location = Z_loc)


y <- c(rep(1,5), rep(0,5))
x <- c(1,5,4,4,2,6,6,4,3,5)
x2 <- c(2,3,3,1,1,2,1,2,8,6)

y = rnorm(10000)
x = rexp(10000)
reg = lm(y ~ x)

forest = causal_forest(X = X, Y = rnorm(n), W = rbinom(n, 1, 0.5))
