n = 1000
X <- matrix(rnorm(4 * n), ncol = 4,
            dimnames = list(NULL, c("x.1", "x.2", "x.3", "x.4")))
gender <- rbinom(n , 1, 0.5)
gender <- factor(gender, levels = c(0,1), labels = c("men", "women"))
native <- rbinom(n , 1, 0.5)
native <- factor(native, levels = c(0,1), labels = c("yes", "no"))
#X <- data.frame(x.gender = gender, x.native = native)
Z <- matrix(rnorm(4 * n), ncol = 4,
            dimnames = list(NULL, c("z.1", "z.2", "z.3", "z.4")))
X_loc <- matrix(runif(2 * n, 40, 42), ncol = 2)
Z_loc <- matrix(runif(2 * n, 40, 42), ncol = 2)
w = rexp(n)

bla <- suppressMessages(OOI(~ x_ * z_ , X = X, Z = Z, sim.factor = 4))#,
           X.location = X_loc, Z.location = Z_loc)


