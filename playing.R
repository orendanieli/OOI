n = 100
X <- matrix(rnorm(4 * n), ncol = 4,
            dimnames = list(NULL, c("x.1", "x.2", "x.3", "x.4")))
gender <- rbinom(n , 1, 0.5)
gender <- factor(gender, levels = c(0,1), labels = c("women", "men"))
native <- rbinom(n , 1, 0.5)
native <- factor(native, levels = c(0,1), labels = c("yes", "no"))
X <- data.frame(x.gender = gender, x.native = native)
######
Z <- matrix(rbinom(4 * n, 20, 0.3), ncol = 4,
            dimnames = list(NULL, c("z.11", "z.12", "z.21", "z.22")))
Z <- as.data.frame(Z)
Z$z.11 <- as.factor(Z$z.11)
Z$z.21 <- as.factor(Z$z.21)
####
X_loc <- matrix(runif(2 * n, 40, 42), ncol = 2)
Z_loc <- matrix(runif(2 * n, 40, 42), ncol = 2)
w = rexp(n)


bla <- OOI(~ x_ * z.1_ + z.2_ + x_*d + z.1_ * d, X = X, Z = Z,
           X.location = X_loc, Z.location = Z_loc, sim.factor = 1, dist.fun = two_dim_function,
           pred = F, dist.order = c(2, 2))
############################
dif = rep(NA,10)
for(j in 1:10){
  seed = runif(1, 0, .Machine$integer.max)
  set.seed(seed)
  print(seed)
  n <- 200
  men <-rbinom(n, 1, 0.5)
  native <- rbinom(n, 1, 0.5)
  size <- rbinom(n, 3, 0.5)
  wage <- rnorm(n, 100, 10)
  X_df <- data.frame(x.men = factor(men, c(0,1), c("yes", "no")),
                     x.native = factor(native, c(0,1), c("yes", "no")))
  X_mat <- as.matrix(cbind(men, native))
  X_mat <- add_prefix(X_mat, "x.")
  Z_df <- data.frame(z.size = factor(size, c(0,1,2,3), c("A", "B", "C", "D")),
                     z.wage = wage)
  Z_mat <- cbind(wage, A = 1*(size == 1), B = 1*(size == 0),
                 C = 1*(size == 3))
  Z_mat <- add_prefix(Z_mat, "z.")
  X_loc <- matrix(runif(2*n, 40, 42), ncol = 2)
  Z_loc <- matrix(runif(2*n, 40, 42), ncol = 2)
  #compare results:
  mat_results <- suppressWarnings(OOI(~  x_*z_ + x_*d + z_*d, X = X_mat,
                                      Z = Z_mat, X_loc, Z_loc, sim.factor = 3,
                                      seed = seed, dist.order = 3))
  df_results <- suppressWarnings(OOI(~  x_*z_ + x_*d + z_*d, X = X_df,
                                     Z = Z_df, X_loc, Z_loc, sim.factor = 3,
                                     seed = seed, dist.order = 3))
  print(summary(abs(mat_results$ooi - df_results$ooi)))
  dif[j] = max(abs(mat_results$ooi - df_results$ooi))
}
summary(dif)


two_dim_function = function(x, z){
  x = x[1]; z = z[1]
  c(abs(x-z), (x-z)^2)
}

calc_dist(X_loc[1:10,], Z_loc[1:10,], fun = two_dim_function, dist.order = c(2, 3))


X = cbind(X, rnorm(n))
colnames(X)[2] = "x.random"
ooi_obj <- suppressWarnings(OOI(~ x_ * d + x.random*x.men, X = X, X.location = X_loc, Z.location = Z_loc,
                                dist.fun = dis_function, dist.order = c(2,1), sim.factor = 1))



