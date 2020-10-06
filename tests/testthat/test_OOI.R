library(OOI)

set.seed(123)

#simulate data (for which we can calculate things theoretically)
n <- 2000
men_rate <- 0.1
women_rate <- 0.2
men <- rbinom(n, 1, 0.5)
X_loc <- matrix(runif(n, 0, 100), ncol = 1)
dist <- rep(NA, n)
men_inc <- men == 1
dist[men_inc] <- rexp(n = sum(men_inc), rate = men_rate) #distance for men
dist[!men_inc] <- rexp(n = sum(!men_inc), rate = women_rate) #distance for women
direction <- sample(c(1, -1), size = n, T)
Z_loc <- matrix(rep(NA, n), ncol = 1)
Z_loc[men_inc,] <- X_loc[men_inc,] + dist[men_inc] * direction[men_inc]
Z_loc[!men_inc,] <- X_loc[!men_inc,] + dist[!men_inc] * direction[!men_inc]
X <- matrix(men, ncol = 1, dimnames = list(NULL, "x.men"))
#define simple distance function
dis_function <- function(x, y){abs(x - y)}
ooi_obj <- suppressWarnings(OOI(~ x_ * d, X = X, X.location = X_loc, Z.location = Z_loc,
                                dist.fun = dis_function, dist.order = 1, sim.factor = 1))

test_that("OOI returns correct output", {
  ooi <- ooi_obj$ooi
  #choose workers who are far enough from the edges
  q25 <- quantile(X_loc[,1], probs = 0.3)
  q75 <- quantile(X_loc[,1], probs = 0.7)
  central <- (X_loc[,1] > q25) & (X_loc[,1] < q75)
  ooi_men <- ooi[men_inc & central]
  ooi_women <- ooi[!men_inc & central]
  #theoretical results
  theo_ooi_men <- 1 - log(50 * men_rate)
  theo_ooi_women <- 1 - log(50 * women_rate)
  ooi_men_hat <- mean(ooi_men)
  ooi_women_hat <- mean(ooi_women)
  expect_true(abs(ooi_men_hat - theo_ooi_men) < 0.1 &
                abs(ooi_women_hat - theo_ooi_women) < 0.1)
})

test_that("OOI returns the right job_worker probabilities", {
  p_hat <- ooi_obj$job_worker_prob
  p_theo <- 0.5 * men_rate * exp(-dist * men_rate) #men
  p_theo[!men_inc] <- 0.5 * women_rate * exp(-dist[!men_inc] * women_rate) #women
  diff <- abs(p_hat - p_theo)
  expect_true(max(diff) < 0.1)
})

test_that("OOI returns the same results for matrices and data frames with factors", {
  #simulate data (matrices and data frames)
  n <- 50
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
  mat_results <- suppressWarnings(OOI(~  x_*z_ + z_*d + x_*d, X = X_mat,
                                      Z = Z_mat, X_loc, Z_loc, sim.factor = 3,
                                      seed = 2))
  df_results <- suppressWarnings(OOI(~  x_*z_ + z_*d + x_*d, X = X_df,
                                      Z = Z_df, X_loc, Z_loc, sim.factor = 3,
                                     seed = 2))
  expect_true(max(abs(mat_results$ooi - df_results$ooi)) < 0.01)
})


