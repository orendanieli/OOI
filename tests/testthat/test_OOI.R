library(OOI)

set.seed(123)

test_that("OOI returns correct output", {
  #simulate data
  n <- 1000
  men <- rbinom(n, 1, 0.5)
  X_loc <- matrix(runif(n, 0, 100), ncol = 1)
  dist <- rep(NA, n)
  men_inc <- men == 1
  dist[men_inc] <- rexp(n = sum(men_inc), rate = 1) #distance for men
  dist[!men_inc] <- rexp(n = sum(!men_inc), rate = 2) #distance for women
  direction <- sample(c(1, -1), size = n, T)
  Z_loc <- matrix(rep(NA, n), ncol = 1)
  Z_loc[men_inc,] <- X_loc[men_inc,] + dist[men_inc] * direction[men_inc]
  Z_loc[!men_inc,] <- X_loc[!men_inc,] + dist[!men_inc] * direction[!men_inc]
  X <- matrix(men, ncol = 1, dimnames = list(NULL, "x.men"))
  #define simple distance function
  dis_function <- function(x, y){abs(x - y)}
  ooi_obj <- OOI(~ x_ * d, X = X, X.location = X_loc, Z.location = Z_loc,
                 dist.fun = dis_function, dist.order = 1, sim.factor = 3)
  ooi <- ooi_obj$ooi
  #choose workers who are far enough from the edges
  q25 <- quantile(X_loc[,1], probs = 0.25)
  q75 <- quantile(X_loc[,1], probs = 0.75)
  central <- (X_loc[,1] > q25) & (X_loc[,1] < q75)
  ooi_men <- ooi[men_inc & central]
  ooi_women <- ooi[!men_inc & central]
  #theoretical results
  theo_ooi_men <- 1 - log(50)
  theo_ooi_women <- 1 - log(100)
  ooi_men_hat <- mean(ooi_men)
  ooi_women_hat <- mean(ooi_women)
  expect_true(abs(ooi_men_hat - theo_ooi_men) < 0.1 &
                abs(ooi_women_hat - theo_ooi_women) < 0.1)
})
