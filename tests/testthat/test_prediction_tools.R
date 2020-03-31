library(OOI)

set.seed(123)

test_that("coef_reshape returns correct coefficient matrices", {
  coeffs <- round(rexp(12), 2)
  names(coeffs) <- c("x1", "x2", "z1", "z2", "d", "x1:z1",
                     "x1:z2", "x2:z1", "x2:z2", "d:x1", "d:x2", "d:z1")
  b <- coeffs[3:4]
  A <- matrix(coeffs[6:9], ncol = 2, byrow = T)
  A_dist <- coeffs[10:11]
  dz_coef <- coeffs[12]
  output <- coef_reshape(coeffs)
  expect_true(all(output$xz_mat == A) & all(output$z_coef == b) &
                all(output$xd_mat == A_dist) & all(output$dz_mat == dz_coef))
})

test_that("calc_ooi returns the same results as predict", {
  n <- 20
  p <- 4
  w <- rexp(n)
  Xnames <- paste0("x.", 1:p)
  Znames <- paste0("z.", 1:p)
  X <- matrix(rnorm(n*p), ncol = p, dimnames = list(NULL, Xnames))
  Z <- matrix(rnorm(n*p), ncol = p, dimnames = list(NULL, Znames))
  est_data <- prep_data(X, Z, wgt =  w, sim.factor = 4, seed = 4)
  #prepare formula for estimation
  formula <- prep_form(~ x_ * z_ , c(Xnames, Znames))
  #estimate logit
  logit <- glm(as.formula(formula), family = binomial(link='logit'),
               data = est_data, weights = est_data$w)
  coeffs <- logit$coefficients
  coef_matrices <- coef_reshape(coeffs)
  ooi_package <- calc_ooi(coef_matrices, X, Z, wgt = w)
  ooi_predict <- c()
  #calc ooi for each worker with predict
  for(i in 1:n){
    Xi <- matrix(rep(X[i,], n), ncol = p, byrow = T,
                 dimnames = list(NULL, Xnames))
    df <- cbind.data.frame(Xi, Z)
    f <- predict(logit, newdata = df)
    f <- exp(f) * w
    #normalize
    f <- f / sum(f)
    ooi_predict <- c(ooi_predict, -sum(f*log(f)))
  }
  expect_true(min(abs(ooi_predict - ooi_package)) < 0.001)
})


