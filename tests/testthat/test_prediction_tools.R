library(OOI)


test_that("coef_reshape returns correct coefficient matrices", {
  coeffs <- round(rexp(8), 2)
  names(coeffs) <- c("x1", "x2", "z1", "z2", "x1:z1", "x1:z2", "x2:z1", "x2:z2")
  b <- coeffs[c(3,4)]
  A <- matrix(coeffs[1:4], ncol = 2, byrow = T)
  output <- coef_reshape(coeffs)
  expect_true(all(output$xz_mat == A) & all(output$z_coef == b))
})
