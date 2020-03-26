library(OOI)


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
