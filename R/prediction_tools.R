coef_reshape <- function(coeffs){
  coef_names <- names(coeffs)
  inter_pos <- grepl(":", coef_names)
  if(any(inter_pos)){
    left_vars <- lapply(coef_names[inter_pos], function(x){div_inter(x)$left})
    left_vars <- unlist(left_vars)
    right_vars <- lapply(coef_names[inter_pos], function(x){div_inter(x)$right})
    right_vars <- unlist(right_vars)
    #extract first character
    fchar_left <- substr(left_vars, 1, 1)
    fchar_right <- substr(right_vars, 1, 1)
    xz_pos <- (fchar_right == "x" & fchar_left == "z") |
              (fchar_right == "z" & fchar_left == "x")


  }
}
