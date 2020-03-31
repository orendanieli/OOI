#prepares coeffs for prediction process
coef_reshape <- function(coeffs){
  coef_names <- names(coeffs)
  inter_coef <- coef_names[grepl(":", coef_names)]
  if(length(inter_coef) != 0){
    #prepare coef matrix to calculate X*A*Z
    #rows are for x, columns for z
    xz_mat <- prep_matrix(inter_coef, "x", "z", coeffs)
    #prepare coef matrix to calculate X*A*D
    #rows are for x, columns for d
    xd_mat <- prep_matrix(inter_coef, "x", "d", coeffs)
    #prepare coef matrix to calculate D*A*Z
    #rows are for d, columns for z
    dz_mat <- prep_matrix(inter_coef, "d", "z", coeffs)
  }
  sim_coef <- coef_names[!grepl(":", coef_names)]
  fchar <-  substr(sim_coef, 1, 1)
  z_coef <- coeffs[sim_coef[fchar == "z"]]
  if (length(z_coef) == 0)
    z_coef <- NULL
  d_coef <- coeffs[sim_coef[fchar == "d"]]
  if (length(d_coef) == 0)
    d_coef <- NULL
  return(list(xz_mat = xz_mat, xd_mat = xd_mat, dz_mat = dz_mat,
              z_coef = z_coef, d_coef = d_coef))
}

#prepares matrix of coefficients for interaction terms
prep_matrix <- function(inter_coef, term1, term2, coeffs){
  left_vars <- unlist(lapply(inter_coef, function(x){div_inter(x)$left}))
  right_vars <- unlist(lapply(inter_coef, function(x){div_inter(x)$right}))
  #extract first character
  fchar_left <- substr(left_vars, 1, 1)
  fchar_right <- substr(right_vars, 1, 1)
  #position of relevant interactions
  pos <- (fchar_right == term1 & fchar_left == term2) |
          (fchar_right == term2 & fchar_left == term1)
  if(any(pos)){
    vars1 <- unique(c(left_vars[pos & fchar_left == term1],
                      right_vars[pos & fchar_right == term1]))
    vars2 <- unique(c(left_vars[pos & fchar_left == term2],
                      right_vars[pos & fchar_right == term2]))
    #rows are for term1, columns for term2
    mat <- matrix(, nrow = length(vars1), ncol = length(vars2),
                dimnames = list(vars1, vars2))
    #fill mat with coeffs
    for(r in 1:nrow(mat)){
      for(c in 1:ncol(mat)){
        ind <- (left_vars == vars1[r] & right_vars == vars2[c]) |
          (left_vars == vars2[c] & right_vars == vars1[r])
        mat[r, c] <- coeffs[inter_coef[ind]]
      }
    }
    return(mat)
  } else {
    return(NULL)
  }
}

#calculates ooi for the logit method, by looping over workers district
#(district is just grouping of X.location). This is nice because d(i,j) is the
#same for all workers from the same district.
calc_ooi <- function(coef.mat, X, Z, X.location = NULL,
                     Z.location = NULL, wgt = rep(1, nrow(X))){
  #wgt <- wgt / sum(wgt)
  n <- nrow(X)
  A <- coef.mat$xz_mat
  one <- rep(1, n) #will be useful later
  b <- coef.mat$z_coef
  Xnames <- rownames(A)
  Znames <- colnames(A)
  Zcons_names <- names(b)
  X <- X[,Xnames]
  #add column of 1 for Z_cons
  X <- cbind(X, one)
  Z_cons <- b %*% t(Z[,Zcons_names])
  AZ <- A %*% t(Z[,Znames])
  AZ <- rbind(AZ, Z_cons)
  if(is.null(X.location)){
    #generate fake district table
    n_dis <- max(round(n / 50), 2)
    dis <- cut(1:n, n_dis, labels = as.character(1:n_dis))
    dis_table <- cbind.data.frame(worker = 1:n, dis = dis, ooi = rep(NA, n))
    for(i in 1:n_dis){
      workers <- dis_table$worker[dis_table$dis == i]
      chunk <- X[workers,] %*% AZ
      chunk <- chunk #- chunk[,1] - 13
      p <- exp(chunk)
      p <- t(t(p) * wgt)
      #normalize (sum(p(.,j)) = 1)
      sum_p <- as.vector(p %*% one)
      p <- p / sum_p
      chunk <- chunk - log(sum_p) #chunk is in log units, so its just like dividing
      #Sum and get index
      ooi <- -(p * log(p)) %*% one
      dis_table$ooi[dis_table$dis == i] <- ooi
    }
  } else {
    A_dist <- coef.mat$xd_mat
    Xnames <- row.names(A_dist)

  }
  return(dis_table$ooi)
}
