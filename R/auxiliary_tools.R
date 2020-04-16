#calculates geo distance between *two* points
geo_dist <- function(x.loc, z.loc){
  x.long <- x.loc[1]
  x.lat <- x.loc[2]
  z.long <- z.loc[1]
  z.lat <- z.loc[2]
  tmp <- ((z.lat-x.lat)/2 * pi/180)^2
  tmp <- tmp + cos(z.lat * pi/180) * cos(x.lat * pi/180) * sin((z.long-x.long)/2 * pi/180)^2
  tmp <- 2 * atan2(tmp^0.5, (1-tmp)^0.5)
  return(tmp * 3959) #distance in miles
}

#calculates 1:1 distance
calc_dist <- function(X.location, Z.location, fun = geo_dist, dist.order = 2){
  nr_X <- nrow(X.location)
  nr_Z <- nrow(Z.location)
  if(nr_X != nr_Z){
    stop("X.location and Z.location should have same number of rows")
  } else {
    distance <- rep(NA, nr_X)
    for(i in 1:nr_X){
      distance[i] <- fun(X.location[i, ], Z.location[i, ])
    }
  }
  distance <- data.frame(d = distance)
  #add high order distance
  if(dist.order > 1){
    for(i in 2:dist.order){
      distance$tmp <- distance$d^i
      colnames(distance)[names(distance) == "tmp"] <- paste("d", i, sep = "")
    }
  }
  return(distance)
}

#adds simulated matches to real data
prep_data <- function(X, Z = NULL, wgt = rep(1, nrow(X)),
                      sim.factor = 1, seed = NULL){
  set.seed(seed)
  n <- nrow(X)
  n.fake = round(n * sim.factor)
  real_data <- data.frame(worker_id = 1:n, job_id = 1:n, y = rep(1, n))
  #simulate matches
  worker_fake_id <- sample(n, size = n.fake, replace = T, prob = wgt)
  job_fake_id <- sample(n, size = n.fake, replace = T, prob = wgt)
  fake_data <- data.frame(worker_id = worker_fake_id,
                          job_id = job_fake_id, y = rep(0, n.fake))
  res <- rbind(real_data, fake_data)
  res$y <- as.factor(res$y)
  #merge with X Z & weights
  res$w[res$y == 1] <- wgt
  res$w[res$y == 0] <- mean(wgt) #weights for fake matches
  if(is.null(Z)){
    res <- cbind(res, X[res$worker_id, , drop = F])
  } else {
    res <- cbind(res, X[res$worker_id, , drop = F], Z[res$job_id, , drop = F])
  }
  return(res)
}

#prepares formula for glm
#this function also validates that without distance, X*Z must be included
prep_form <- function(formula, var.names, dist.order = 2){
  terms_labels <- labels(terms(formula))
  n <- length(terms_labels)
  #flag for distance
  d_included <- F
  #flag for X*Z
  xz_included <- F
  #initialize formula
  form <- "~ 1"
  for(i in 1:n){
    term <- terms_labels[i]
    #Is this an interaction term?
    inter_term <- grepl(":", term)
    if(!inter_term){
      if(term == "d"){d_included <- T}
      form <- paste(form, ext_names(term, var.names, dist.order), sep = " + ")
    } else {
      #divide to left vars & right vars
      left_vars <- ext_names(div_inter(term)$left, var.names, dist.order)
      right_vars <- ext_names(div_inter(term)$right, var.names, dist.order)
      fchar_l <- substr(left_vars, 1, 1)
      fchar_r <- substr(right_vars, 1, 1)
      if(fchar_r == "x" & fchar_l == "z" | fchar_r == "z" & fchar_l == "x")
        xz_included <- T
      #paste back
      tmp <- paste0("(", left_vars, ")", ":", "(", right_vars, ")")
      form <- paste(form, tmp, sep = " + ")
    }
  }
  if(!d_included & !xz_included)
    stop("either distance ('d') or interaction between X & Z must be included")
  form <- paste("y", form, sep = " ")
  return(form)
}

#divides interaction term into left term and right term
div_inter <- function(inter_term){
  point_pos <- gregexpr(pattern = ":", inter_term)[[1]][1]
  res <- list(left = substr(inter_term, 1, point_pos - 1),
              right = substr(inter_term, point_pos + 1, nchar(inter_term)))
  return(res)
}


#extracts original variable names from term
ext_names <- function(term, var.names, dist.order = 2){
  #this term ends with "_"?
  n.char <- nchar(term)
  last_char <- substr(term, n.char, n.char)
  if(last_char == "_"){
    #take all variables starting with the expression to the left
    term_init <-  substr(term, 1, n.char - 1)
    names_init <- substr(var.names, 1, n.char - 1)
    res <- var.names[term_init == names_init]
    #Is this a distance term?
  } else if (term == "d"){
    dist_terms <- rep(NA, dist.order)
    for(j in 1:dist.order){
      dist_terms[j] <- ifelse(j == 1, "d", paste("d", j, sep = ""))
    }
    res <- dist_terms
  } else {
    res <- term
  }
  return(paste(res, collapse = " + "))
}


#standardizes coefficients
standardize <- function(coeffs, dat, wgt){
  #calculate sd for relevant variables
  coef_names <- names(coeffs)
  inter_pos <- grepl(":", coef_names)
  rel_vars <- coef_names[!inter_pos] #variables without interaction
  sd <- apply(dat[, rel_vars], 2,
              function(x, w = wgt){sqrt(Hmisc::wtd.var(x, w))})
  coeffs[rel_vars] <- coeffs[rel_vars] * sd
  #for interaction terms, we need to mulitply by the SD of each variable
  inter_pos <- which(inter_pos)
  for(i in inter_pos){
    term <- names(coeffs[i])
    var1 <- div_inter(term)$left
    var2 <- div_inter(term)$right
    coeffs[i] <- coeffs[i] * sd[var1] * sd[var2]
  }
  return(coeffs)
}

#converts data.frame to matrix and expands factors to a set of dummy variables
expand_matrix <- function(df){
  if(is.matrix(df) | is.null(df)){
    return(df)
  }
  form <- paste("~", paste0(colnames(df), collapse = "+"))
  form <- as.formula(form)
  df <- model.matrix(form, df)
  #delete intercept
  df <- df[, -1, drop = F]
  return(df)
}
