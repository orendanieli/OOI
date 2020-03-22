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
calc_dist <- function(X.location, Z.location, fun = geo_dist){
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
  return(distance)
}

#adds simulated matches to real data
prep_data <- function(n, wgt = rep(1, n), sim.factor = 1){
  n.fake = round(n * sim.factor)
  real_data <- data.frame(worker_id = 1:n, job_id = 1:n, y = rep(1, n))
  #simulate matches
  worker_fake_id <- sample(n, size = n.fake, replace = T, prob = wgt)
  job_fake_id <- sample(n, size = n.fake, replace = T, prob = wgt)
  fake_data <- data.frame(worker_id = worker_fake_id,
                          job_id = job_fake_id, y = rep(0, n.fake))
  res <- rbind(real_data, fake_data)
  return(res)
}

#prepares formula for glm
prep_form <- function(formula, X.names, Z.names, dist.order = 2){
  terms_labels <- labels(terms(formula))
  n <- length(terms_labels)
  #initialize list
  new_labels <- vector(mode = "list", length = n)
  for(i in 1:n){
    term <- terms_labels[i]
    #Is this an interaction term?
    inter_term <- grepl(":", term)
    if(!inter_term){
      #this term ends with "_"?
      #(=take all variables starting with the expression to the left?)
      n.char <- nchar(term)
      last_char <- substr(term, n.char, n.char)
      if(last_char == "_"){
        new_labels[i] <- f(substr(term, 1, n.char - 1), c(X.names, Z.names))
      } else {

      }
    }
  }

}
