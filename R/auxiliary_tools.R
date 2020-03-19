physical_dist <- function(x.loc, z.loc){
  x.long <- x.loc[1]
  x.lat <- x.loc[2]
  z.long <- z.loc[1]
  z.lat <- z.loc[2]
  tmp <- ((z.lat-x.lat)/2 * pi/180)^2
  tmp <- tmp + cos(z.lat * pi/180) * cos(x.lat * pi/180) * sin((z.long-x.long)/2 * pi/180)^2
  tmp <- 2 * atan2(tmp^0.5, (1-tmp)^0.5)
  return(tmp * 3959) #distance in miles
}

calc_dist <- function(pairwise){
  if(!pairwise){
    distance <- distHaversine(X.location, Z.location)
  } else {
    distance <- distm(X.location, Z.location, fun = distHaversine)
  }
  return(distance)
}

prepare_data <- function(n, wgt = rep(1, n), sim.factor = 1){
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
