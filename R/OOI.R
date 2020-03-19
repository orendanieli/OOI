OOI <- function(formula = NULL,
                X,
                Z,
                X.location = NULL,
                Z.location = NULL,
                wgt = rep(1, nrow(X)),
                method = "logit",
                sim.factor = 1,
                dist.func = "physical_dist",
                dist.order = 2){
  #validate input
  #validate_data
  #validate_formula
  #prepare data for logit estimation
  n <- nrow(X)
  est_data <- prepare_data(n, wgt, sim.factor)
  #merge with distance
  if(!is.null(X.location)){
    x_loc <- X.location[est_data$worker_id,]
    z_loc <- Z.location[est_data$job_id,]
    est_data$distance <- calc_dist(x_loc, z_loc)
  }
  #merge with X Z & weights
  est_data <- cbind(est_data, w = wgt[est_data$job_id],
                    X[est_data$worker_id,], Z[est_data$job_id,])
  return(est_data)
}
