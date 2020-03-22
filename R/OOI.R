#' Outside Option Index
#'
#' calculates OOI
#'
#' @param formula
#' @param X matrix or data frame with workers characteristics.
#' @param Z matrix or data frame with jobs characteristics.
#' @param X.location an optional matrix or data frame with location for workers. could be
#'                   geographical location (i.e., geo-coordinates) or any other
#'                   feature that can be used in order to measure distance between
#'                   worker and job using 'dist.fun'.
#' @param Z.location same as 'X.location' but for jobs.
#' @param wgt an optional vector of weights.
#' @param method
#' @param sim.factor
#' @param dist.fun
#' @param dist.order

OOI <- function(formula = NULL,
                X,
                Z,
                X.location = NULL,
                Z.location = NULL,
                wgt = rep(1, nrow(X)),
                method = "logit",
                sim.factor = 1,
                dist.fun = geo_dist,
                dist.order = 2){
  #validate input
  #validate_data
  #validate_formula
  #prepare data for logit estimation
  n <- nrow(X)
  est_data <- prep_data(n, wgt, sim.factor)
  #merge with distance
  if(!is.null(X.location)){
    x_loc <- X.location[est_data$worker_id,]
    z_loc <- Z.location[est_data$job_id,]
    est_data$distance <- calc_dist(x_loc, z_loc, dist.fun)
  }
  #merge with X Z & weights
  est_data <- cbind(est_data, w = wgt[est_data$job_id],
                    X[est_data$worker_id,], Z[est_data$job_id,])
  #prepare formula for estimation
  X_names <- colnames(X)
  Z_names <- colnames(Z)
  formula <- prep_form(formula, X_names, Z_names, dist.order)
  return(est_data)
}
