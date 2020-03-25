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
    est_data$d <- calc_dist(x_loc, z_loc, dist.fun)
  }
  #add high order distance
  if(dist.order > 1){
    for(i in 2:dist.order){
      est_data$tmp <- est_data$D^i
      colnames(est_data)[names(est_data) == "tmp"] <- paste("d", i, sep = "")
    }
  }
  #merge with X Z & weights
  est_data$w[est_data$y == 1] <- wgt
  est_data$w[est_data$y == 0] <- mean(wgt) #weights for fake matches
  est_data <- cbind(est_data, X[est_data$worker_id,], Z[est_data$job_id,])
  #prepare formula for estimation
  var_names <- c(colnames(X), colnames(Z))
  formula <- prep_form(formula, var_names, dist.order)
  #estimate logit
  logit <- glm(as.formula(formula), family = binomial(link='logit'),
               data = est_data, weights = est_data$w)
  coeffs <- logit$coefficients
  coeffs_sd <- sqrt(diag(vcov(logit)))
  pseudo_r2 <- round(1 - logit$deviance / logit$null.deviance, 3)
  #calculate standardized coefficients
  stand_coeffs <- standardize(coeffs[-1], est_data, wgt)

  output <- list(coeffs, coeffs_sd, pseudo_r2, stand_coeffs)
  return(output)
}
