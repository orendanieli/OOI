#' Outside Option Index
#'
#' calculates OOI
#'
#' @param formula a formula describing the model to be fitted in order to
#'                estimate P(Z|X). This formula uses a syntax similar to STATA, and so
#'                "x_" refers to all variables with the prefix "x". Similarly, "d" refers
#'                to the distance polynomial (see the example below).
#' @param X matrix or data frame with workers characteristics.
#' @param Z matrix or data frame with jobs characteristics.
#' @param X.location an optional matrix or data frame with location for workers. could be
#'                   geographical location (i.e., geo-coordinates) or any other
#'                   feature that can be used in order to measure distance between
#'                   worker and job using 'dist.fun'.
#' @param Z.location same as 'X.location' but for jobs.
#' @param wgt an optional vector of weights.
#' @param method a method for estimating P(Z|X). currently not in use.
#' @param sim.factor a variable that determines how much fake data to simulate
#'                   (relative to real data).
#' @param dist.fun a distance function to calculate the distance between X.location and
#'                 Z.location.
#' @param dist.order the order of the distance polynomial.
#' @param seed the seed of the random number generator.

OOI <- function(formula = NULL,
                X,
                Z,
                X.location = NULL,
                Z.location = NULL,
                wgt = rep(1, nrow(X)),
                method = "logit",
                sim.factor = 1,
                dist.fun = geo_dist,
                dist.order = 2,
                seed = runif(1, 0, .Machine$integer.max)){
  #validate input
  #validate_data
  #validate_formula
  #prepare data for logit estimation
  est_data <- prep_data(X, Z, wgt, sim.factor, seed)
  #merge with distance
  if(!is.null(X.location)){
    x_loc <- X.location[est_data$worker_id,]
    z_loc <- Z.location[est_data$job_id,]
    D <- calc_dist(x_loc, z_loc, dist.fun, dist.order)
    est_data <- cbind(est_data, D)
  }
  #prepare formula for estimation
  var_names <- c(colnames(X), colnames(Z))
  formula <- prep_form(formula, var_names, dist.order)
  #estimate logit
  logit <- glm(as.formula(formula), family = binomial(link='logit'),
               data = est_data, weights = est_data$w)
  coeffs <- logit$coefficients
  coeffs_sd <- sqrt(diag(vcov(logit)))
  pseudo_r2 <- round(1 - (logit$deviance / logit$null.deviance), 3)
  #calculate standardized coefficients
  stand_coeffs <- standardize(coeffs[-1], est_data, wgt)
  #reshape coefficients (necessary for prediction)
  coef_matrices <- coef_reshape(coeffs)
  #predict OOI
  ooi <- predict_ooi(coef_matrices, X, Z, X.location, Z.location,
                     wgt, dist.fun, dist.order)
  output <- list(coeffs = coeffs,
                 coeffs_sd = coeffs_sd,
                 pseudo_r2 = pseudo_r2,
                 standardized_coeffs = stand_coeffs,
                 ooi = ooi)
  return(output)
}
