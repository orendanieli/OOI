#' Outside Option Index
#'
#' calculates the 'outside option index' (defined as
#' \eqn{-\sum P(Z|X) * log(P(Z|X) / P(Z))})
#'  for workers, using employer-employee data.
#'
#' @param formula a formula describing the model to be fitted in order to
#'                estimate P(Z|X) / P(Z). This formula uses a syntax similar to STATA, and so
#'                "x_" refers to all variables with the prefix "x", while
#'                "z_" refers to all variables with the prefix "z". Similarly, "d" refers
#'                to the distance polynomial (see the example below).
#' @param X matrix or data frame with workers characteristics. Note that all column names
#'          should start with "x" (necessary for \code{\link{coef_reshape}}).
#' @param Z an optional matrix or data frame with jobs characteristics. Note that all column names
#'          should start with "z" (necessary for \code{\link{coef_reshape}}).
#' @param X.location an optional matrix or data frame with location for workers. could be
#'                   geographical location (i.e., geo-coordinates) or any other
#'                   feature that can be used in order to measure distance between
#'                   worker and job using 'dist.fun'.
#' @param Z.location same as 'X.location' but for jobs.
#' @param wgt an optional numeric vector of weights.
#' @param pred logical. If TRUE (default), predicts the ooi for the provided data.
#' @param method a method for estimating P(Z|X) / P(Z). currently not in use.
#' @param sim.factor a variable that determines how much fake data to simulate
#'                   (relative to real data).
#' @param dist.fun a distance function to calculate the distance between X.location and
#'                 Z.location. The default function \code{\link{geo_dist}} is suitable
#'                 for data with geo-coordinates.
#' @param dist.order the order of the distance polynomial.
#' @param seed the seed of the random number generator.
#'
#' @return an object of class "ooi". This object is a list containing
#' the following components:
#'  \item{coeffs}{coefficients from the estimated logit.}
#'  \item{coeffs_sd}{coefficients SE.}
#'  \item{pseudo_r2}{McFadden's pseudo-R squared for the estimated logit.}
#'  \item{standardized_coeffs}{standardized coefficients.}
#'  \item{ooi}{the outside option index.}
#'  \item{orig_arg}{a list containing the original arguments (necessary
#'  for \code{\link{predict.ooi}}).}

#' @export

OOI <- function(formula = NULL,
                X,
                Z = NULL,
                X.location = NULL,
                Z.location = NULL,
                wgt = rep(1, nrow(X)),
                pred = T,
                method = "logit",
                sim.factor = 1,
                dist.fun = geo_dist,
                dist.order = 2,
                seed = runif(1, 0, .Machine$integer.max)){
  formula <- validate_input(formula, X, Z, X.location, Z.location, wgt)
  #prepare formula for estimation
  var_names <- c(colnames(X), colnames(Z))
  formula <- prep_form(formula, var_names, dist.order)
  #prepare data for logit estimation
  est_data <- prep_data(X, Z, wgt, sim.factor, seed)
  #merge with distance
  if(!is.null(X.location)){
    if(all(apply(X.location, 2, is.numeric))){
      #rounding is necessary for gen_dist()
      X.location = round(X.location, 3)
    }
    x_loc <- X.location[est_data$worker_id, , drop = F]
    z_loc <- Z.location[est_data$job_id, , drop = F]
    D <- calc_dist(x_loc, z_loc, dist.fun, dist.order)
    est_data <- cbind(est_data, D)
  }
  #estimate logit
  logit <- glm(as.formula(formula), family = binomial(link='logit'),
               data = est_data, weights = est_data$w)
  coeffs <- logit$coefficients
  if(logit$rank < length(coeffs)){
    warning(paste("Logit was estimated on a singular matrix.",
                  "Replacing NA coefficients with 0"))
    coeffs <- replace(coeffs, is.na(coeffs), 0)
  }
  coeffs_se <- round(sqrt(diag(vcov(logit))), 4)
  pseudo_r2 <- round(1 - (logit$deviance / logit$null.deviance), 3)
  #calculate standardized coefficients
  est_data <- expand_matrix(est_data[est_data$y == 1,])
  stand_coeffs <- standardize(coeffs[-1], est_data, wgt)
  #prepare output
  orig_arg <- list(X = X,
                    Z = Z,
                    X.location = X.location,
                    Z.location = Z.location,
                    wgt = wgt,
                    dist.fun = dist.fun,
                    dist.order = dist.order)
  output <- list(coeffs = coeffs,
                 coeffs_se = coeffs_se,
                 pseudo_r2 = pseudo_r2,
                 standardized_coeffs = stand_coeffs,
                 orig_arg = orig_arg,
                 formula = formula)
  #predict OOI
  if(pred){
    #reshape coefficients (necessary for prediction)
    coef_matrices <- coef_reshape(coeffs)
    ooi <- predict_ooi(coef_matrices, X, Z, X.location, Z.location,
                       wgt, dist.fun, dist.order, F)
    output[["ooi"]] <- ooi
  }
  class(output) <- "ooi"
  return(output)
}


#' Predict Outside Option Index
#'
#' predicts the OOI for new coefficients (for counterfactual analysis)
#' and/or new data. In addition, the HHI (Herfindahl-Hirschman Index) can be predicted instead
#' of the OOI.
#'
#' @param object an ooi object.
#' @param new.coef a new *named* vector of coefficients. check the coefficients produced by
#'                 the main function to see the right format for this vector.
#' @param new.X a new X matrix / data frame.
#' @param new.Z a new Z matrix / data frame.
#' @param new.X.location a new X.location matrix / data frame.
#' @param new.Z.location a new Z.location matrix / data frame.
#' @param new.wgt a new vector of weights
#' @param hhi whether to calculate the HHI (Herfindahl-Hirschman Index) instead of the OOI.
#'            default is FALSE.
#'
#' @return if there are no new arguments and hhi = FALSE, returns the original ooi. otherwise,
#' returns a vector of ooi/hhi calculated using the new/original arguments.
#' @export
#add example
predict.ooi <- function(object,
                        new.coef = NULL,
                        new.X = NULL,
                        new.Z = NULL,
                        new.X.location = NULL,
                        new.Z.location = NULL,
                        new.wgt = NULL,
                        hhi = F){
  #if there are no new inputs, return pre-computed ooi.
  if(is.null(new.coef) & is.null(new.X) & is.null(new.Z) & !hhi){
    return(object$ooi)
  } else {
    #prepare inputs
    x <- object$orig_arg
    coeffs <- if(is.null(new.coef)) object$coeffs else new.coef
    X <- if(is.null(new.X)) x$X else new.X
    Z <- if(is.null(new.Z)) x$Z else new.Z
    X.location <- if(is.null(new.X.location)) x$X.location else new.X.location
    Z.location <- if(is.null(new.Z.location)) x$Z.location else new.Z.location
    wgt <- if(is.null(new.wgt)) x$wgt else new.wgt
    #reshape coefficients (necessary for prediction)
    coef_matrices <- coef_reshape(coeffs)
    #predict OOI
    ooi <- predict_ooi(coef_matrices, X, Z, X.location, Z.location,
                       wgt, x$dist.fun, x$dist.order, hhi)
    return(ooi)
  }
}

#' Add prefix
#'
#' Adds a prefix to the column names of a matrix
#'
#' @param df data.frame or matrix
#' @param prefix a prefix to be added
#'
#' @return matrix / data.frame with new column names.
#' @export
add_prefix <- function(df, prefix){
  col_names <- colnames(df)
  col_names <- paste0(prefix, col_names)
  colnames(df) <- col_names
  return(df)
}
