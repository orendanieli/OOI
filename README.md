
OOI: Outside Option Index
=========================

An R package that calculates the Outside Option Index proposed by Caldwell and Danieli (2018). This index uses the cross- sectional concentration of similar workers across job types to quantify the availability of outside options as a function of workers’ characteristics (e.g. commuting costs, preferences, and skills.)

Installation
------------

Currently, only the GitHub version is available through:

``` r
#install.packages("devtools")
devtools::install_github("eladg9/OOI")
```

Usage Example
-------------

``` r
library(OOI)
#generate data
#worker and job characteristics:
n <- 1000
men <- rbinom(n, 1, 0.5)
size <- 1 + rgeom(n, 0.1) 
size[men == 0] <- size[men == 0] + 2
worker_resid <- matrix(round(runif(n, 0, 20), 1), ncol = 1)
job_location <- matrix(round(runif(n, 20, 40), 1), ncol = 1)
#prepare data
#define distance function:
dist_metric <- function(x, y){abs(y - x)}
X <- data.frame(men = men)
Z <- data.frame(size = size)
#add "x" / "z" to column names:
X <- add_prefix(X, "x.")
Z <- add_prefix(Z, "z.")
ooi_object <- OOI(formula = ~ x_*z_ + x_*d + z_*d, X = X, Z = Z,
                  X.location = worker_resid, Z.location = job_location,
                  sim.factor = 3, dist.fun = dist_metric, dist.order = 3)
#we can extract the ooi using predict():
ooi <- predict(ooi_object)
#we can also estimate the ooi with different coefficients (for counterfactual analysis):
coeffs <- ooi_object$coeffs
coeffs[names(coeffs) == "x.men"] <- 0
new_ooi <- predict(ooi_object, new.coef = coeffs)
#and there is also an option to predict the ooi for new data. see ?predict.ooi
#for further details.
```

References
----------
