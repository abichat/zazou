#' Compute objective function
#'
#' @param Y A vector of size m.
#' @param X A vector of size m*(n+m).
#' @param sigma Needed if \code{type} is set to \code{"scaledlasso"}.
#' @param type Character. If you want the \code{"lasso"} or
#' \code{"scaledlasso"} objective function.
#' @inheritParams estimate_shifts
#'
#' @return a function that take Delta (a vector of size n+m) as argument
#' and returns a scalar
#' @export
compute_objective_function <- function(Y, X, lambda, sigma,
                                       type = c("lasso", "scaledlasso")) {
  objfun <- switch (type,
                    "lasso" = compute_OF_lasso(Y, X, lambda),
                    "scaledlasso" = compute_OF_scaledlasso(Y, X, lambda, sigma)
  )

  return(objfun)
}

compute_OF_lasso <- function(Y, X, lambda){
  objective_function <- function(Delta) {
    YXD <- Y - X %*% Delta
    shrinkage <- lambda * sum(abs(Delta))
    (sum(YXD^2) / 2 + shrinkage)
  }
  return(objective_function)
}

compute_OF_scaledlasso <- function(Y, X, lambda, sigma){
  objective_function <- function(Delta) {
    YXD <- Y - X %*% Delta
    shrinkage <- lambda * sum(abs(Delta))
    (sum(YXD^2) / 2 / sigma / length(Y) + sigma / 2 + shrinkage)
  }
  return(objective_function)
}


#' Compute gradient function
#'
#' @inheritParams compute_objective_function
#'
#' @return a function that take Delta (a vector of size n+m) as argument
#' and returns a vector of size n+m
#' @export
compute_gradient_function <- function(Y, X, lambda){
  gradient_function <- function(Delta){
    -crossprod(X, Y) + crossprod(X) %*% Delta - rep(lambda, length(Delta))
  }
  return(gradient_function)
}
