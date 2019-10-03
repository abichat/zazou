#' Compute objective function
#'
#' @param Y a vector of size m
#' @param X a vector of size m*(n+m)
#' @param lambda a positive regularization parameter
#'
#' @return a function that take Delta (a vector of size n+m) as argument
#' and returns a scalar
#' @export
compute_objective_function <- function(Y, X, lambda) {
  objective_function <- function(Delta) {
    YXD <- Y - X %*% Delta
    shrinkage <- lambda * sum(Delta)
    crossprod(YXD) - shrinkage
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
