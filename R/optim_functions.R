#' Compute objective function
#'
#' @param Y a vector of size m
#' @param X a vector of size m*(n+m)
#' @param lambda a positive regularization parameter
#'
#' @return a function that take Delta, a vector of size n+m, as argument
#' @export
#'
#' @examples
compute_objective_function <- function(Y, X, lambda) {
  objective_function <- function(Delta) {
    YXD <- Y - X %*% Delta
    shrinkage <- lambda * sum(Delta)
    crossprod(YXD) - shrinkage
  }
  return(objective_function)
}
