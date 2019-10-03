#' Solve unidirectional constrained problem
#'
#' This function minimizes \eqn{\beta} in the 1D problem \eqn{1/2 * ||y - x \beta||_2^2 + \lambda |\beta|} subject to \eqn{\beta < 0}.
#'
#' The analytical solution of this problem is given by \deqn{\beta* = min(0, (y'x + \lambda) / x'x ).}
#'
#' @param y a vector of size n
#' @param x a vector of size n
#' @param lambda a positive regularization parameter
#'
#' @return The scalar solution of the 1D optimization problem
#' @export
#'
#' @examples
#' solve_univariate(1:4, -(4:1), 2)
solve_univariate <- function(y, x, lambda = 0) {
  unconstrained <- (crossprod(y, x) + lambda) / crossprod(x)
  min(0, unconstrained)
}
