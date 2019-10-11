#' Solve unidirectional constrained problem
#'
#' This function minimizes \eqn{\beta} in the 1D problem
#' \eqn{1/2 * ||y - x \beta||_2^2 + \lambda |\beta|} subject to \eqn{\beta < 0}.
#'
#' The analytical solution of this problem is given by
#' \deqn{\beta* = min(0, (y'x + \lambda) / x'x ).}
#'
#' @param y a vector of size n.
#' @param x a vector of size n.
#' @inheritParams estimate_shifts
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

#' @rdname solve_univariate
#'
#' @param beta the current value of beta.
#' @param coord the coordinate to be updated.
#' @param X a matrix of size n x p.
#'
#' @return the new value of beta
#' @export
update_univariate <- function(beta, coord, y, X, lambda){
  beta_i <- beta[-coord]
  X_i <- X[, -coord]
  xi <- X[, coord]
  y_i <- y - X_i %*% beta_i

  beta_next <- beta
  beta_next[coord] <- solve_univariate(y_i, xi, lambda)

  beta_next
}

#' @rdname solve_univariate
#'
#' @param beta0 the initial position of beta.
#'
#' @return the estimate value of beta
#' @export
solve_multivariate <- function(beta0, y, X, lambda) {
  p <- ncol(X)
  beta <- beta0
  for(i in 1:150){
    coord <- sample(p, size = 1)
    beta <- update_univariate(beta, coord, y, X, lambda)
  }
  return(beta)
}
