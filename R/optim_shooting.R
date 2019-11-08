#' Solve unidirectional constrained problem
#'
#' This function minimizes \eqn{\beta} in the 1D problem
#' \eqn{1/2 * ||y - x \beta||_2^2 + \lambda |\beta|} subject to \eqn{\beta <= 0} (or \eqn{x\beta <= 0} coordinate wise)
#'
#' The analytical solution of this problem is given by
#' \deqn{\beta* = min(0, (y'x + \lambda) / x'x ).}
#' for the first constraint and is slightly more complicated for the second constraint (refer to the corresponding vignette)
#'
#' @param y a vector of size n.
#' @param x a vector of size n.
#' @param allow_positive Logical. Default FALSE. Allow positive values for \eqn{\beta} (but still enforce the constraint \eqn{x\beta <= 0})
#' @inheritParams estimate_shifts
#'
#' @return The scalar solution of the 1D optimization problem
#' @export
#'
#' @examples
#' solve_univariate(1:4, -(4:1), 2)
solve_univariate <- function(y, x, lambda = 0, allow_positive = FALSE) {
  ytx <- crossprod(y, x)
  ## In all cases, return 0 if abs(ytx) is too small
  if (abs(ytx) < lambda) return(0)
  ## In all cases, return (ytx + lambda) / crossprod(x) if ytx < -lambda
  if (ytx < -lambda) return( drop((ytx + lambda) / crossprod(x)) )
  ## Remaining cases: ytx > lambda
  ## If any x[i] > 0, allow_positive is void, return 0
  if (!allow_positive || any(x > 0)) return(0)
  ## Current case: ytx > lambda, allow_positive and all x[i] < 0
  ##               mitigate (ytx - lambda) / crossprod(x) by beta_max
  beta_max <- min(pmin(y, 0) / x, na.rm = TRUE) ## min_i (y[i]_- / x[i]_-)
  min(beta_max, drop((ytx - lambda) / crossprod(x)) )
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
#' @param prob a vector of probability weights for obtaining the coordinates
#' to be sampled
#' @param ... additional parameters
#'
#' @return the estimated value of beta
#' @export
solve_multivariate <- function(beta0, y, X, lambda, prob = NULL, ...) {
  p <- length(beta0)
  beta <- beta0

  yhat <- X %*% beta0
  fn_obj <- compute_objective_function(y, X, lambda)

  update_univariate2 <- function(beta, coord){
    betai <- beta[coord]
    beta_i <- beta[-coord]
    xi <- X[, coord]
    yi <- y - yhat + betai * xi

    beta_next <- beta
    beta_next[coord] <- solve_univariate(yi, xi, lambda)

    # mise à jour de Yhat
    diff_beta <- beta_next - beta
    yhat <<- yhat + diff_beta[coord] * xi

    beta_next
  }


  i <- 0
  eps <- 10 ^ -8
  progress <- +Inf
  obj <- fn_obj(beta0)

  while (i < p || progress > eps) {
    for (j in 1:p) {
      coord <- sample(p, size = 1, prob = prob)
      beta <- update_univariate2(beta, coord)
      i <- i + 1
    }
    new_obj <- fn_obj(beta)
    progress <- abs(new_obj - obj) / obj
    obj <- new_obj
  }


  value <- fn_obj(beta)

  list(par = beta, value = value, method = "shooting", iterations = i,
       last_progress = progress)
}
