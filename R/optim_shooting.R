#' Solve unidirectional constrained problem
#'
#' This function minimizes \eqn{\beta} in the 1D problem
#' \eqn{1/2 * ||y - z - x \beta||_2^2 + \lambda |\beta|} subject to
#' either \eqn{\beta <= 0} or \eqn{z + x\beta <= 0} (coordinate wise).
#'
#' The analytical solution of this problem is given by
#' \deqn{\beta* = min(0, ((y-z)'x + \lambda) / x'x ).}
#' when using the first constraint and is slightly more
#' complex when using the second constraint
#' (refer to the corresponding vignette)
#'
#' @param y a vector of size n.
#' @param x a vector of size n.
#' @param z a vector of size n.
#' @param use_constraint Logical. Default TRUE. If \code{TRUE}, the
#' return value \eqn{\beta} satisfies either \eqn{\beta <= 0} or
#' \eqn{z + x\beta <= 0} coordinate wise.
#' @param constraint_type Either "beta" (default) or "yhat". Ensures that
#' all coordinates of \eqn{\beta} (for \code{constraint_type = "beta"}) or
#' \eqn{z + x\beta} (for \code{constraint_type = "yhat"}) are negative.
#' Not used if \code{use_constraint} is set to \code{FALSE}.
#' @inheritParams estimate_shifts
#'
#' @return The scalar solution \eqn{\beta} of the 1D optimization problem
#' @export
#'
#' @examples
#' solve_univariate(1:4, -(4:1), 2)
solve_univariate <- function(y, x, z = rep(0, length(y)), lambda = 0,
                             use_constraint = TRUE,
                             constraint_type = c("beta", "yhat"),
                             ...) {
  constraint_type <- match.arg(constraint_type)
  allow_positive <- switch(constraint_type,
    "beta" = FALSE,
    "yhat" = TRUE
  )
  ytx <- crossprod(y - z, x)
  ## In all cases, return 0 if abs(ytx) is too small
  if (abs(ytx) < lambda) {
    return(0)
  }
  ## In all cases, return (ytx + lambda) / crossprod(x) if ytx < -lambda
  if (ytx < -lambda) {
    return(drop((ytx + lambda) / crossprod(x)))
  }
  ## Remaining cases: ytx > lambda
  ## - Without negativity constraint return shrinked estimate
  if (!use_constraint) {
      return(drop((ytx - lambda) / crossprod(x)))
  }
  ## - If any x[i] & z[i] > 0, allow_positive is void, return 0
  if (!allow_positive || any( (x > 0) & (z > 0))) {
    return(0)
  }
  ## - Else, ytx > lambda and constraint on yhat
  ##               mitigate (ytx - lambda) / crossprod(x) by feasible set
  ## Upper bound of feasible set: min_{i: x[i]>0} (-z[i] / x[i])
  x_plus <- x > sqrt(.Machine$double.eps)
  if (any(x_plus)) {
    beta_max <- min(-z[x_plus] / x[x_plus])
  } else {
    beta_max <- Inf
  }
  ## Lower bound of feasible set: max_{i: x[i]<0} (-z[i] / x[i])
  x_minus <- x < -sqrt(.Machine$double.eps)
  if (any(x_minus)) {
    beta_min <- max(-z[x_minus] / x[x_minus])
  } else {
    beta_min <- -Inf
  }
  ## Check that z + x * beta <= 0 is feasible.
  if (beta_min > beta_max) {
    stop("The constraint is not feasible. Consider changing the constraint.")
  }
  min(beta_max, drop((ytx - lambda) / crossprod(x)))
}

#' @rdname solve_univariate
#'
#' @param beta0 the initial position of beta.
#' @param X a matrix of size  x p.
#' @param max.it Maximum number of iterations.
#' @param prob a vector of probability weights for obtaining the coordinates
#' to be sampled.
#' @param ... further arguments passed to or from other methods.
#'
#' @return the estimated value of beta
#' @export
solve_multivariate <- function(beta0, y, X, lambda, prob = NULL,
                               max.it = 100, ...) {
  p <- length(beta0)
  beta <- beta0

  yhat <- X %*% beta0
  # fn_obj <- compute_objective_function(y, X, lambda)
  ## Fast alternative to compute_objective_function leveraging
  ## the fact that we have access to yhat (= X %*% beta)
  fn_obj <- function(beta) {
    sum( (y - yhat)^2 ) / 2 + lambda * sum(abs(beta))
  }

  ## update_coord only has side effects
  update_coord <- function(coord, ...){
    betai <- beta[coord]
    xi <- X[, coord]
    zi <- yhat - betai * xi ## X[ , -coord] %*% beta[-coord]
    # update betai
    betai <- solve_univariate(y = y, x = xi, z = zi, lambda = lambda, ...)
    # update beta
    beta[coord] <<- betai
    # update yhat
    yhat <<- zi + betai * xi ## X %*% beta
  }

  it <- 1
  eps <- 10 ^ -8
  progress <- +Inf

  ## Keep track of objective values throughout the iterations
  obj_vals <- rep(NA, max.it)
  obj_vals[it] <- fn_obj(beta0)

  while (it < max.it && progress > eps) {
    ## Update coordinates in random order (rather than random coordinates)
    coord_order <- sample(p)
    for (coord in coord_order) {
      update_coord(coord, ...)
    }
    ## Store current objective value and compute progress
    it <- it + 1
    obj_vals[it] <- fn_obj(beta)
    progress <- abs(obj_vals[it] - obj_vals[it - 1]) / obj_vals[it - 1]
  }

  list(par = beta, value = obj_vals[it], method = "shooting", iterations = it,
       objective_values = obj_vals[1:it], last_progress = progress)
}
