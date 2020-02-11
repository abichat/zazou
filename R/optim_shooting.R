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
#' @param u a vector of size n.
#' @param t a vector of size n. Column of the incidence matrix.
#' @param constraint_type Character. "beta" (default), "yhat" or "none".
#' Ensures that all coordinates of \eqn{\beta}
#' (for \code{constraint_type = "beta"}) or \eqn{z + x\beta}
#' (for \code{constraint_type = "yhat"}) are negative.
#' @inheritParams estimate_shifts
#'
#' @return The scalar solution \eqn{\beta} of the 1D optimization problem
#' @export
#'
#' @examples
#' solve_univariate(1:4, -(4:1), 2)
solve_univariate <- function(y, x, z = rep(0, length(y)), u, t, lambda = 0,
                             constraint_type = c("beta", "yhat", "none"),
                             ...) {
  constraint_type <- match.arg(constraint_type)
  ## Compute unconstrained estimator
  ytx <- sum((y - z) * x) ## crossprod(y -z, x)
  if (abs(ytx) <= lambda) beta <- 0
  if (ytx < -lambda) beta <- (ytx + lambda) / sum(x^2)
  if (ytx > lambda)  beta <- (ytx - lambda) / sum(x^2)
  ## Check constraint
  if (constraint_type == "none") return(beta)
  ## Compute feasible set
  if (constraint_type == "beta") {
    beta_min <- -Inf
    beta_max <- 0
  } else { ## constraint_type == yhat
    ## Rounding to avoid numerical errors
    u[abs(u) < 10 * .Machine$double.eps] <- 0
    t[abs(t) < 10 * .Machine$double.eps] <- 0
    ## Upper bound of feasible set: min_{i: t[i]>0} (-u[i] / t[i])
    t_plus <- t > 0
    if (any(t_plus)) {
      beta_max <- min(-u[t_plus] / t[t_plus])
    } else {
      beta_max <- Inf
    }
    ## Lower bound of feasible set: max_{i: t[i]<0} (-u[i] / t[i])
    t_minus <- t < 0
    if (any(t_minus)) {
      beta_min <- max(-u[t_minus] / t[t_minus])
    } else {
      beta_min <- -Inf
    }
  }
  ## Check that feasible set (u + t * beta <= 0) is not empty
  if (beta_min - beta_max > sqrt(.Machine$double.eps)) {
    stop("The constraint is not feasible. Consider changing the constraint.")
  }
  max(beta_min, min  ## Project unconstrained estimate beta to feasible set [beta_min, beta_min]
(beta, beta_max))
}

#' @rdname solve_univariate
#'
#' @param beta0 the initial position of beta.
#' @param X a matrix of size  x p.
#' @param max_it Maximum number of iterations.
#' @param prob a vector of probability weights for obtaining the coordinates
#' to be sampled.
#' @param ... further arguments passed to or from other methods.
#'
#' @return the estimated value of beta
#' @export
solve_multivariate <- function(beta0, y, X, lambda, prob = NULL,
                               max_it = 500, incidence_mat,
                               constraint_type = c("beta", "yhat", "none"),
                               ...) {
  constraint_type <- match.arg(constraint_type)
  p <- length(beta0)
  beta <- beta0

  yhat <- X %*% beta
  if(constraint_type == "yhat") {
    I <- incidence_mat
    J <- I %*% beta
  }
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

    if(constraint_type == "yhat") {
      Ji <- J[coord]
      Ii <- I[, coord]
      ui <- J - betai * Ii
    }

    # update betai
    if(constraint_type != "yhat") {
      betai <- solve_univariate(y = y, x = xi, z = zi, lambda = lambda, ...)
    } else {
      betai <- solve_univariate(y = y, x = xi, z = zi,
                                u = ui, t = Ii, lambda = lambda, ...)
      J <<- ui + betai * Ii
    }

    # update beta
    beta[coord] <<- betai
    # update yhat
    yhat <<- zi + betai * xi ## X %*% beta
  }

  it <- 1
  eps <- 10 ^ -8
  progress <- +Inf

  ## Keep track of objective value
  obj_vals <- fn_obj(beta0)

  while (it < max_it && progress > eps) {
    ## Update coordinates in random order (rather than random coordinates)
    coord_order <- sample(p)
    for (coord in coord_order) {
      update_coord(coord, ...)
    }
    ## Store current objective value and compute progress
    it <- it + 1
    new_obj_vals <- fn_obj(beta)
    progress <- abs(new_obj_vals - obj_vals) / obj_vals
    obj_vals <- new_obj_vals
  }
  beta <- data.frame(estimate = beta)
  list(par = beta, value = obj_vals, method = "shooting",
       iterations = it, last_progress = progress)
}
