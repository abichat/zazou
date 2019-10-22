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
#' @param prob a vector of probability weights for obtaining the coordinates
#' to be sampled
#'
#' @return the estimate value of beta
#' @export
solve_multivariate <- function(beta0, y, X, lambda, prob = NULL) {
  p <- length(beta0)
  beta <- beta0
  for(i in 1:200){
    coord <- sample(p, size = 1, prob = prob)
    beta <- update_univariate(beta, coord, y, X, lambda)
  }

  fn_obj <- compute_objective_function(y, X, lambda)
  value <- as.vector(fn_obj(beta))

  list(par = beta, value = value, method = "shooting", iterations = i)
}

#' @rdname solve_univariate
#'
#' @return the estimate value of beta
#' @export
solve_multivariate2 <- function(beta0, y, X, lambda, prob = NULL) {
  p <- length(beta0)
  beta <- beta0

  yhat <- X %*% beta0

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

  for(i in 1:200){
    coord <- sample(p, size = 1, prob = prob)
    beta <- update_univariate2(beta, coord)
  }

  fn_obj <- compute_objective_function(y, X, lambda)
  value <- fn_obj(beta)

  list(par = beta, value = value, method = "shooting", iterations = i)
}

#' @rdname solve_univariate
#'
#' @return the estimate value of beta
#' @export
solve_multivariate3 <- function(beta0, y, X, lambda, prob = NULL) {
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
  never_null <- TRUE

  while (i < p || progress > eps) {
    for (j in 1:(ceiling(p / 10) + 20)) {
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
