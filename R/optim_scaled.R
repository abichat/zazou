#' Update sigma in the scaled lasso algorithm
#'
#' Compute the exact value of sigma that nullify the gradient.
#'
#' @param y A vector of size m.
#' @param X A vector of size m*(n+m).
#' @param beta The current value of beta.
#'
#' @return Updated version of sigma.
update_sigma <- function(y, X, beta){
  sqrt(sum((y - X %*% beta)^2) / nrow(X))
}

#' Scaled lasso
#'
#' Scaled lasso to compute \code{beta_init} and \code{hsigma}.
#'
#' @param beta0 Initial value for beta.
#' @param y A vector of size m.
#' @param X A vector of size m*(n+m).
#' @param lambda A regularization parameter. If \code{NULL}, a universal
#' value is chosen.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list composed by \code{beta_init} and \code{hsigma}.
#'
scaled_lasso <- function(beta0, y, X, lambda = NULL, ...){

  n <- length(y)
  p <- ncol(X)
  beta <- beta0
  sigma <- var(y) / n

  if(is.null(lambda)){
    lambda <- sqrt(2 * log(p)/n) # As in Sun and Zhang (2012)
  }

  i <- 0
  eps <- 10 ^ -9
  progress <- +Inf
  obj <- compute_objective_function(y, X, n * sigma * lambda,
                                    type = "lasso")(beta)

  while (i < p || progress > eps) {
    beta <- solve_multivariate(beta, y, X, n * sigma * lambda, ...)$par$estimate
    sigma <- update_sigma(y, X, beta)
    new_obj <- compute_objective_function(y, X, n * sigma * lambda,
                                          type = "lasso")(beta)
    progress <- abs(new_obj - obj) / obj
    obj <- new_obj
    # cat(paste(c("Iteration:", i, round(sigma, 5), round(new_obj, 5))), "\n")
    i <- i + 1
  }
  list(beta_init = beta, hsigma = sigma)
}
