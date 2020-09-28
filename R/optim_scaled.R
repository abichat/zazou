#' Scaled lasso
#'
#' Scaled lasso to compute \code{beta_init} and \code{hsigma}.
#'
#' @param y A vector of size m.
#' @param X A vector of size m*(n+m).
#' @param lambda A regularization parameter. If \code{NULL}, a universal
#' value is chosen.
#' @param beta0 Initial value for beta.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list composed by \code{par} (beta estimate) and
#' \code{sigma_scaledlasso} among others.
#' @export
#'
solve_scaled_lasso <- function(y, X, lambda = NULL, beta0, ...){

  n <- length(y)
  p <- ncol(X)
  beta <- beta0
  sigma <- var(y) / n

  if(is.null(lambda)){
    lambda <- sqrt(2 * log(p)/n) # As in Sun and Zhang (2012)
  }

  it <- 0
  eps <- 10 ^ -9
  progress <- +Inf
  obj <- compute_objective_function(Y = y, X = X,
                                    lambda = lambda, sigma = sigma,
                                    type = "scaled lasso")(beta)

  while (it < p || progress > eps) {
    beta <- solve_multivariate(y = y, X = X, beta0 = beta,
                               lambda = n * sigma * lambda, ...)$par
    sigma <- update_sigma(y, X, beta)
    new_obj <- compute_objective_function(Y = y, X = X,
                                          lambda = lambda, sigma = sigma,
                                          type = "scaled lasso")(beta)
    progress <- abs(new_obj - obj) / obj
    obj <- new_obj
    # cat(paste(c("Iteration:", it, round(sigma, 5), round(new_obj, 5))), "\n")
    it <- it + 1
  }

  list(par = beta, sigma_scaledlasso = sigma,
       method = "scaled lasso", value = obj,
       iterations = it, last_progress = progress)
}


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
