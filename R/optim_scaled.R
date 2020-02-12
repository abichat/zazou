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
#' @return A list composed by \code{par} (beta estimate) and
#' \code{sigma_scaledlasso} among others.
#' @export
#'
scaled_lasso <- function(beta0, y, X, lambda = NULL, ...){

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
                                    type = "scaledlasso")(beta)

  while (it < p || progress > eps) {
    beta <- solve_multivariate(beta0 = beta, y = y, X = X,
                               lambda = n * sigma * lambda, ...)$par$estimate
    sigma <- update_sigma(y, X, beta)
    new_obj <- compute_objective_function(Y = y, X = X,
                                          lambda = lambda, sigma = sigma,
                                          type = "scaledlasso")(beta)
    progress <- abs(new_obj - obj) / obj
    obj <- new_obj
    # cat(paste(c("Iteration:", it, round(sigma, 5), round(new_obj, 5))), "\n")
    it <- it + 1
  }

  list(par = data.frame(estimate = beta), sigma_scaledlasso = sigma,
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
