#' Scaled lasso
#'
#' Scaled lasso to compute \code{beta_init} and \code{hsigma}.
#'
#' @param y A vector of size m.
#' @param X A vector of size m*(n+m).
#' @param projected Logical. Should each coordinate of \code{beta_init} be projected on \eqn{\mathbb{R}_{-}}?
#' @param ... Not used, here for avoiding errors when passing unknown arguments
#' from previous ellipsis.
#'
#' @return A list composed by \code{beta_init} and \code{hsigma}.
#' @export
#' @importFrom scalreg scalreg
#'
scaled_lasso <- function(y, X, projected = TRUE, ...){
  object <- scalreg(X, y)
  beta_init <- object$coefficients
  if(projected){
    beta_init[beta_init > 0] <- 0
  }
  list(beta_init = beta_init, hsigma = object$hsigma)
}

update_sigma <- function(y, X, beta){
  sqrt(sum((y - X %*% beta)^2) / nrow(X))
}

scaled_lasso2 <- function(beta0, y, X, lambda, ...){

  n <- length(y)
  p <- ncol(X)
  beta <- beta0
  sigma <- var(y) / n

  i <- 0
  eps <- 10 ^ -9
  progress <- +Inf
  obj <- compute_objective_function(y, X, n * sigma * lambda)(beta)

  while (i < p || progress > eps) {
    beta <- solve_multivariate(beta, y, X, n * sigma * lambda, ...)$par
    sigma <- update_sigma(y, X, beta)
    new_obj <- compute_objective_function(y, X, n * sigma * lambda)(beta)
    progress <- abs(new_obj - obj) / obj
    obj <- new_obj
    # cat(paste(c("Iteration:", i, round(sigma, 5), round(new_obj, 5))), "\n")
    i <- i + 1
  }
  list(beta = beta, sigma = sigma, obj = obj)
}
