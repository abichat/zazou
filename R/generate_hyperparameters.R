#' Generate a logarithmic lambda grid
#'
#' @param x Design matrix.
#' @param y Response vector.
#' @param n_lambda Length of the lambda grid.
#' @param min_ratio Ratio between the minimum and maximum values in
#' the lambda grid.
#' @param ... Not used, here for avoiding errors when passing unknown arguments
#' from previous ellipsis.
#'
#' @export
#' @return a numeric vector of length \code{nlambda}
#'
lambda_grid <- function(x, y, n_lambda = 10, min_ratio = 1e-2, ...) {
  ## center y and x
  y <- scale(y, center = TRUE, scale = FALSE)
  x <- scale(x, center = TRUE, scale = FALSE)
  xty <- drop(crossprod(x, y))
  lambda_max <- max(abs(xty))
  return(10 ^ seq(from = log10(lambda_max), to = log10(min_ratio * lambda_max),
                  len = n_lambda))
}


#' Gamma
#'
#' Generate a gamma for the column-wise inverse procedure.
#'
#' @param X Design matrix, size \code{m*(n+m)}.
#' @param factor Numeric.
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return A suggested value of \eqn{\gamma}.
#' @export
#'
generate_gamma <- function(X, factor = 2, ...){
  stopifnot(factor > 0)
  factor * sqrt(log(ncol(X)) / nrow(X))
}
