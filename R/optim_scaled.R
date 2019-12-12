#' Scaled lasso
#'
#' Scaled lasso to compute \code{beta_init} and \code{hsigma}.
#'
#' @param Y A vector of size m.
#' @param X A vector of size m*(n+m).
#' @param projected Logical. Should \code{beta_init} be projected on \eqn{\mathbb{R}_{-}} ?
#' @param ... Not used, here for avoiding errors when passing unknown arguments
#' from previous ellipsis.
#'
#' @return A list composed by \code{beta_init} and \code{hsigma}.
#' @export
#' @importFrom scalreg scalreg
#'
scaled_lasso <- function(Y, X, projected = TRUE, ...){
  object <- scalreg(X, Y)
  beta_init <- object$coefficients
  if(projected){
    beta_init[beta_init > 0] <- 0
  }
  list(beta_init = beta_init, hsigma = object$hsigma)
}
