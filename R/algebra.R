#' Inverse of the square root of a matrix
#'
#' Compute \eqn{R}, the inverse of the square root of a matrix \eqn{M}
#' such that \deqn{M^{-1} = R'R.}
#'
#' @param M a symmetric positive definite matrix
#'
#' @return a lower triangular matrix
#' @export
#'
#' @examples
#' M <- matrix(2, ncol = 3, nrow = 3) + diag(5, ncol = 3, nrow = 3)
#' R <- inverse_sqrt(M)
#' t(R) %*% R %*% M ## Should be the Identity
inverse_sqrt <- function(M){
  L <- t(chol(M))
  solve(L)
}

#' Norm 0
#'
#' Compute the number of non-zero values
#'
#' @param x Numeric.
#'
#' @return The number of non-zero values in x (\eqn{\|x\|_0}).
#' @export
#'
#' @examples
#' x <- c(1, 2, 3, 0, 4)
#' norm0(x)
norm0 <- function(x){
  stopifnot(is.numeric(x))
  sum(x != 0)
}

