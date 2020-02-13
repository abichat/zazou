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
#' @param reverse If \code{TRUE}, returns the number of zeros in \code{x}.
#'
#' @return By default, the number of non-zero values in \code{x}
#' (\eqn{\|x\|_0}). If \code{reverse} is set to \code{TRUE},
#' the number of zero in \code{x}.
#' @export
#'
#' @examples
#' x <- c(1, 2, 3, 0, 4)
#' norm0(x)
norm0 <- function(x, reverse = FALSE) {
  stopifnot(is.numeric(x))
  n <- sum(x != 0)

  if (isTRUE(reverse)) {
    return(length(x) - n)
  } else {
    return(n)
  }
}

