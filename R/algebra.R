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
#' Compute the number of non-zero values.
#'
#' @param x Numeric.
#' @param rev If \code{TRUE}, returns the number of zeros in \code{x}.
#' @param prop If \code{TRUE}, returns the proportion of
#' (non-)zeros in \code{x}.
#'
#' @return By default, the number of non-zero values in \code{x}
#' (\eqn{\|x\|_0}).
#' If \code{rev} is set to \code{TRUE}, the number of zero in \code{x}.
#' If \code{prop} is set to \code{TRUE}, i return the proportion
#' instead of the absolute count.
#' @export
#'
#' @examples
#' x <- c(1, 2, 3, 0, 4)
#' norm0(x)
#' norm0(x, rev = TRUE)
#' norm0(x, prop = TRUE)
norm0 <- function(x, rev = FALSE, prop = FALSE) {
  stopifnot(is.numeric(x))
  stopifnot(length(x) >= 1)

  n <- sum(x != 0)
  l <- length(x)

  if (isTRUE(rev)) n <- l - n
  if (isTRUE(prop)) n <- n / l

  return(n)
}

