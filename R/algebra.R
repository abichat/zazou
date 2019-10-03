#' Inverse of the square root of a matrix
#'
#' Compute \eqn{R}, the inverse of the square root of a matrix \eqn{M} such that \deqn{M^{-1} = R'R.}
#'
#' @param M a symmetric positive definite matrix
#'
#' @return a lower triangular matrix
#' @export
#'
#' @examples
#' M <- matrix(2, ncol = 3, nrow = 3) + diag(5, ncol = 3, nrow = 3)
#' R <- inverse_sqrt(M)
#' t(R) %*% R %*% M
inverse_sqrt <- function(M){
  L <- t(chol(M))
  solve(L)
}

