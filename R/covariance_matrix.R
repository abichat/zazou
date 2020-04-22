#' Covariance matrix
#'
#' Compute normalized covariance matrix
#'
#' The coefficient is proportional to
#' \eqn{exp(- \alpha * d_{i,j}) - exp(-2 * \alpha * h) / (2 * \alpha)}
#'
#' @param tree a phylo object
#' @param alphaOU hyperparameter
#'
#' @return the covariance matrix of the OU process on the tree
#' @export
#'
#' @importFrom stats cophenetic
#'
#' @examples
#' tree <- ape::rtree(3)
#' covarianceOU_matrix(tree, alpha = 1)
covarianceOU_matrix <- function(tree, alphaOU){
  distance_matrix <- cophenetic(tree)
  h <- max(distance_matrix) / 2
  covar <-
    (exp(- alphaOU * distance_matrix) - exp(-2 * alphaOU * h)) / (2 * alphaOU)
  covar / covar[1, 1]
}
