#' Covariance matrix
#'
#' Compute normalized covariance matrix
#'
#' The coefficient is proportional to
#' \eqn{exp(- \alpha * d_{i,j}) - exp(-2 * \alpha * h) / (2 * \alpha)}
#'
#' @param tree a phylo object
#' @param alpha hyperparameter
#'
#' @return the covariance matrix of the OU process on the tree
#' @export
#'
#' @importFrom stats cophenetic
#'
#' @examples
#' tree <- ape::rtree(3)
#' covariance_matrix(tree, alpha = 1)
covariance_matrix <- function(tree, alpha){
  distance_matrix <- cophenetic(tree)
  h <- max(distance_matrix) / 2
  covar <- exp(- alpha * distance_matrix) - exp(-2 * alpha * h) / (2 * alpha)
  covar / covar[1, 1]
}
