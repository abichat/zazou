#' Simulation of z-scores on branches
#'
#' Simulates z-scores
#'
#' Shifts are sampled from "statistically significant values of z-scores"
#'
#' @param tree tree.
#' @param alpha the parameter for the OU process.
#' @param shifts vector of shifts.
#' @param Nshifts number of shifts. Must be specified if \code{shifts} is not.
#'
#' @return a vector of simulated z-scores.
#' @export
#'
#' @importFrom stats rnorm runif
#' @importFrom ape Ntip
#' @examples
#' tree <- ape::rcoal(7)
#' shift = c(0, 0, -3, 0, 0, 0, 0, 0, 0, -2, 0, 0)
#' simu_zscores(tree, alpha = 0.2, shift = c(0, 0, -3, 0, 0, 0, 0, 0, 0, -2, 0, 0))
#'
simu_zscores <- function(tree, alpha, shifts = NULL, Nshifts = NULL){

  if(is.null(Nshifts)  && is.null(shifts)){
    stop("You must specified at least 'shifts' or 'Nshifts'.")
  }

  Nbranches <- nrow(tree$edge)

  if(is.null(shifts)){
    shifts <- rep(0, Nbranches)
    shifts[sample(length(shifts), Nshifts)] <- p2z(runif(Nshifts, max = 0.05))
  }

  if(nrow(tree$edge) != length(shifts)){
    stop("'shift' length must be length of branches in 'tree'.")
  }

  N <- Ntip(tree)

  incidence_mat <- incidence_matrix(tree)

  covar_mat <- covariance_matrix(tree, alpha = alpha)
  sqrtcovar <- inverse_sqrt(covar_mat)

  incidence_mat %*% shifts + sqrtcovar %*% rnorm(N)
}
