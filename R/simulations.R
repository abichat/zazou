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
#' @param graph show the tree, shifts and z-scores be plotted?
#' @param ... additional parameters to pass to \code{plot_shifts()}.
#'
#' @return a vector of simulated z-scores.
#' @export
#'
#' @importFrom stats rnorm runif
#' @importFrom graphics plot
#' @importFrom ape Ntip
#'
#' @examples
#' set.seed(42)
#' tree <- ape::rcoal(5)
#' shift <- c(0, 0, -3, 0, 0, 0, 0, 0)
#' simu_zscores(tree, alpha = 0.2, shift = shift, graph = TRUE)
#' simu_zscores(tree, alpha = 0.2, Nshifts = 2, graph = TRUE)
simu_zscores <- function(tree, alpha, shifts = NULL, Nshifts = NULL,
                         graph = FALSE, ...){

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

  true_zs <- incidence_mat %*% shifts

  obs_zs <- true_zs + sqrtcovar %*% rnorm(N)

  if(graph){
    plot(plot_shifts(tree, shifts, true_scores = true_zs,
                     obs_scores = obs_zs, ...))
  }

  return(obs_zs)
}
