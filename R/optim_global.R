#' Estimate shifts
#'
#' @param Delta0 initial position (size n+m)
#' @param zscores z-scores (size m)
#' @param incidence_mat incidence matrix (size m*(n+m))
#' @param covar_mat covariance matrix (size m*m)
#' @param lambda a positive regularization parameter
#'
#' @return a list
#' @export
#' @importFrom stats optim
estimate_shifts <- function(Delta0, zscores, incidence_mat, covar_mat, lambda){
  R <- inverse_sqrt(covar_mat)
  Y <- R %*% zscores
  X <- R %*% incidence_mat
  optim(par = Delta0,
        fn = compute_objective_function(Y, X, lambda),
        gr = compute_gradient_function(Y, X, lambda),
        upper = 0, method = "L-BFGS-B")
}
