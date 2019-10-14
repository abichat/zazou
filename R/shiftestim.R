
#' 'shiftestim' object
#'
#' @param listopt an output of 'estimate_shifts()'
#' @param tree tree
#' @param zscores vector of observed z-scores
#' @param lambda regularization parameter
#' @param alpha parameter
#' @param covar_mat the covariance matrix
#'
#' @return a 'shiftestim' object
#' @export
as_shiftestim <- function(listopt, tree, zscores, lambda, alpha, covar_mat) {

  required_names <- c("par", "value", "method")

  if (!is.list(listopt) ||
      !all(required_names %in% names(listopt))) {
        stop("'listopt' must be the output of a optimisation function.")
  }

  zscores_est <- incidence_matrix(tree) %*% listopt$par

  obj <- list(zscores_obs = zscores,
              zscores_est = zscores_est,
              shift_est = listopt$par,
              objective_value = listopt$value, lambda = lambda,
              tree = tree, alpha = alpha, covar_mat = covar_mat,
              method = listopt$method,
              optim_info = listopt[setdiff(names(listopt), required_names)])
  class(obj) <- "shiftestim"
  return(obj)
}




#' @rdname as_shiftestim
#'
#' @name print.shiftestim
#'
#' @param x a 'shiftestim' object.
#' @param ... further arguments to be passed to or from other methods.
#' @inheritParams plot_shifts
#'
#' @importFrom utils head
#'
#' @export
print.shiftestim <- function(x, digits = 3, ...){

  if(is.null(x$alpha)){
    txt_alpha <- "Covariance matrix has been manually specified."
  } else {
    h <- tree_height(x$tree)
    sigma <- sqrt(2 * x$alpha) / (1 - exp(- 2 * x$alpha * h))
    txt_alpha <- paste0("Covariance matrix has been estimated from an OU",
                        " with alpha = ", round(x$alpha, digits),
                        " and sigma = ", round(sigma, digits), ".")
  }



  cat(txt_alpha, "\n")
  cat("Optimisation algorithm: ", x$method, ".\n", sep = "")
  cat("---\n")
  cat("Estimated shifts:", head(round(x$shift_est, digits), 10), "...\n")
  cat(sum(x$shift_est != 0), "shifts have been identified (ie",
      100 * round(mean(x$shift_est == 0), digits), "% of sparsity).\n")
  cat("---\n")
  cat("Estimated z-scores:", head(round(x$zscores_est, digits), 10), "...\n")
  cat("Observed z-scores: ", head(round(x$zscores_obs, digits), 10), "...\n")
}


#' @rdname as_shiftestim
#' @inheritParams plot_shifts
#'
#' @export
plot.shiftestim <- function(x, digits = 3, ...){
  plot_shifts(tree = x$tree, shifts = x$shift_est,
              obs_scores = x$zscores_obs, est_scores = x$zscores_est,
              digits = digits, ...)
}


