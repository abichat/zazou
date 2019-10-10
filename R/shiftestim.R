
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
  if (!is.list(listopt) ||
      !all(names(listopt) == c("par", "value", "counts",
                               "convergence", "message"))) {
        stop("'listopt' must be the output of 'estimate_shifts()'")
  }
  obj <- list(zscores_obs = zscores, shift_est = listopt$par,
              objective_value = listopt$value, lambda = lambda,
              tree = tree, alpha = alpha, covar_mat = covar_mat,
              optim_info = listopt[c("counts", "convergence", "message")])
  class(obj) <- "shiftestim"
  return(obj)
}




#' @rdname as_shiftestim
#'
#' @name print.shiftestim
#'
#' @param x a 'shiftestim' object.
#' @param digits number of digits to round to.
#' @param ... further arguments to be passed to or from other methods.
#'
#' @importFrom utils head
#'
#' @export
print.shiftestim <- function(x, digits = 4, ...){

  if(is.null(x$alpha)){
    txt_alpha <- "Covariance matrix has been manually specified."
  } else {
    h <- tree_height(x$tree)
    sigma <- sqrt(2 * x$alpha) / (1 - exp(- 2 * x$alpha * h))
    txt_alpha <- paste0("Covariance matrix has been estimated from an OU",
                        " with alpha = ", round(x$alpha, digits),
                        " and sigma = ", round(sigma, digits), ".")
  }



  zscores_est <- incidence_matrix(x$tree) %*% x$shift_est

  cat(txt_alpha, "\n")

  cat("Estimated shifts:", head(round(x$shift_est, digits), 10), "...\n")
  cat(sum(x$shift_est != 0), "shifts have been identified (ie",
      100 * round(mean(x$shift_est == 0), digits), "% of sparsity).\n")
  cat("Estimated z-scores:", head(round(zscores_est, digits), 10), "...\n")
  cat("Observed z-scores: ", head(round(x$zscores_obs, digits), 10), "...\n")
}


#' @rdname as_shiftestim
#'
#' @export
plot.shiftestim <- function(x, digits = 4, ...){
  plot_shifts(shifts = round(x$shift_est, digits),
              tree = x$tree, scores = x$zscores_obs)
}


