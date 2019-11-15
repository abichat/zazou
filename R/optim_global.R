#' Estimate shifts
#'
#' @rdname estimate_shifts
#'
#' @param Delta0 initial position (size n+m)
#' @param zscores z-scores (size m)
#' @param lambda a grid of positive regularization parameters
#' @param tree the tree used to compute the incidence and covariance matrices
#' (if not specified)
#' @param alpha a grid of positive alpha parameters
#' @param method method to use for the optimization. One of \code{L-BFGS-B}
#' or \code{shooting}.
#' @param ... additional parameters
#'
#' @return an object of class \code{shiftestim}
#'
#' @export
estimate_shifts <- function(Delta0, zscores, tree, alpha, lambda = NULL,
                             method = c("L-BFGS-B", "shooting"), ...){

  method <- match.arg(method)

  ## local estimation routine (for a single lambda)
  fitting_procedure <- function(Delta0, X, Y, lambda, ...) {
    if (method == "L-BFGS-B") {
      opt <- optim(par = Delta0,
                   fn = compute_objective_function(Y, X, lambda),
                   gr = compute_gradient_function(Y, X, lambda),
                   upper = 0, method = "L-BFGS-B")
      opt <- c(opt, method = "L-BFGS-B")
    }
    if (method == "shooting") {
      opt <- solve_multivariate(Delta0, Y, X, lambda)
    }
    return(opt)
  }

  ## Bookkeeping variables
  best_bic <- Inf
  bic_df <- data.frame(alpha  = numeric(0),
                       lambda = numeric(0),
                       bic    = numeric(0))
  best_model <- NULL
  shifts <- list()

  ## Outer loop on alpha
  for (alp in alpha) {
    ## Compute covariance matrices, design matrix and response vector
    covar_mat <- covariance_matrix(tree, alp)
    incidence_mat <- incidence_matrix(tree)
    R <- inverse_sqrt(covar_mat)
    Y <- R %*% zscores
    X <- R %*% incidence_mat

    ## Set lambda grid for inner loop on lambda
    if (is.null(lambda)) {
      current_lambda <- lambda_grid(X, Y)
    } else {
      current_lambda <- lambda
    }

    ## Inner loop on lambda
    for (lam in current_lambda) {
      ## Compute current model
      opt <- fitting_procedure(Delta0, X, Y, lambda = lam, ...)
      current_model <- as_shiftestim(
        listopt = opt, tree = tree, zscores = zscores,
        lambda = lam, alpha = alp, covar_mat = covar_mat
      )
      shifts <- c(shifts, list(current_model$shift_est))
      ## Update bic table
      bic_df <- rbind(bic_df, data.frame(alpha = alp,
                                         lambda = lam,
                                         bic = current_model$bic))
      ## Update best model
      if (current_model$bic < best_bic) {
        best_model <- current_model
        best_bic   <- current_model$bic
      }
    } ## Close lambda loop
  } ## Close alpha loop

  if(any(c(length(alpha), length(current_lambda)) >= 2)){
    bic_df$shift_est <- shifts
    best_model$optim_info$bic_selection <- bic_df
    best_model$method <- paste(method, "with model selection")
  }

  return(best_model)
}
