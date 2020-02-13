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
#' @param criterion criterion on which the selection is done.
#' @param constraint_type Constrains on shifts.
#' @param ... further arguments passed to or from other methods.
#'
#' @return an object of class \code{shiftestim}
#'
#'
#' @export
#' @importFrom stats optim
estimate_shifts <- function(Delta0, zscores, tree, alpha, lambda = NULL,
                            method = c("L-BFGS-B", "shooting",
                                       "scaledlasso", "desparsifiedlasso"),
                            criterion = c("bic", "pbic"),
                            constraint_type = c("beta", "yhat", "none"), ...){

  method <- match.arg(method)
  criterion <- match.arg(criterion)
  constraint_type <- match.arg(constraint_type)

  if(method %in% c("scaledlasso", "desparsifiedlasso")){
    if(any(c(length(alpha), length(lambda)) != 1)){
      stop(paste0("No model selection can be done with ", method, "."))
    }
  }

  ## local estimation routine (for a single lambda)
  fitting_procedure <- function(Delta0, X, Y, lambda, ...) {

    opt <- switch (method,
                   "L-BFGS-B" = solve_lbfgsb(Delta0 = Delta0, X = X, Y = Y,
                                             lambda = lambda,
                                             constraint_type = constraint_type, ...),
                   "shooting" = solve_multivariate(beta0 = Delta0, y = Y, X = X,
                                                   lambda = lambda,
                                                   constraint_type = constraint_type, ...),
                   "scaledlasso" = scaled_lasso(beta0 = Delta0, y = Y, X = X,
                                                lambda = lambda,
                                                constraint_type = constraint_type, ...),
                   "desparsifiedlasso" =
                     solve_desparsified(Delta0 = Delta0, Y = Y, X = X,
                                        lambda = lambda,
                                        constraint_type = constraint_type, ...))
    return(opt)
  }

  ## Bookkeeping variables
  best_criterion <- Inf
  bic_df <- data.frame(alpha  = numeric(0),
                       lambda = numeric(0),
                       objective_value = numeric(0),
                       bic    = numeric(0),
                       pbic   = numeric(0))
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
      current_lambda <- lambda_grid(x = X, y = Y, ...)
    } else {
      current_lambda <- lambda
    }

    ## Inner loop on lambda
    for (lam in current_lambda) {
      ## Compute current model
      opt <- fitting_procedure(Delta0 = Delta0, X = X, Y = Y, lambda = lam,
                               incidence_mat = incidence_mat, ...)
      current_model <- as_shiftestim(
        listopt = opt, tree = tree, zscores = zscores,
        lambda = lam, alpha = alp)
      shifts <- c(shifts, list(current_model$shift_est))

      ## Update bic table
      bic_df <- rbind(bic_df,
                      data.frame(alpha = alp, lambda = lam,
                                 objective_value =
                                   current_model$objective_value,
                                 bic = current_model$bic,
                                 pbic = current_model$pbic))

      ## Update best model
      current_criterion <- getElement(current_model, criterion)
      if (is.finite(current_criterion) &&
          current_criterion < best_criterion) {
        best_model <- current_model
        best_criterion <- current_criterion
      }
    } ## Close lambda loop
  } ## Close alpha loop

  if(any(c(length(alpha), length(current_lambda)) >= 2)){
    bic_df$shift_est <- shifts
    best_model$optim_info$bic_selection <- bic_df
    best_model$optim_info$criterion <- criterion
    best_model$method <- paste(method, "with model selection")
  }

  best_model$optim_info$supp_arg <- list(...)

  return(best_model)
}
