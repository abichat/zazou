#' Estimate shifts
#'
#' @rdname estimate_shifts
#'
#' @param zscores z-scores (size m)
#' @param lambda a grid of positive regularization parameters
#' @param tree the tree used to compute the incidence and covariance matrices
#' (if not specified)
#' @param alphaOU a grid of positive alpha parameters
#' @param method method to use for the optimization. One of \code{lasso},
#' \code{scaled lasso}, or \code{L-BFGS-B}.
#' @param criterion criterion on which the selection is done.
#' @param beta0 initial position (size n+m)
#' @param constraint_type Constrains on shifts.
#' @param ... further arguments passed to or from other methods.
#'
#' @return an object of class \code{shiftpunct}
#'
#'
#' @export
#' @importFrom stats optim
estimate_shifts <- function(zscores, tree, alphaOU, lambda = NULL,
                            method = c("lasso", "scaled lasso", "L-BFGS-B"),
                            criterion = c("bic", "pbic"),
                            beta0 = rep(0, length(tree$edge.length)),
                            constraint_type = c("beta", "yhat", "none"),
                            ...){

  method <- match.arg(method)
  criterion <- match.arg(criterion)
  constraint_type <- match.arg(constraint_type)

  ## local estimation routine (for a single lambda)
  fitting_procedure <- function(beta0, X, Y, lambda, ...) {

    opt <- switch(method,
                  "lasso" = solve_multivariate(
                    y = Y, X = X, lambda = lambda, beta0 = beta0,
                    constraint_type = constraint_type, ...),

                  "scaled lasso" = solve_scaled_lasso(
                    y = Y, X = X, lambda = lambda, beta0 = beta0,
                    constraint_type = constraint_type, ...),

                  "L-BFGS-B" = solve_lbfgsb(
                    X = X, Y = Y, lambda = lambda, beta0 = beta0,
                    constraint_type = constraint_type, ...)

                   )
    return(opt)
  }

  ## Bookkeeping variables
  best_criterion <- Inf
  bic_df <- data.frame(alphaOU  = numeric(0),
                       lambda = numeric(0),
                       objective_value = numeric(0),
                       bic    = numeric(0),
                       pbic   = numeric(0))
  best_model <- NULL
  shifts <- list()

  ## Outer loop on alphaOU
  for (alp in alphaOU) {
    ## Compute covariance matrices, design matrix and response vector
    mat_covarOU <- covarianceOU_matrix(tree, alp)
    mat_incidence <- incidence_matrix(tree)
    R <- inverse_sqrt(mat_covarOU)
    Y <- R %*% zscores
    X <- R %*% mat_incidence

    ## Set lambda grid for inner loop on lambda
    if (is.null(lambda)) {
      current_lambda <- lambda_grid(x = X, y = Y, ...)
    } else {
      current_lambda <- lambda
    }

    ## Inner loop on lambda
    for (lam in current_lambda) {
      ## Compute current model
      opt <- fitting_procedure(beta0 = beta0, X = X, Y = Y, lambda = lam,
                               mat_incidence = mat_incidence, ...)
      current_model <- as_shiftpunct(
        listopt = opt, tree = tree, zscores = zscores,
        lambda = lam, alphaOU = alp)
      shifts <- c(shifts, list(current_model$shifts_est))

      ## Update bic table
      bic_df <- rbind(bic_df,
                      data.frame(alphaOU = alp, lambda = lam,
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
  } ## Close alphaOU loop

  if (any(c(length(alphaOU), length(current_lambda)) >= 2)) {
    bic_df$shifts_est <- shifts
    best_model$optim_info$bic_selection <- bic_df
    best_model$optim_info$criterion <- criterion
    best_model$method <- paste(method, "with model selection")
  }

  best_model$optim_info$supp_arg <- list(...)
  best_model$optim_info$supp_arg$alpha_conf <- NULL

  return(best_model)
}
