#' Estimate shifts
#'
#' @param Delta0 initial position (size n+m)
#' @param zscores z-scores (size m)
#' @param incidence_mat incidence matrix (size m*(n+m))
#' @param covar_mat covariance matrix (size m*m) between z-scores
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


#' @rdname estimate_shifts
#' @inheritParams estimate_shifts
#' @param tree the tree form which is calculated the incidence matrix and the
#' covariance matrix if not specified
#' @param alpha parameter > 0
#' @param method method to use. One of \code{L-BFGS-B} or \code{shooting}.
#'
#' @return a shiftestim object
#' @export
#' @importFrom stats optim
estimate_shifts2 <- function(Delta0, zscores, tree, lambda = 0,
                             alpha = NULL, covar_mat = NULL,
                             method = c("L-BFGS-B", "shooting",
                                        "shooting2", "shooting3"), ...){

  method <- match.arg(method)

  if(is.null(alpha) && is.null(covar_mat)){
    stop("Either 'alpha' or 'covar_mat' must be specified.")
  }

  specified_covar <- TRUE
  if(is.null(covar_mat)){
    covar_mat <- covariance_matrix(tree, alpha)
    specified_covar <- FALSE
  }

  incidence_mat <- incidence_matrix(tree)
  R <- inverse_sqrt(covar_mat)
  Y <- R %*% zscores
  X <- R %*% incidence_mat

  if(method == "L-BFGS-B"){
    opt <- optim(par = Delta0,
                 fn = compute_objective_function(Y, X, lambda),
                 gr = compute_gradient_function(Y, X, lambda),
                 upper = 0, method = "L-BFGS-B")
    opt <- c(opt, method = "L-BFGS-B")
  }
  if(method == "shooting"){
    opt <- solve_multivariate(Delta0, Y, X, lambda)
  }
  if(method == "shooting2"){
    opt <- solve_multivariate2(Delta0, Y, X, lambda)
  }
  if(method == "shooting3"){
    opt <- solve_multivariate3(Delta0, Y, X, lambda, ...)
  }

  if(specified_covar) alpha <- NULL

  as_shiftestim(listopt = opt, tree = tree, zscores = zscores,
                lambda = lambda, alpha = alpha, covar_mat = covar_mat)
}

#' @rdname estimate_shifts
#' @param ... additional parameters
#' @export
estimate_shifts3 <- function(Delta0, zscores, tree, alpha, lambda = 0,
                             method = c("L-BFGS-B", "shooting",
                                        "shooting2", "shooting3"), ...){

  method <- match.arg(method)

  if(length(lambda) == 1 && length(alpha) == 1){

    covar_mat <- covariance_matrix(tree, alpha)
    incidence_mat <- incidence_matrix(tree)
    R <- inverse_sqrt(covar_mat)
    Y <- R %*% zscores
    X <- R %*% incidence_mat

    if(method == "L-BFGS-B"){
      opt <- optim(par = Delta0,
                   fn = compute_objective_function(Y, X, lambda),
                   gr = compute_gradient_function(Y, X, lambda),
                   upper = 0, method = "L-BFGS-B")
      opt <- c(opt, method = "L-BFGS-B")
    }
    if(method == "shooting"){
      opt <- solve_multivariate(Delta0, Y, X, lambda)
    }
    if(method == "shooting2"){
      opt <- solve_multivariate2(Delta0, Y, X, lambda)
    }
    if(method == "shooting3"){
      opt <- solve_multivariate3(Delta0, Y, X, lambda, ...)
    }

    return(as_shiftestim(listopt = opt, tree = tree, zscores = zscores,
                         lambda = lambda, alpha = alpha, covar_mat = covar_mat))

  } else {

    df_bic <- expand.grid(alpha = alpha, lambda = lambda)
    best_bic <- Inf
    best_model <- NULL
    i <- 1

    for (lam in lambda) {
      for (alp in alpha) {
        model <- estimate_shifts3(Delta0, zscores, tree, alpha = alp,
                                  lambda = lam, method = method, ...)
        df_bic[i, 3] <- model$bic
        if (model$bic < best_bic) {
          best_model <- model
          best_bic <- model$bic
        }
        i <- i + 1
      }
    }

    if(is.null(best_model)){
      warning("No model can be selectionned by BIC. Please consider ")
    }

    colnames(df_bic)[3] <- "bic"
    best_model$optim_info$bic_selection <- df_bic
    best_model$method <- paste(best_model$method, "with model selection")
    return(best_model)
  }
}
