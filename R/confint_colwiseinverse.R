#' Desparsified lasso based on column-wise inverse matrix
#'
#' @details If \code{gamma} is missing, it is chose with
#' \code{generate_gamma} whose \code{factor} argument could be passed to with
#' \code{...}.
#'
#' @param x A 'shiftpunct' object.
#' @param alpha_conf The confidence level.
#' @param ... Further arguments to be passed to or from other methods.
#' @inheritParams solve_colwiseinverse
#'
#' @return a list
#' @export
#'
confint_colwiseinverse <- function(x, alpha_conf = 0.05, gamma,
                                   silent_on_errors = TRUE, ...){
  stopifnot(inherits(x, "shiftpunct"))

  mat_covarOU <- covarianceOU_matrix(x$tree, x$alpha)
  mat_incidence <- incidence_matrix(x$tree)
  R <- inverse_sqrt(mat_covarOU)
  Y <- R %*% x$zscores_obs
  X <- R %*% mat_incidence
  # hsigma <- x$optim_info$sigma_scaledlasso

  XTXn <- crossprod(X) / nrow(X)

  if(missing(gamma)){
    gamma <- generate_gamma(X, ...)
  }

  M <- try(solve_colwiseinverse(A = XTXn, gamma = gamma,
                                silent_on_errors = silent_on_errors, ...),
           silent = silent_on_errors)

  if(!is.matrix(M)){
    warning("Constrains are not feasible. M is set to identity matrix.")
    M <- diag(ncol(X))
  }


  new_beta <- update_beta_colwiseinverse(X = X, y = Y, beta_init = x$shifts_est,
                                         colwiseinverse = M)

  # tau <- noise_factor_colwiseinverse(X, colwiseinverse = M, XTXn = XTXn)
  mat_covar_noise <-
    covariance_noise_matrix_colwiseinverse(X, colwiseinverse = M, XTXn = XTXn)

  # shifts_est <- df_confint_pvalue(estimate = new_beta, sigma = hsigma,
  #                                tau = tau, alpha_conf = alpha_conf)
  #
  # zscores_est <- df_conf_leaves(shifts = shifts_est$estimate,
  #                               covariance_noise_mat = mat_covar_noise,
  #                               mat_incidence = mat_incidence, sigma = hsigma,
  #                               alpha_conf = alpha_conf)


  list(shifts_est = data.frame(estimate = new_beta),
       zscores_est = NA, covariance_noise_matrix = mat_covar_noise,
       method = "colwise inverse",
       colwiseinverse = M, gamma = gamma)
}


#' @rdname update_beta_scoresystem
#' @param colwiseinverse The columnwise inverse.
#'
#' @export
update_beta_colwiseinverse <- function(X, y, beta_init, colwiseinverse){
  res <- y - X %*% beta_init
  correction <- colwiseinverse %*% t(X) %*% res / nrow(X)
  as.numeric(beta_init + correction)
}

#' @rdname noise_factor_scoresystem
#' @inheritParams update_beta_colwiseinverse
#' @param XTXn Matrix.
#'
#' @export
noise_factor_colwiseinverse <- function(X, colwiseinverse, XTXn) {
  MXTXxM <- colwiseinverse %*% XTXn %*% t(colwiseinverse)
  sqrt(diag(MXTXxM) / nrow(X))
}

#' @rdname covariance_noise_matrix_scoresystem
#' @inheritParams noise_factor_colwiseinverse
#' @export
#'
covariance_noise_matrix_colwiseinverse <- function(X, colwiseinverse, XTXn){
  colwiseinverse %*% XTXn %*% t(colwiseinverse) / nrow(X)
}
