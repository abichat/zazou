#' Desparsified lasso based on column-wise inverse matrix
#'
#' @param x a 'shiftestim' object.
#' @param alpha_conf the confidence level.
#' @param ... further arguments to be passed to or from other methods.
#'
#' @return a list
#' @export
#'
confint_colwiseinverse <- function(x, alpha_conf = 0.05, ...){
  stopifnot(inherits(x, "shiftestim"))

  mat_covar <- covariance_matrix(x$tree, x$alpha)
  mat_incidence <- incidence_matrix(x$tree)
  R <- inverse_sqrt(mat_covar)
  Y <- R %*% x$zscores_obs
  X <- R %*% mat_incidence
  hsigma <- x$optim_info$sigma_scaledlasso

  XTXn <- crossprod(X) / nrow(X)

  gamma <- 2 * sqrt(log(ncol(X))/nrow(X))

  M <- try(solve_colwiseinverse(XTXn, gamma), silent = TRUE)
  # Error in if (lb > up) { : valeur manquante là où TRUE / FALSE est requis ??
  i <- 0

  while(!is.matrix(M) && i < 20){
    M <- try(solve_colwiseinverse(XTXn, gamma), silent = TRUE)
    i <- i + 1
  }

  if(!is.matrix(M)){
    warning("Constrains are not feasible. M is set to identity matrix.")
    M <- diag(ncol(X))
  }


  new_beta <- update_beta_colwiseinverse(X = X, y = Y, beta_init = x$shifts_est,
                                         colwiseinverse = M)

  tau <- noise_factor_colwiseinverse(X, colwiseinverse = M, XTXn = XTXn)

  shifts_est <- df_confint_pvalue(estimate = new_beta, sigma = hsigma,
                                 tau = tau, alpha_conf = alpha_conf)


  list(shifts_est = shifts_est, zscores_est = x$zscores_est,
       noise_factor = tau, alpha_conf = alpha_conf,
       method = "colwiseinverse")
}


update_beta_colwiseinverse <- function(X, y, beta_init, colwiseinverse){
  res <- y - X %*% beta_init
  correction <- colwiseinverse %*% t(X) %*% res / nrow(X)
  beta_init + correction
}

noise_factor_colwiseinverse <- function(X, colwiseinverse, XTXn) {
  MXTXxM <- colwiseinverse %*% XTXn %*% t(colwiseinverse)
  sqrt(diag(MXTXxM) / nrow(X))
}


