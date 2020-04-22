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

  mat_covarOU <- covarianceOU_matrix(x$tree, x$alpha)
  mat_incidence <- incidence_matrix(x$tree)
  R <- inverse_sqrt(mat_covarOU)
  Y <- R %*% x$zscores_obs
  X <- R %*% mat_incidence
  hsigma <- x$optim_info$sigma_scaledlasso

  XTXn <- crossprod(X) / nrow(X)

  gamma <- generate_gamma(X)

  M <- try(solve_colwiseinverse(XTXn, gamma), silent = TRUE)

  ntry_max_for_matrix <- 100
  ntry_mat <- 0

  while(!is.matrix(M) && ntry_mat < ntry_max_for_matrix){
    M <- try(solve_colwiseinverse(XTXn, gamma), silent = TRUE)
    ntry_mat <- ntry_mat + 1
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


