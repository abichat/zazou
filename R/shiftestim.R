
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
  zscores_est <- zscores_est[, 1]

  obj <- list(zscores_obs = zscores,
              zscores_est = zscores_est,
              shift_est = listopt$par,
              objective_value = listopt$value, lambda = lambda,
              tree = tree, alpha = alpha, covar_mat = covar_mat,
              method = listopt$method,
              optim_info = listopt[setdiff(names(listopt), required_names)])

  ### Compute accessory slots --------------
  ## Binary tree
  if (is.binary(tree)){
    obj$is_bin <- TRUE
  } else {
    obj$is_bin <- FALSE
  }
  ## Parsimony score
  obj$pars_score <- fitch(multi2di(obj$tree), obj$zscores_est) # parsimony score

  ## Sigma
  if (is.null(obj$alpha)){
    obj$sigma <- NULL
    obj$bic <- NULL
  } else {
    h <- tree_height(obj$tree)
    obj$sigma <- sqrt(2 * obj$alpha) / (1 - exp(- 2 * obj$alpha * h))
    obj$bic <- bic(obs_zscores = obj$zscores_obs, est_zscores = obj$zscores_est,
                   est_shifts = obj$shift_est, sigma = obj$sigma)
  }

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
#' @importFrom stringr str_pad
#' @importFrom ape is.binary multi2di
#'
#' @export
print.shiftestim <- function(x, digits = 3, ...){

  if(is.null(x$alpha)){
    txt_alpha <- "Covariance matrix has been manually specified."
    txt_bic <- ""
  } else {
    txt_alpha <- paste0("Covariance matrix has been estimated from an OU",
                        " with alpha = ", round(x$alpha, digits),
                        " and sigma = ", round(x$sigma, digits), "")
    txt_bic <- paste0("BIC: ", round(x$bic, digits), "\n")
  }

  txt_tree1 <- paste0("Tree is", ifelse(x$is_bin, " ", " not "), "binary")
  tree <- x$tree

  txt_tree2 <- paste0("with ", length(x$tree$tip.label), " leafs and ",
                      nrow(x$tree$edge), " branches\n")

  zobs <- as.character(head(round(x$zscores_obs, digits), 10))
  zest <- as.character(head(round(x$zscores_est, digits), 10))
  nchar <- pmax(nchar(zobs), nchar(zest))

  zobs <- str_pad(zobs, width = nchar, side = "left")
  zest <- str_pad(zest, width = nchar, side = "left")

  if(length(x$zscores_obs) <= 10){
    dots_z <- "\n"
  } else {
    dots_z <- "...\n"
  }

  cat(txt_tree1, txt_tree2)
  cat(txt_alpha, "\n")
  cat("---\n")
  cat("Optimisation algorithm: ", x$method, "\n", sep = "")
  cat("Regularization parameter: lambda =", round(x$lambda, digits), "\n")
  cat("Objective value: ", round(x$objective_value, digits), "\n", sep = "")
  cat(txt_bic)
  cat("---\n")
  cat("Estimated shifts:", head(round(x$shift_est, digits), 10), "...\n")
  cat(sum(x$shift_est != 0), "shifts have been identified (ie",
      100 * round(mean(x$shift_est == 0), digits), "% of sparsity)\n")
  cat("A parsimonious solution would involve", x$pars_score, "shifts\n")
  cat("---\n")
  cat("Observed z-scores: ", zobs, dots_z)
  cat("Estimated z-scores:", zest, dots_z)
  cat(sum(x$zscores_est != 0), "z-scores have been shifted (ie",
      100 * round(mean(x$zscores_est == 0), digits), "% of sparsity)\n")
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


