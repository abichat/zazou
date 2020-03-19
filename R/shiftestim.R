
#' 'shiftestim' object
#'
#' @param listopt an output of 'estimate_shifts()'
#' @param tree tree
#' @param zscores vector of observed z-scores
#' @param lambda regularization parameter
#' @param alpha parameter
#'
#' @return a 'shiftestim' object
#' @export
as_shiftestim <- function(listopt, tree, zscores, lambda, alpha) {

  required_names <- c("par", "value", "method")

  if (!is.list(listopt) ||
      !all(required_names %in% names(listopt))) {
        stop("'listopt' must be the output of a optimisation function.")
  }

  zscores_est <- incidence_matrix(tree) %*% listopt$par$estimate
  zscores_est <- data.frame(leaf = rownames(zscores_est),
                            estimate = zscores_est[, 1],
                            stringsAsFactors = FALSE)
  rownames(zscores_est) <- NULL

  if(listopt$method %in% "desparsified lasso"){
    hciz <-
      size_half_confint_zscores(
        covariance_noise_mat = listopt$covariance_noise_matrix,
        incidence_mat = incidence_matrix(tree),
        hsigma = listopt$hsigma_scaledlasso,
        alpha_conf = listopt$alpha_confint
      )
    zscores_est$lower <- zscores_est$estimate - hciz
    zscores_est$upper <- zscores_est$estimate + hciz
  }


  obj <- list(zscores_obs = zscores,
              zscores_est = zscores_est,
              shift_est = listopt$par,
              objective_value = listopt$value, lambda = lambda,
              tree = tree, alpha = alpha,
              method = listopt$method,
              optim_info = listopt[setdiff(names(listopt), required_names)])

  ### Compute accessory slots --------------
  ## Binary tree
  obj$is_bin <- is.binary(tree)

  ## Parsimony score
  states <- obj$zscores_est$estimate
  names(states) <- obj$zscores_est$leaf
  obj$pars_score <- fitch(multi2di(obj$tree), states) # parsimony score

  ## Sigma
  h <- tree_height(obj$tree)
  obj$sigma <- sqrt(2 * obj$alpha) / (1 - exp(- 2 * obj$alpha * h))

  ## BIC & pBIC
  obj$bic <- bic(obs_zscores = obj$zscores_obs,
                 est_zscores = obj$zscores_est$estimate,
                 est_shifts = obj$shift_est$estimate, sigma = obj$sigma)
  obj$pbic <- pbic(obs_zscores = obj$zscores_obs,
                   est_zscores = obj$zscores_est$estimate,
                  est_shifts = obj$shift_est$estimate, sigma = obj$sigma,
                  alpha = alpha, tree = tree)

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

  txt_alpha <- paste0("Covariance matrix has been estimated from an OU",
                      " with alpha = ", round(x$alpha, digits),
                      " and sigma = ", round(x$sigma, digits), "")

  txt_tree1 <- paste0("Tree is", ifelse(x$is_bin, " ", " not "), "binary")
  tree <- x$tree

  txt_tree2 <- paste0("with ", length(x$tree$tip.label), " leafs and ",
                      nrow(x$tree$edge), " branches\n")

  zobs <- as.character(head(round(x$zscores_obs, digits), 10))
  zest <- as.character(head(round(x$zscores_est$estimate, digits), 10))
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
  cat(paste0("BIC: ", round(x$bic, digits), "\n"))
  cat(paste0("pBIC: ", round(x$pbic, digits), "\n"))
  cat("---\n")
  cat("Estimated shifts:", head(round(x$shift_est$estimate, digits), 10),
      "...\n")
  cat(norm0(x$shift_est$estimate), "shifts have been identified (ie",
      100 * round(norm0(x$shift_est$estimate, rev = TRUE, prop = TRUE), digits),
      "% of sparsity)\n")
  cat("A parsimonious solution would involve", x$pars_score, "shifts\n")
  cat("---\n")
  cat("Observed z-scores: ", zobs, dots_z)
  cat("Estimated z-scores:", zest, dots_z)
  cat(norm0(x$zscores_est$estimate),  "z-scores have been shifted (ie",
      100 * round(norm0(x$zscores_est$estimate, rev = TRUE, prop = TRUE),
                  digits),
      "% of sparsity)\n")
}


#' @rdname as_shiftestim
#' @inheritParams plot_shifts
#'
#' @export
plot.shiftestim <- function(x, digits = 3, ...){
  plot_shifts(tree = x$tree, shifts = x$shift_est$estimate,
              obs_scores = x$zscores_obs, est_scores = x$zscores_est$estimate,
              digits = digits, ...)
}


