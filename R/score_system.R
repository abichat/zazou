#' Score system
#'
#' Compute the score system of a
#'
#' @param X A vector of size m*(n+m).
#' @param y A vector of size m.
#' @param beta_init Initial value of beta found with scaled lasso.
#' @param hsigma Estimate value of sigma, found with scaled lasso.
#'
#' @return The matrix of score system, same dimension as X.
#' @export
#'
#' @importFrom hdi lasso.proj
score_system <- function(X, y, beta_init, hsigma) {
  obj <-
    suppressWarnings(suppressMessages(
      hdi:::calculate.Z(x = X, parallel = FALSE, ncores = 1,
                        verbose = FALSE, Z = NULL, do.ZnZ = FALSE)
      # lasso.proj(x = X, y = y, betainit = beta_init,
      #            sigma = hsigma, return.Z = TRUE)
    ))
  sco <- obj$Z
  attr(sco, "scaled:scale") <- NULL
  sco
}


calculate_Z <- function(X){
  nodewiselasso.out <- score_nodewise_lasso(X = X)
  Z <- nodewiselasso.out$out$Z
  scaleZ <- nodewiselasso.out$out$scaleZ
}


score_nodewise_lasso <- function(X){
  lambdas <- hdi:::nodewise.getlambdasequence(x = X)

  cvlambdas <- cv_nodewise_bestlambda(lambdas = lambdas, X = X)
}

cv_nodewise_bestlambda <- function(lambdas, X){
  K <- 10
  n <- nrow(X)
  p <- ncol(X)
  l <- length(lambdas)

  ## Based on code from cv.glmnet for sampling the data
  dataselects <- sample(rep(1:K, length = n))

  totalerr <- mapply(cv_nodewise_err_unitfunction,
                     c = 1:p,
                     K = K,
                     dataselects = list(dataselects = dataselects),
                     X = list(X = X),
                     lambdas = list(lambdas = lambdas),
                     SIMPLIFY = FALSE)
}

cv_nodewise_err_unitfunction <- function(c, K, dataselects, X, lambdas) {
  cv_nodewise_totalerr(c = c,
                       K = K,
                       dataselects = dataselects,
                       X = X,
                       lambdas = lambdas)
}


cv_nodewise_totalerr <- function(c, K, dataselects, X, lambdas) {
  totalerr <- matrix(nrow = length(lambdas), ncol = K)

  for(i in 1:K){ ## loop over the test sets
    whichj <- dataselects == i ##the test part of the data

    glmnetfit <- glmnet::glmnet(x = X[!whichj,-c, drop = FALSE],
                        y = X[!whichj, c, drop = FALSE],
                        lambda = lambdas, alpha = 0)
    predictions  <- predict(glmnetfit, newx = X[whichj, -c, drop = FALSE],
                            s = lambdas)
    totalerr[, i] <- apply((X[whichj, c] - predictions)^2, 2, mean)
  }

  totalerr
}
