
#' Score system
#'
#' Compute the score system of a matrix.
#'
#' @param X Matrix to pseudo orthogonalize of size m*(n+m).
#'
#' @return The matrix of score system, same dimension as X.
#' @export
calculate_Z <- function(X){

  lambdas <- get_lambda_sequence(X = X)
  best_lambda <- choose_best_lambda(lambdas = lambdas, X = X)$lambda.min
  nodewiselasso_out <- calculate_Z_with_best_lambda(X = X, lambda = best_lambda)

  nodewiselasso_out$Z
}


#' Vector of lambdas to test
#'
#' @inheritParams calculate_Z
#'
#' @importFrom glmnet glmnet
#' @importFrom stats quantile
#' @return 100 values for lambdas
get_lambda_sequence <- function(X) {
  nlambda <- 100
  p <- ncol(X)
  lambdas <- c() # Unknow length for glmnet()$lambda
  for (c in 1:p) {
    lambdas <- c(lambdas, glmnet(X[,-c], X[, c])$lambda)
  }
  lambdas <- quantile(lambdas, probs = seq(0, 1, length.out = nlambda))
  lambdas <- sort(lambdas, decreasing = TRUE)
  return(lambdas)
}

#' Choose best lambda by cross-validation
#'
#' @inheritParams calculate_Z
#' @param lambdas Numeric vector of lambda to test.
#'
#' @return The best lambda, which minimize RMSE prediction
choose_best_lambda <- function(lambdas, X){
  K <- 10
  n <- nrow(X)
  p <- ncol(X)
  l <- length(lambdas)

  ## Based on code from cv.glmnet for sampling the data
  dataselects <- sample(rep(1:K, length = n))

  totalerr <- mapply(rmse_glmnet,
                     c = 1:p,
                     K = K,
                     dataselects = list(dataselects = dataselects),
                     X = list(X = X),
                     lambdas = list(lambdas = lambdas),
                     SIMPLIFY = FALSE)

  err.array  <- array(unlist(totalerr), dim = c(length(lambdas), K, p))
  err.mean   <- apply(err.array, 1, mean) ## 1 mean for each lambda

  pos.min    <- which.min(err.mean)
  lambda.min <- lambdas[pos.min]

  list(lambda.min = lambda.min)
}

#' Compute prediction errors
#'
#' @param c Column to use.
#' @param K Number of batches.
#' @param dataselects Batches.
#' @inheritParams choose_best_lambda
#'
#' @return Matrix of errors per lambda and batch.
#'
#' @importFrom glmnet glmnet
#' @importFrom stats predict
rmse_glmnet <- function(c, K, dataselects, X, lambdas) {
  totalerr <- matrix(nrow = length(lambdas), ncol = K)

  for(i in 1:K){ ## loop over the test sets
    whichj <- dataselects == i ##the test part of the data

    if(var(X[!whichj, c, drop = FALSE])>0){
      glmnetfit <- glmnet(x = X[!whichj, -c, drop = FALSE],
                                  y = X[!whichj,  c, drop = FALSE],
                                  lambda = lambdas)
      predictions  <- predict(glmnetfit, newx = X[whichj, -c, drop = FALSE],
                              s = lambdas)
      totalerr[, i] <- apply((X[whichj, c] - predictions)^2, 2, mean)
    } else {
      totalerr[, i] <- mean((X[whichj, c] - mean(X[!whichj, c]))^2)
    }
  }

  totalerr
}


#' Compute the score system once lambda is chosen
#'
#' @inheritParams choose_best_lambda
#' @param lambda Numeric.
#'
#' @return The score system.
calculate_Z_with_best_lambda <- function(X, lambda) {
  n <- nrow(X)
  p <- ncol(X)
  # Z <- matrix(numeric(n * p), n)
  Z <- mapply(calculate_Z_column_with_best_lambda, i = 1:p,
              X = list(X = X), lambda = lambda)
  Z <- rescale_score_system(Z, X)
  return(Z)
}


#' Compute a column of the score system
#'
#' @param i Integer. The column to use.
#' @inheritParams calculate_Z_with_best_lambda
#'
#' @return A column of the score system.
#'
#' @importFrom stats predict
calculate_Z_column_with_best_lambda <- function(i, X, lambda) {
  glmnetfit  <- glmnet(X[, -i], X[, i])
  prediction <- predict(glmnetfit, X[, -i], s = lambda)
  return(X[, i] - prediction)
}

#' Rescale the score system
#'
#' @param Z The non-scaled score system.
#' @inheritParams calculate_Z
#'
#' @return The scaled score system.
rescale_score_system <- function(Z, X) {
  scaleZ <- diag(crossprod(Z, X)) / nrow(X)
  Z      <- scale(Z, center = FALSE, scale = scaleZ)
  return(list(Z = Z, scaleZ = scaleZ))
}
