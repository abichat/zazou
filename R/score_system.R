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
  nodewiselasso_out <- score_nodewise_lasso(X = X)
  Z <- nodewiselasso_out$out$Z
  scaleZ <- nodewiselasso_out$out$scaleZ
  list(Z = Z, scaleZ = scaleZ)
}


score_nodewise_lasso <- function(X){
  lambdas <- nodewise_getlambdasequence(X = X)

  cvlambdas <- cv_nodewise_bestlambda(lambdas = lambdas, X = X)

  bestlambda <- cvlambdas$lambda.min

  Z <- score_getZforlambda(x = X, lambda = bestlambda)

  out <- Z

  return_out <- list(out = out,
                     bestlambda = bestlambda)
  return(return_out)
}

nodewise_getlambdasequence <- function(X) {
  nlambda <- 100
  p <- ncol(X)
  lambdas <- c()
  for (c in 1:p) {
    lambdas <- c(lambdas, glmnet::glmnet(X[,-c], X[, c])$lambda)
  }
  lambdas <- quantile(lambdas, probs = seq(0, 1, length.out = nlambda))
  lambdas <- sort(lambdas, decreasing = TRUE)
  return(lambdas)
}

cv_nodewise_bestlambda <- function(lambdas, X){
  K <- 10
  n <- nrow(X)
  p <- ncol(X)
  l <- length(lambdas)

  ## Based on code from cv.glmnet for sampling the data=
  dataselects <- sample(rep(1:K, length = n))

  totalerr <- mapply(cv_nodewise_err_unitfunction,
                     c = 1:p,
                     K = K,
                     dataselects = list(dataselects = dataselects),
                     X = list(X = X),
                     lambdas = list(lambdas = lambdas),
                     SIMPLIFY = FALSE)

  err.array  <- array(unlist(totalerr), dim = c(length(lambdas), K, p))
  err.mean   <- apply(err.array, 1, mean) ## 1 mean for each lambda

  err.se     <- apply(apply(err.array, c(1, 2), mean), 1, sd) / sqrt(K)

  pos.min    <- which.min(err.mean)
  lambda.min <- lambdas[pos.min]

  stderr.lambda.min <- err.se[pos.min]

  list(lambda.min = lambda.min,
       lambda.1se = max(lambdas[err.mean < (min(err.mean) + stderr.lambda.min)]))
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

    if(var(X[!whichj, c, drop = FALSE])>0){
      glmnetfit <- glmnet::glmnet(x = X[!whichj, -c, drop = FALSE],
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


score_getZforlambda <- function(x, lambda) {
  n <- nrow(x)
  p <- ncol(x)
  Z <- matrix(numeric(n * p), n)
  Z <- mapply(score_getZforlambda_unitfunction, i = 1:p,
              x = list(x = x), lambda = lambda)
  Z <- score_rescale(Z, x)
  return(Z)
}


score_getZforlambda_unitfunction <- function(i, x, lambda) {
  glmnetfit  <- glmnet::glmnet(x[, -i], x[, i])
  prediction <- predict(glmnetfit, x[, -i], s = lambda)
  return(x[, i] - prediction)
}

score_rescale <- function(Z, x) {
  scaleZ <- diag(crossprod(Z, x)) / nrow(x)
  Z      <- scale(Z, center = FALSE, scale = scaleZ)
  return(list(Z = Z, scaleZ = scaleZ))
}
