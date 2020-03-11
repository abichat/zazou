# # devtools::load_all(".")
#
# data(chlamydiae)
#
# pval_obs <- test_kruskalwallis(chlamydiae$X, chlamydiae$Y)$p.value
# zsco_obs <- p2z(pval_obs)
#
# tree <- force_ultrametric(chlamydiae$tree)
# N_branch <- length(tree$edge.length)
#
# # estDL <- estimate_shifts(Delta0 = rep(0, N_branch), zscores = zsco_obs,
# #                          lambda = 1, tree = tree,
# #                          alpha = 1, method = "desparsifiedlasso")
#
# alp <- 1
#
#
# mat_covar <- covariance_matrix(tree, alp)
# incidence_mat <- incidence_matrix(tree)
# R <- inverse_sqrt(mat_covar)
# Y <- R %*% zsco_obs
# X <- R %*% incidence_mat
#
# # calculate_Z(X)
#
# lambdas <- hdi:::nodewise.getlambdasequence(x = X)
#
# cv_nodewise_bestlambda(lambdas = lambdas, X = X)
# calculate_Z(X)
#
#
# K <- 10
# n <- nrow(X)
# p <- ncol(X)
# l <- length(lambdas)
#
# ## Based on code from cv.glmnet for sampling the data
# set.seed(42)
# dataselects <- sample(rep(1:K, length = n))
#
# totalerr <- mapply(cv_nodewise_err_unitfunction,
#                    c = 1:p,
#                    K = K,
#                    dataselects = list(dataselects = dataselects),
#                    X = list(X = X),
#                    lambdas = list(lambdas = lambdas),
#                    SIMPLIFY = FALSE)
#
# for(i in 1:p){
#   print(i)
#   cv_nodewise_totalerr(c = i, K = K, dataselects = dataselects,
#                                X = X, lambdas = lambdas)
# }
#
# # i <- 17
# # cv_nodewise_totalerr(c = i, K = K, dataselects = dataselects,
# #                      X = X, lambdas = lambdas)
#
#
# c <- 17
#
# totalerr <- matrix(nrow = length(lambdas), ncol = K)
#
# for(i in 1:K){ ## loop over the test sets
#   print(i)
#   whichj <- dataselects == i ##the test part of the data
#
#   glmnetfit <- glmnet::glmnet(x = X[!whichj,-c, drop = FALSE],
#                               y = X[!whichj, c, drop = FALSE],
#                               lambda = lambdas, alpha = 1)
#   predictions  <- predict(glmnetfit, newx = X[whichj, -c, drop = FALSE],
#                           s = lambdas)
#   totalerr[, i] <- apply((X[whichj, c] - predictions)^2, 2, mean)
# }
#
# i <- 8
#
# X
#
# # glmnetfit <- glmnet::glmnet(x = X[!whichj,-c, drop = FALSE],
# #                             y = X[!whichj, c, drop = FALSE],
# #                             lambda = lambdas, alpha = 1)
#
# X[!whichj, c, drop = FALSE]
#
# which(whichj)
# which(X[, c] != 0)
#
#
# totalerr
#
# err_array <- array(unlist(totalerr), dim = c(length(lambdas), K, p))
# err_mean <- apply(err_array, 1, mean)
#
# which.min(err_mean)
# lambdas[which.min(err_mean)]
#
