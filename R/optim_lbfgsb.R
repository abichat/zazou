solve_lbfgsb <- function(Delta0, X, Y, lambda,
                         constraint_type = c("beta", "yhat", "none"), ...){

  constraint_type <- match.arg(constraint_type)
  if(constraint_type == "yhat"){
    stop("The constraint 'yhat' is not available for L-BFGS-B solving.")
  }
  upbbound <- ifelse(constraint_type == "beta", 0, Inf)

  opt <- optim(par = Delta0,
               fn = compute_objective_function(Y, X, lambda, type = "lasso"),
               gr = compute_gradient_function(Y, X, lambda),
               upper = upbbound, method = "L-BFGS-B", ...)
  opt <- c(opt, method = "L-BFGS-B")
  opt$par <- data.frame(estimate = opt$par)
  opt
}
