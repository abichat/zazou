solve_colwiseinverse <- function(A, gamma){
  dim <- ncol(A)
  M <- matrix(NA, nrow = dim, ncol = dim)
  for(col in seq_len(dim)){
    M[, col] <- solve_colwiseinverse_col(col, A, gamma)
  }
  return(M)
}


solve_colwiseinverse_col <- function(col, A, gamma){
  dim <- ncol(A)
  m <- rep(0, dim)
  for(iter in 1:1000){
    sampled_ind <- sample(dim)
    for(ind in sampled_ind){
      m[ind] <- solve_colwiseinverse_col_cell(m, ind, col, A, gamma)
    }
  }
  return(m)
}

solve_colwiseinverse_col_cell <- function(m, ind, col, A, gamma){
  vec_coef_p2 <- coef_p2(m, ind, A)
  argmin_unconstrained <- minimun_p2(vec_coef_p2)
  bounds <- compute_bounds(m, ind, col, A, gamma)
  if(length(argmin_unconstrained) == 1){
    argmin <- min(bounds$upper_bound,
                  max(bounds$lower_bound, argmin_unconstrained))
  } else {
    argmin <- argmin_between_two(m, ind, A,
                                 bounds$lower_bound,
                                 bounds$upper_bound)
  }
  return(argmin)
}

coef_p2 <- function(m, ind, A){
  a <- A[ind, ind]
  b <- sum(m[-ind] * A[-ind, ind]) + sum(m[-ind] * A[ind, -ind])
  c <- as.vector(m[-ind] %*% A[-ind, -ind] %*% m[-ind])
  return(list(a = a, b = b, c = c))
}

minimun_p2 <- function(coef) {
  if (coef$a > 0) {
    return(-coef$b / (2 * coef$a))
  } else if (coef$a == 0) {
    return(0)
  } else {
    return(c(-Inf, Inf))
  }
}

compute_bounds <- function(m, ind, col, A, gamma){
  dim <- length(m)
  Am_minus_ind <- A[, -ind] %*% m[- ind]
  # Am_ind <- A[, ind] * m[ind]
  ecol <- rep(0, dim)
  ecol[col] <- 1
  bounds <- matrix(c(-gamma, gamma), byrow = TRUE,
                   ncol = 2, nrow = dim)
  bounds <- (bounds + ecol - as.numeric(Am_minus_ind)) / A[, ind]
  lb <- max(bounds[, 1])
  up <- min(bounds[, 2])
  if(lb > up){
    stop("Constrains are not feasible")
  }
  return(list(lower_bound = lb, upper_bound = ub))
}

argmin_between_two <- function(m, ind, A, prop1, prop2){
  m_1 <- m_2 <- m
  m_1[ind] <- prop1
  m_2[ind] <- prop2
  mAm_1 <- m_1 %*% A %*% m_1
  mAm_2 <- m_2 %*% A %*% m_2
  if(mAm_1 < mAm_2) {
    return(prop1)
  } else {
    return(prop2)
  }
}


