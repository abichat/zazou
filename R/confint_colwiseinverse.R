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
  for(iter in 1:100){
    sampled_ind <- sample(dim)
    for(ind in sampled_ind){
      m[ind] <- solve_colwiseinverse_col_cell(m, ind, col, A, gamma)
    }
  }
  return(m)
}

solve_colwiseinverse_col_cell <- function(m, ind, col, A, gamma){
  vec_coef_p2 <- coef_p2(m, ind, A)
  min_unconstrained <- minimun_p2(vec_coef_p2)
  bounds <- compute_bounds(m, ind, col, A, gamma)

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
  return(list(lower_bound = max(bounds[, 1]),
              upper_bound = min(bounds[, 2])))
}


