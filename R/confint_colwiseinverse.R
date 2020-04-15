solve_colwiseinverse <- function(A, gamma){
  dim <- ncol(A)
  M <- matrix(NA, nrow = dim, ncol = dim)
  for(j in seq_len(dim)){
    M[, j] <- solve_colwiseinverse_col(A, j, gamma)
  }
  return(M)
}


solve_colwiseinverse_col <- function(A, j, gamma){
  dim <- ncol(A)
  m <- rep(0, dim)
  for(iter in 1:100){
    sampled_ind <- sample(dim)
    for(ind in sampled_ind){
      m[ind] <- solve_colwiseinverse_col_cell(A, m, ind, gamma)
    }
  }
  retrun(m)
}

solve_colwiseinverse_col_cell <- function(A, m, ind, gamma){
  vec_coef_p2 <- coef_p2(A, m, ind)
}

coef_p2 <- function(A, m, ind){
  a <- A[ind, ind]
  b <- sum(m[-ind] * A[-ind, ind]) + sum(m[-ind] * A[ind, -ind])
  c <- m[-ind] %*% A[-ind, -ind] %*% m[-ind]
  return(c(a = a, b = b, c = c))
}


