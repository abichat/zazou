#' Solve column-wise inverse
#'
#' @param A A square matrix to column-wise inverse.
#' @param gamma Non-negative.
#' @param silent_on_errors Logical, default to TRUE.
#'
#' @return The column-wise, inverse, same size as \code{A}.
#' @export
#'
solve_colwiseinverse <- function(A, gamma, silent_on_errors = TRUE){
  dim <- ncol(A)
  M <- matrix(NA, nrow = dim, ncol = dim)

  for(col in seq_len(dim)){

    # cat("\n###############\n### col =", col, "###\n###############\n")

    new_col <- try(solve_colwiseinverse_col(col, A, gamma), silent = silent_on_errors)

    ntry_max_per_column <- 500
    ntry <- 0

    while(!is.numeric(new_col) && ntry < ntry_max_per_column){

      # cat("\n####### ntry for col", col, "=", ntry, "#######\n\n")

      new_col <- try(solve_colwiseinverse_col(col, A, gamma), silent = silent_on_errors)
      ntry <- ntry + 1
    }

    if(is.numeric(new_col)){
      M[, col] <- new_col
    } else {
      stop("Constrains are not feasible.")
    }


  }
  return(M)
}


#' Compute a column for the column-wise inverse.
#'
#' @param col The column number.
#' @inheritParams solve_colwiseinverse
#'
#' @return The column, wile respecting constrains.
#' @export
#'
solve_colwiseinverse_col <- function(col, A, gamma){
  dim <- ncol(A)
  # m <- rep(0, dim)
  # m <- rnorm(dim, sd = 1 / sqrt(dim))
  # m <- rnorm(dim, sd = 1 / sqrt(dim))
  # m <- rep(1, dim)
  m <- 2 * (rbinom(dim, 1, 0.5) - 0.5) * rnorm(dim, mean = 1, sd = 1 / sqrt(dim))

  # cat("# First m =", m, "\n")


  max_it <- 100
  eps <- 10 ^ -8

  iter <- 0
  obj <- m %*% A %*% m
  progress <- +Inf

  while(iter < max_it && progress > eps){
    sampled_ind <- sample(dim)

    for(ind in sampled_ind){
      # cat("| > cell:", ind, "\n")
      m[ind] <- solve_colwiseinverse_col_cell(m, ind, col, A, gamma)
    }

    new_obj <- m %*% A %*% m
    progress <- abs((obj - new_obj) / obj)
    if(is.nan(progress)) progress <- 0
    obj <- new_obj
    iter <- iter + 1
    # cat("##--------------#\n")
    # cat("## iter:", iter, "| progress:", progress, "| obj:", obj, "\n")
    # cat("## m: ", m, "\n")
    # cat("##--------------#\n")
  }
  # cat("###################\n")
  return(m)
}

#' Update column cell
#'
#' @inheritParams solve_colwiseinverse_col
#' @param m The column whose a cell need to be updated.
#' @param ind The number of the cell in the column.
#'
#' @return The updated column.
#' @export
solve_colwiseinverse_col_cell <- function(m, ind, col, A, gamma){
  vec_coef_p2 <- coef_p2(m, ind, A)
  argmin_unconstrained <- minimun_p2(vec_coef_p2)
  # cat("|   >   argmin unconstrained:", argmin_unconstrained, "\n")
  bounds <- compute_bounds(m, ind, col, A, gamma)
  if(length(argmin_unconstrained) == 1){
    argmin <- min(bounds$upper_bound,
                  max(bounds$lower_bound, argmin_unconstrained))
  } else {
    argmin <- argmin_between_two(m, ind, A,
                                 bounds$lower_bound,
                                 bounds$upper_bound)
  }
  # cat("|   >   argmin:", argmin, "\n")
  if(abs(argmin) < 10^-15) argmin <- 0
  return(argmin)
}

#' Compute polynomial coefficients.
#'
#' Compute coefficients of the polynomial \eqn{m'Am} with the cell
#' \code{ind} of \code{m} as the indeterminate.
#'
#' @inheritParams solve_colwiseinverse_col_cell
#'
#' @return A list with \code{a}, \code{b} and \code{c} components.
#' @export
#' @examples
#'   m <- c(-1, 4)
#'   A <- matrix(c(3, 2, 2, 5), nrow = 2)
#'   coef_p2(A = A, m = m, ind = 1)
#'   coef_p2(A = A, m = m, ind = 2)
coef_p2 <- function(m, ind, A){
  a <- A[ind, ind]
  b <- sum(m[-ind] * A[-ind, ind]) + sum(m[-ind] * A[ind, -ind, drop = FALSE])
  c <- as.vector(m[-ind] %*% A[-ind, -ind, drop = FALSE] %*% m[-ind])
  return(list(a = a, b = b, c = c))
}

#' Compute the minimum of a second order polynomial
#'
#' @param coef A list with \code{a}, \code{b} and \code{c} components.
#'
#' @return A numeric of length one (eventually \code{Inf} or \code{-Inf})
#' if the polynomial if convex. \code{c(-Inf, Inf)} if the polynomial
#' is strictly concave.
#' @export
#' @seealso \code{\link{coef_p2}}
#' @examples
#' minimun_p2(list(a = 1, b = 2, c = 0))
#' minimun_p2(list(a = -1, b = 2, c = 0))
#' minimun_p2(list(a = 0, b = 2, c = 0))
#' minimun_p2(list(a = 0, b = 0, c = 3))
minimun_p2 <- function(coef) {
  if (coef$a > 0) {
    return(-coef$b / (2 * coef$a))
  } else if (coef$a == 0) {
    if(coef$b > 0){
      return(-Inf)
    } else if (coef$b < 0) {
      return(Inf)
    } else {
      return(0)
    }
  } else {
    return(c(-Inf, Inf))
  }
}

#' Compute constrains bounds
#'
#' @inheritParams solve_colwiseinverse_col_cell
#'
#' @return A list with \code{lower_bound} and \code{upper_bound} components.
#' @export
#'
compute_bounds <- function(m, ind, col, A, gamma){
  dim <- length(m)
  Am_minus_ind <- A[, -ind, drop = FALSE] %*% m[- ind, drop = FALSE]
  # Am_ind <- A[, ind] * m[ind]
  ecol <- rep(0, dim)
  ecol[col] <- 1
  bounds <- matrix(c(-gamma, gamma), byrow = TRUE,
                   ncol = 2, nrow = dim)
  bounds <- (bounds + ecol - as.numeric(Am_minus_ind)) / A[, ind]
  # correct the case when A[, ind] < 0
  for(i in seq_len(dim)){
    if(A[i, ind] < 0){
      bounds[i, ] <- rev(bounds[i, ])
    }
  }
  lb <- max(bounds[, 1])
  ub <- min(bounds[, 2])
  # cat("|   >   lower bound:", lb, "upper bound:", ub, "\n")
  if(lb > ub){
    stop("Constrains are not feasible")
  }
  return(list(lower_bound = lb, upper_bound = ub))
}

#' Find the argmin between two propositions
#'
#' @details The function to minimize is \eqn{m'Am}.
#'
#' @inheritParams solve_colwiseinverse_col_cell
#' @param prop1 First proposition for the \code{ind} cell in \code{m}.
#' @param prop2 Second proposition for the \code{ind} cell in \code{m}
#'
#' @return \code{prop1} or \code{prop2}.
#' @export
#'
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


