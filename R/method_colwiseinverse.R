#' Solve column-wise inverse
#'
#' @param A A semi positive definite matrix
#' @param gamma Numeric. Non-negative
#' @param ntry_max Integer. Maximum umber of try for each column.
#' @param silent_on_errors Logical, default to TRUE.
#' @param silent_on_try Logical, default to TRUE.
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return The column-wise, inverse, same size as \code{A}.
#' @export
#'
solve_colwiseinverse <- function(A, gamma, ntry_max = 500000,
                                 silent_on_tries = TRUE, silent_on_errors = TRUE,
                                 ...){
  ## Validate that A is semi definite positive
  stopifnot(
    isSymmetric(A),
    all(eigen(A, symmetric = TRUE, only.values = TRUE)$values >= 0)
  )

  dim <- ncol(A)
  gen_inv <- solve(A + diag(gamma, nrow = dim))
  M <- matrix(NA, nrow = dim, ncol = dim)

  for(col in seq_len(dim)){

    # cat("\n###############\n### col =", col, "###\n###############\n")

    col0 <- gen_inv[, col]

    new_col <- NULL
    ntry <- 0

    while(!is.numeric(new_col) && ntry < ntry_max){

      # cat("\n####### ntry for col", col, "=", ntry, "#######\n\n")
      col0 <- col0 + rnorm(dim, sd = 1 / sqrt(dim))

      new_col <- try(solve_colwiseinverse_col(col, A, gamma,
                                              m0 = col0),
                     silent = silent_on_errors)
      ntry <- ntry + 1
    }

    try_ies <- ifelse(ntry >= 2, "tries", "try")

    if (is.numeric(new_col)) {
      M[, col] <- new_col
      if (!silent_on_tries)
        cat("Column ", col, " succeeded after ", ntry, " ", try_ies, ".\n",
            sep = "")
    } else {
      if (!silent_on_tries)
        cat("Column ", col, " failed after ", ntry, " ", try_ies, ".\n",
            sep = "")
      stop(paste0("Algorithm fails to converge for column ", col, "."))

    }


  }
  return(M)
}


#' Compute a column for the column-wise inverse.
#'
#' @param col The column number.
#' @param m0 Startup column.
#' @inheritParams solve_colwiseinverse
#'
#' @importFrom stats rbinom
#'
#' @return The column, wile respecting constrains.
#' @export
#'
solve_colwiseinverse_col <- function(col, A, gamma, m0){
  dim <- ncol(A)

  if(missing(m0)){
    # m <- rep(0, dim)
    # m <- rnorm(dim, sd = 1 / sqrt(dim))
    # m <- rep(1, dim)
    m <- sample(x = c(-1, 1), size = dim, replace = TRUE) * rnorm(dim, mean = 1,
                                                 sd = 1 / sqrt(dim))
  } else {
    m <- m0
  }

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
    stop("Constrains are not feasible.")
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






#' Compute a column for the column-wise inverse.
#'
#' @param col The column number.
#' @param m Startup column.
#' @param max_it Integer. Maximum number of iterations
#' @inheritParams solve_colwiseinverse
#'
#' @return The column, while respecting constrains.
#' @export
#'
fast_solve_colwiseinverse_col <- function(col, A, gamma, m, max_it = 5000) {

  ### Bookkeeping and transformation (should be move to outer loop)
  dim <- ncol(A)
  svd <- svd(A, nv = 0)
  U <- svd$u; d <- svd$d
  B <- U %*% diag(d)

  ##  Initialize m
  if (missing(m)) m <- rnorm(
    n = dim,
    mean = c(-1, 1)[sample.int(2, dim, TRUE)],
    sd = 1 / sqrt(dim)
  )
  ## constraint vectors: Bm - e_i
  constraint <- B %*% m - rep(c(0, 1, 0), times = c(col - 1, 1, dim - col))

  if (any(abs(constraint) > gamma)) warning("The starting point is not in the feasible set. Updates may be meaningless.")

  ## Objective function
  fn_obj <- function(x) { sum(svd$d*x^2) }
  obj_vals <- fn_obj(m)

  # Helper functions:
  ## update coord only has side effects
  update_coord <- function(coord){
    mi <- m[coord]
    b <- B[, coord]
    c <- constraint - mi * b ## B[ , -coord] %*% m[-coord] - e_i

    # update cell mi
    mi <- try(update_cell(c, b, gamma), silent = TRUE)
    if (is.numeric(mi)) {
      # update m and constraints only if mi is numeric
      m[coord] <<- mi
      constraint <<- c + mi * b
    }
  }

  it <- 1
  eps <- 10 ^ -8
  progress <- +Inf

  while (it < max_it && progress > eps && obj_vals > 0) {
    ## Update coordinates in random order (rather than random coordinates)
    coord_order <- sample.int(dim)
    for (coord in coord_order) {
      # cat(coord, sep = "\n")
      update_coord(coord)
      # cat(paste(m, collapse = ","), sep = "\n")
    }
    ## Store current objective value and compute progress
    it <- it + 1
    new_obj_vals <- fn_obj(m)
    progress <- abs(new_obj_vals - obj_vals) / obj_vals
    obj_vals <- new_obj_vals
  }

  ## Checkfor problems
  if (it == max_it) warning("Convergence not reached")

  ## return solution
  t(U) %*% m
}


update_cell <- function(c, b, gamma) {
  ## c = c - e
  ## Rounding to avoid numerical errors
  b[abs(b) < 10 * .Machine$double.eps] <- 0
  active_set <- b != 0
  ## Upper bound of feasible set: min_{i: b[i] != 0} (c[i]+gamma)/b[i]
  if (any( (-c + gamma < - sqrt(.Machine$double.eps)) & b == 0)) {
    upper_bound <- -Inf
  } else {
    upper_bound <- min( (-c + gamma)[active_set] / b[active_set] )
  }
  ## Lower bound of feasible set: max_{i: b[i] != 0} (c[i]-gamma)/b[i]
  if (any( (-c - gamma > sqrt(.Machine$double.eps)) & b == 0)) {
    lower_bound <- Inf
  } else {
    lower_bound <- max( (-c - gamma)[active_set] / b[active_set] )
  }
  ## Check that feasible set (u + v * beta <= 0) is not empty
  if (upper_bound - lower_bound < sqrt(.Machine$double.eps)) {
    stop("The constraint is not feasible. Consider changing the constraint.")
  }
  ## Project 0 on feasibility set and return result
  max(lower_bound, min(0, upper_bound))
}
