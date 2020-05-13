#' Solve column-wise inverse
#'
#' @param A A semi positive definite matrix
#' @param gamma Numeric. Non-negative
#' @param ntry_max Integer. Maximum umber of try for each column.
#' @param silent_on_errors Logical, default to TRUE.
#' @param silent_on_tries Logical, default to TRUE.
#' @param ... Further arguments to be passed to or from other methods.
#'
#' @return The column-wise, inverse, same size as \code{A}.
#' @export
#'
solve_colwiseinverse <- function(A, gamma, ntry_max = 10000,
                                 silent_on_tries = TRUE,
                                 silent_on_errors = TRUE, ...){
  ## Validate that A is semi definite positive through SVD decompisition
  if (!isSymmetric(A)) stop("Matrix A should be symmetric.")
  dim <- ncol(A)
  svd <- svd(A, nu = dim, nv = 0)
  U <- svd$u
  d <- svd$d
  if (any(d < 0)) stop("Matrix A should be semi positive definite.")

  ## Compute generalized inverse through svd
  gen_inv <- t(U) %*% diag(1 / (d + gamma)) %*% U
  # gen_inv <- solve(A + diag(gamma, nrow = dim))
  M <- matrix(NA, nrow = dim, ncol = dim)

  for(col in seq_len(dim)){

    # cat("\n###############\n### col =", col, "###\n###############\n")

    col0 <- gen_inv[, col]

    new_col <- NULL
    ntry <- 0

    while(!is.numeric(new_col) && ntry < ntry_max){

      # cat("\n####### ntry for col", col, "=", ntry, "#######\n\n")
      col0 <- col0 #+ rnorm(dim, sd = 1 / sqrt(dim))

      new_col <- try(fast_solve_colwiseinverse_col(col = col, svdA = svd, A = A,
                                                   gamma = gamma, m0 = col0),
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
#' @param svdA The SVD decomposition of \code{A} as e list with \code{d} and
#' \code{u} components.
#' @param col The column number.
#' @param m0 Startup column.
#' @param max_it Integer. Maximum number of iterations.
#' @inheritParams solve_colwiseinverse
#'
#' @return The column, while respecting constrains.
#' @export
#'
fast_solve_colwiseinverse_col <- function(col, svdA, A, gamma, m0, max_it = 5000, ...) {

  ### Bookkeeping and transformation (should be move to outer loop)
  if(missing(svdA)){
    svdA <- svd(A, nv = 0)
  }
  U <- svdA$u
  d <- svdA$d
  dim <- ncol(U)
  B <- U * d[col(U)] ## equivalent to but faster than U %*% diag(d)
  e_i <- rep(c(0, 1, 0), times = c(col - 1, 1, dim - col))

  ##  Initialize m
  ## TO DO:
  ## - truncate d to d[1:K] where d[(K+1):n] = 0
  ## - keep only m_small = m[1:K] and optimize only on m_small
  ## - return m = c(m_small, rep(0, n-K))
  if (missing(m0)) {
    ## Not too costly initialization, should be in the feasible set
    ## If d is ill-conditioned, increase eigenvalues before inversion
    d_inv <- d
    if (d[length(d)] / d[1] < 1e-5){
      d_inv <- d + max(1e-5, max(d) / 1e5)
    }
    m <- (1/d_inv) * U[col, ]
    ## Equivalent to but faster than
    # m <- diag(1 / d_inv) %*% t(U) %*% e_i
  } else {
    m <- m0
  }
  ## constraint vectors: Bm - e_i
  constraint <- B %*% m - e_i

  # print(m)
  if (any(abs(constraint) > gamma)) {
    # cat("Starting point is not in the feasible set.\n")
    warning(paste("The starting point is not in the feasible set.",
                  "Updates may be meaningless."))
  }

  ## Objective function
  fn_obj <- function(x) { sum(d*x^2) }
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
    } else {
      n_skip <<- n_skip + 1
    }
    # cat("nskip = ", n_skip, "\n")
  }

  ## update best: Update coord that leads to smallest l2 norm
  ##
  update_smallest <- function() {
    best <- list(obj = Inf, i = NULL, mi = NA)
    ## Try all coordinates in turn and record the one with max decrease in the
    ## objective function
    for (i in 1:dim) {
      # cat("Small ", i, "\n")
      mi <- try(update_cell(constraint - m[i]*B[, i], B[, i], gamma),
                silent = TRUE)
      if (is.numeric(mi)) {
        m_temp <- m
        m_temp[i] <- mi
        obj <- sum((U %*% m_temp) ^ 2)
        if (obj <= best$obj) {
          best$obj <- obj
          best$i <- i
          best$mi <- mi
        }
      }
    }
    if (!is.null(best$i)) {
      # cat("We are in the feasible region !!!! ##############\n")
      constraint <<- constraint + (best$mi - m[best$i]) * B[, best$i]
      m[best$i] <<- best$mi
    } else {
      # cat("Smallest NULL: EXIT\n")
      stop("No cell could be updated for the first time in column ", col, ".")
    }
  }

  it <- 1
  eps <- 10 ^ -8
  progress <- +Inf

  while (it < max_it && progress > eps && obj_vals > 0) {
    # cat("Iteration", it, "\n")
    ## Try all coordinates and move along best one to avoid being stuck
    ## in local optimum
    if (it == 1) {
      ## Hacky: start with update leading to smallest l2 norm
      update_smallest()
    }
    ## Update coordinates in random order (rather than random coordinates)
    coord_order <- sample.int(dim)
    n_skip <- 0
    # cat("Reinitialization of nskip.\n")
    for (coord in coord_order) {
      # cat(coord, sep = "\n")
      update_coord(coord)
      # cat(paste(m, collapse = ","), sep = "\n")
    }
    ## Stop if no coord has been updated
    if(n_skip == dim){
      cat("## Exit ##\n")
      stop("No cell could be updated in column ", col, ".")
    }
    ## Store current objective value and compute progress
    it <- it + 1
    new_obj_vals <- fn_obj(m)
    progress <- abs(new_obj_vals - obj_vals) / obj_vals
    obj_vals <- new_obj_vals
    # cat("iteration", it, "done\n")
  }

  ## Checkfor problems
  if (it == max_it) warning("Convergence not reached for column ", col, ".")
  # cat("Column ", col, " succeeded after ", it, " iterations.\n", sep = "")
  ## return solution
  U %*% m
}


update_cell <- function(c, b, gamma) {
  ## c = c - e
  ## Rounding to avoid numerical errors
  b[abs(b) < 10 * .Machine$double.eps] <- 0
  active_set <- b != 0
  # cat("Active set:\n")
  # print(active_set)
  ## If |c_i| > gamma for any i such that b_i = 0, the feasible set is empty
  if (any(abs(c[!active_set]) > gamma)) {
    stop("Feasible set is empty")
  }
  if (any(active_set)) {
    ## Bounds of feasible sets are (-c_i +/- gamma) / b_i
    bound_1 <- (-c + gamma)[active_set] / b[active_set]
    bound_2 <- (-c - gamma)[active_set] / b[active_set]
    # cat("Bounds:\n")
    # print(bound_1)
    # print(bound_2)
    ## Min of all upper bounds
    upper_bound <- min(ifelse(bound_1 > bound_2, bound_1, bound_2))
    ## Max of all lower bounds
    lower_bound <- max(ifelse(bound_1 > bound_2, bound_2, bound_1))
    ## Check that feasible set is not empty
    if (upper_bound - lower_bound < sqrt(.Machine$double.eps)) {
      # cat("Column ", col, " failed after ", it, " iterations.\n", sep = "")
      stop("Feasible set is empty.")
    }
  } else {# b is null and c <= gamma
    # upper_bound <- Inf
    # lower_bound <- -Inf
    return(0)
  }
  ## Project 0 on feasibility set and return result
  max(lower_bound, min(0, upper_bound))
}


#' Find a feasible solution for a set of linear constraints
#'
#' @param B Matrix coding for a set of linear constraints
#' @param tol Tolerance allowed for the constraints
#' (mostly used to avoid floating point errors)
#' @inheritParams fast_solve_colwiseinverse_col
#'
#' @return A feasible point
#' @export
#'
#' @examples
#' B <- diag(1, 3)
#' find_feasible(B, 1, 0) ## c(1, 0, 0)
#' find_feasible(B, 1, 0.5) ## c(0.5, 0, 0)
#' B <- matrix(c(1, 1, 0, 0), 2)
#' \dontrun{
#' find_feasible(B, 1, gamma = 0.4) ## No solution for gamma lower than 0.5}
#' find_feasible(B, 1, gamma = 0.6) ## Plenty of solutions for gamma higher than 0.5
find_feasible <- function(B, col, gamma, m0 = rep(0, ncol(B)), max_it = 1e5, tol = 1e-14) {
  ## bookkeeping variables
  it <- 1
  n <- nrow(B) ## Ensures that procedures works for non square matrices
  e_i <- rep.int(c(0, 1, 0), times = c(col - 1, 1 , n - col))
  constraint <- B %*% m0
  B_normed <- B / rowSums(B^2) ## B[j, ] / \| B[j, ] \|_2^2
  B_normed[rowSums(B^2) == 0, ] <- 0 ## By construction, if B[j, ] is a null row, set B_normed[j, ] to a null row
  BtB_normed <- tcrossprod(B, B_normed) ## B %*% t(B_normed)

  ## Helper functions (called only for their side effects)
  is_feasible <- function() {
    all(abs(constraint - e_i) <= gamma + tol)
  }
  update_m0 <- function() {  ## used only for its side effect
    ## find most violated constraint
    j <- which.max(abs(constraint - e_i))
    ## Compute update for x
    step_length <- abs(constraint[j] - e_i[j]) - gamma
    direction <- -sign(constraint[j] - e_i[j])
    update <- direction * step_length
    ## Update m0 and constraint
    m0 <<- m0 + update * B_normed[j, ]
    constraint <<- constraint + update * BtB_normed[, j] ## BtB_normed[, j] = B %*% t(B_normed[j, ])
  }

  while (it <= max_it & !is_feasible()) {
    update_m0()
    it <- it + 1
  }

  if (!is_feasible()) warning(paste("No feasible solution found after", max_it, "iterations. Aborting"))

  return(m0)
}
