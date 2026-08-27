# ============================================================================
# Contrast matrices
# ============================================================================

# Build general contrast matrix based on the Jacobian of the working model, the contrast matrix A, and the covariance matrix Sigma.
build_contrast_matrix <- function(jacobian, A, Sigma) {
  # B_contrast <- t(A) %*% solve(A %*% Sigma %*% t(A)) %*% A %*% jacobian
  B_contrast <- t(A) %*% MASS::ginv(A %*% Sigma %*% t(A)) %*% A %*% jacobian
  Matrix::t(B_contrast)
}


# Contrast matrix for omnibus test with a single outcome. The contrast matrix is of the form:
# | 0 -1 0 ... 0 | 0 1 0 ... 0 |
# | 0 0 -1 ... 0 | 0 0 1 ... 0 |
# | 0 0 0 ... -1 | 0 0 0 ... 1 |
build_omnibus_contrast_single_outcome <- function(K) {
  contrast <- matrix(0, nrow = K, ncol = 2 * (K + 1))
  contrast[, 2:(K + 1)] <- diag(-1, K)
  contrast[, (K + 3):(2 * (K + 1))] <- diag(1, K)
  contrast
}

# Contrast matrix for omnibus test with multiple outcomes.
build_omnibus_contrast_multi_outcome <- function(J, K, times = NULL) {
  if (!is.null(times)) {
    if (!is.list(times)) {
      stop("times must be a list of numeric vectors.")
    }
    # number of outcomes
    J <- length(times)
    no_of_measurements <- purrr::map_dbl(times, ~ length(.x)) %>%
      sum()
    contrast_matrix <- matrix(0, nrow = no_of_measurements - J, ncol = 2 * no_of_measurements)
    
    for (outcome_idx in seq_len(J)) {
      K <- length(times[[outcome_idx]]) - 1
      if (outcome_idx == 1) {
        start_row <- 1
        end_row <- K
        start_col <- 1
        end_col <- 2 * (K + 1)
      } else {
        start_row <- sum(purrr::map_dbl(times[1:(outcome_idx - 1)], ~ length(.x) - 1)) + 1
        end_row <- start_row + K - 1
        start_col <- sum(purrr::map_dbl(times[1:(outcome_idx - 1)], ~ 2 * length(.x))) + 1
        end_col <- start_col + 2 * (K + 1) - 1
      }
      
      contrast_matrix[start_row:end_row, start_col:end_col] <- build_omnibus_contrast_single_outcome(K)
    }
    
    return(contrast_matrix)
  } else {
    diag(1, J) %x% build_omnibus_contrast_single_outcome(K)
  }
}

build_summing_contrast_multi_outcome <- function(J, K, times = NULL) {
  # The time points in times should agree; otherwise, summing across outcomes is not
  # meaningful.
  if (!is.null(times)) {
    if (!is.list(times)) {
      stop("times must be a list of numeric vectors.")
    }
    # number of outcomes
    J <- length(times)
    no_of_measurements_vec <- purrr::map_dbl(times, ~ length(.x))
    if (!all(no_of_measurements_vec == no_of_measurements_vec[1])) {
      stop("All elements of times must have the same length.")
    }
    K = no_of_measurements_vec[1] - 1
  }
  # Create a contrast matrix that sums the treatment effects across all outcomes. This is done by summing the
  # rows of the omnibus contrast matrix for each time point across all outcomes. The resulting contrast matrix will
  # have K rows (one for each time point) and 2 * J * (K + 1) columns (one for each outcome and time point).
  omnibus_contrast <- build_omnibus_contrast_multi_outcome(J, K)
  # Matrix to sum the treatment effects across all outcomes for each time point.
  sum_matrix <- matrix(0, nrow = K, ncol = J * K)
  for (time_idx in seq_len(K)) {
    sum_matrix[time_idx, seq(time_idx, J * K, by = K)] <- 1
  }
  # Multiply the omnibus contrast matrix by the sum matrix to get the summing contrast matrix.
  summing_contrast <- sum_matrix %*% omnibus_contrast
  summing_contrast
}

# Contrast matrix for linear working model with multiple outcomes.
build_linear_contrast_multi_outcome <- function(J, K, times, Sigma) {
  if (!is.list(times) & !is.numeric(times)) {
    stop("times must be a numeric vector or a list of numeric vectors.")
  }
  if (is.numeric(times)) {
    times <- rep(list(times), J)
  }
  # number of outcomes
  J <- length(times)
  
  # number of measurements
  no_of_measurements_vec <- purrr::map_dbl(times, ~ length(.x))
  no_of_measurements <- sum(no_of_measurements_vec)
  
  jacobian <- matrix(0, nrow = J, ncol = 2 * no_of_measurements)
  for (outcome_idx in seq_len(J)) {
    K <- no_of_measurements_vec[outcome_idx] - 1
    if (outcome_idx == 1) {
      start_col <- 1
    } else {
      start_col <- sum(purrr::map_dbl(times[1:(outcome_idx - 1)], ~ length(.x) * 2)) + 1
    }
    end_col <- start_col + K
    
    jacobian[outcome_idx, start_col:end_col] <- times[[outcome_idx]]
  }
  jacobian = t(jacobian)
  
  A_contrast <- build_omnibus_contrast_multi_outcome(times = times)
  
  build_contrast_matrix(jacobian, A_contrast, Sigma)
}

# Contrast matrix for linear working model with multiple outcomes, assuming a
# common slope across outcomes.
build_linear_contrast_common <- function(times, Sigma) {
  if (!is.list(times)) {
    stop("times must be a list of numeric vectors.")
  }
  # number of outcomes
  J <- length(times)
  
  # number of measurements
  no_of_measurements_vec <- purrr::map_dbl(times, ~ length(.x))
  no_of_measurements <- sum(no_of_measurements_vec)
  
  jacobian <- matrix(0, nrow = 1, ncol = 2 * no_of_measurements)
  for (outcome_idx in seq_len(J)) {
    K <- no_of_measurements_vec[outcome_idx] - 1
    if (outcome_idx == 1) {
      start_col <- 1
    } else {
      start_col <- sum(purrr::map_dbl(times[1:(outcome_idx - 1)], ~ length(.x) * 2)) + 1
    }
    end_col <- start_col + K
    
    jacobian[1, start_col:end_col] <- times[[outcome_idx]]
  }
  jacobian = t(jacobian)
  
  A_contrast <- build_omnibus_contrast_multi_outcome(times = times)
  
  build_contrast_matrix(jacobian, A_contrast, Sigma)
}