# ============================================================================
# Contrast matrices
# ============================================================================

#' Build a plug-in contrast matrix
#'
#' Constructs the contrast matrix `B_n` used in the targeted test from a
#' working model's Jacobian, a user-specified contrast matrix `A`, and the
#' first-stage covariance matrix `Sigma`.
#'
#' @param jacobian numeric matrix, the Jacobian of the working model's mean
#'   function evaluated at the null estimate.
#' @param A numeric matrix, the contrast matrix defining the omnibus null
#'   hypothesis (see [build_omnibus_contrast_multi_outcome()]).
#' @param Sigma numeric matrix, the first-stage covariance estimate.
#'
#' @returns numeric matrix, the plug-in contrast matrix `B_n`.
#' @export
build_contrast_matrix <- function(jacobian, A, Sigma) {
  # B_contrast <- t(A) %*% solve(A %*% Sigma %*% t(A)) %*% A %*% jacobian
  B_contrast <- t(A) %*% MASS::ginv(A %*% Sigma %*% t(A)) %*% A %*% jacobian
  Matrix::t(B_contrast)
}


#' Omnibus contrast matrix for a single outcome
#'
#' Builds the contrast matrix testing equality of the control and
#' experimental-group means at every post-baseline time point, for a single
#' outcome. The contrast matrix is of the form:
#' | 0 -1 0 ... 0 | 0 1 0 ... 0 |
#' | 0 0 -1 ... 0 | 0 0 1 ... 0 |
#' | 0 0 0 ... -1 | 0 0 0 ... 1 |
#'
#' @param K integer, the number of post-baseline time points.
#'
#' @returns numeric matrix with `K` rows and `2 * (K + 1)` columns.
#' @export
build_omnibus_contrast_single_outcome <- function(K) {
  contrast <- matrix(0, nrow = K, ncol = 2 * (K + 1))
  contrast[, 2:(K + 1)] <- diag(-1, K)
  contrast[, (K + 3):(2 * (K + 1))] <- diag(1, K)
  contrast
}

# Contrast matrix for omnibus test with multiple outcomes.
#' Omnibus contrast matrix for multiple outcomes
#'
#' Block-diagonal stacking of [build_omnibus_contrast_single_outcome()] across
#' outcomes, testing equality of the control and experimental-group means at
#' every post-baseline time point for every outcome.
#'
#' @param J integer, the number of outcomes. Ignored (and recomputed) when
#'   `times` is supplied.
#' @param K integer, the number of post-baseline time points, shared across
#'   outcomes. Ignored (and recomputed per-outcome) when `times` is supplied.
#' @param times optional list of numeric time vectors, one per outcome, used
#'   when outcomes have differing numbers of time points.
#'
#' @returns numeric contrast matrix.
#' @export
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

#' Summing contrast matrix for multiple outcomes
#'
#' Builds a contrast matrix that sums the treatment effects across all
#' outcomes at each shared time point. All outcomes must share the same
#' number of measurements.
#'
#' @inheritParams build_omnibus_contrast_multi_outcome
#'
#' @returns numeric contrast matrix with `K` rows.
#' @export
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

#' Linear-slope contrast matrix for multiple outcomes
#'
#' Builds the plug-in contrast matrix for testing an outcome-specific linear
#' slope difference between treatment groups.
#'
#' @param J integer, the number of outcomes.
#' @param K integer, the number of post-baseline time points (used only when
#'   `times` is a plain numeric vector shared across outcomes).
#' @param times numeric vector (shared across outcomes) or list of numeric
#'   vectors, one per outcome.
#' @param Sigma numeric matrix, the first-stage covariance estimate.
#'
#' @returns numeric contrast matrix with `J` rows.
#' @export
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

#' Common linear-slope contrast matrix across outcomes
#'
#' Builds the plug-in contrast matrix for testing a treatment-effect slope
#' shared across all outcomes.
#'
#' @param times list of numeric time vectors, one per outcome.
#' @param Sigma numeric matrix, the first-stage covariance estimate.
#'
#' @returns numeric contrast matrix with 1 row.
#' @export
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
