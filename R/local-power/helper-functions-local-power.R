# Helper functions for local power.


# Compute local power from non-centrality parameter (ncp), degrees of freedom (df), and significance level (alpha),
# using the non-central chi-squared distribution.
power_from_ncp <- function(ncp, df, alpha) {
  critical_value <- qchisq(1 - alpha, df = df)
  1 - pchisq(critical_value, df = df, ncp = ncp)
}

# Compute the power curve for a range of effect sizes (h_grid) given the contrast matrix B,
# significance level (alpha), covariance matrix (Sigma), and base effect size (base_effect).
compute_power_curve <- function(B, h_grid, alpha, Sigma, base_effect) {
  df = nrow(B)

  
  out <- data.frame(
    h = h_grid,
    power = NA_real_,
    ncp = NA_real_,
    df = df,
    stringsAsFactors = FALSE
  )
  BSigmaBinv <- solve(B %*% Sigma %*% t(B))
  for (idx in seq_along(h_grid)) {
    h <- h_grid[idx]
    delta <- h * base_effect
    ncp_value <- as.numeric(t(delta) %*% t(B) %*% BSigmaBinv %*% B %*% delta)
    out$power[idx] <- power_from_ncp(ncp_value, df = df, alpha = alpha)
    out$ncp[idx] <- ncp_value
  }
  
  out
}

# ============================================================================
# Local shift functions
# ============================================================================

# Compute the local effect vector (base_effect elsewhere) given the local
# parametric model's Jacobian.
local_shift_vector <- function(jacobian, h) {

  jacobian_gamma_theta <- matrix(0, nrow = ncol(jacobian), ncol = 1)
  jacobian_gamma_theta[1:ncol(jacobian), 1] <- h
  
  jacobian %*% jacobian_gamma_theta
}

validate_params_list <- function(params_list, J) {
  if (!is.list(params_list)) {
    stop("params_list must be a list.")
  }
  if (length(params_list) != J) {
    stop("params_list must contain exactly one element per outcome (length J).")
  }
}



local_shift_vector_slowing_outcome <- function(ref = "4PL",
                                               J,
                                               K,
                                               times,
                                               h = 1,
                                               params_list,
                                               ...) {
  if (!is.list(times) & !is.numeric(times)) {
    stop("times must be a numeric vector or a list of numeric vectors.")
  }
  if (is.numeric(times)) {
    times <- rep(list(times), J)
  }
  # number of outcomes
  J <- length(times)
  validate_params_list(params_list, J)
  
  # The local shift vector only depends on the Jacobian of the local parametric model.
  
  jacobian_local <- jacobian_slowing_multiple_outcomes(
    J = J,
    K = K,
    times,
    ref = ref,
    params_list = params_list,
    slowing_only = TRUE,
    ...
  )
  
  local_shift_vector(jacobian_local, h)
}

mean_vector <- function(ref = "4PL",
                        J,
                        K,
                        times,
                        params_list = list(),
                        ...) {
  if (!is.list(params_list)) {
    stop("params_list must be a list.")
  }
  J <- length(params_list)
  if (J == 0) {
    stop("params_list must not be empty.")
  }

  if (is.list(times)) {
    if (length(times) != J) {
      stop("When times is a list, it must contain one numeric vector per outcome.")
    }
    times_list <- times
  } else if (is.numeric(times)) {
    times_list <- rep(list(times), J)
  } else {
    stop("times must be a numeric vector or a list of numeric vectors.")
  }

  if (ref == "4PL") {
    mean_vector_dbl <- purrr::map2(times_list, params_list, function(outcome_times, outcome_params) {
      do.call(function_4PL, list(times = outcome_times, outcome_params))
    }) %>%
      do.call(c, .)
  } else if (ref == "nc_spline") {
    mean_vector_dbl <- purrr::map2(times_list, params_list, function(outcome_times, outcome_params) {
      do.call(function_nc_spline, list(times = outcome_times, params = outcome_params, ...))
    }) %>%
      do.call(c, .)
  } else {
    stop("Unknown reference model: ", ref)
  }
  
  mean_vector_dbl
}

# ============================================================================
# Jacobian functions
# ============================================================================

jacobian_slowing_single_outcome_temp <- function(K, times, ref = "4PL", params, ...) {
  jacobian_ref <- jacobian_ref_pm(times, ref = ref, params = params, ...)
  if (is.null(dim(jacobian_ref))) {
    jacobian_ref <- matrix(jacobian_ref, ncol = 1)
  }
  no_params_ref <- ncol(jacobian_ref)
  jacobian <- matrix(0, nrow = 2 * (K + 1), ncol = no_params_ref + 1)

  # Jacobian of the control group (first K+1 rows) with respect to reference parameters.
  jacobian[1:(K + 1), 1:no_params_ref] <- jacobian_ref
  # Jacobian of the treatment group (last K+1 rows) with respect to reference parameters.
  jacobian[(K + 2):(2 * (K + 1)), 1:no_params_ref] <- jacobian_ref
  # Jacobian of the treatment group with respect to the slowing parameter.
  jacobian[(K + 2):(2 * (K + 1)), no_params_ref + 1] <- ref_d(times, ref = ref, params, ...) * times

  jacobian
}

jacobian_ref_pm <- function(times, ref = "4PL", params, ...) {
  if (ref == "4PL") {
    jacobian_4PL(times = times, params = params)
  } else if (ref == "nc_spline") {
    jacobian_nc_spline(times = times, params = params, ...)
  } else {
    stop("Unknown reference model: ", ref)
  }
}

ref_d <- function(times, ref = "4PL", params, ...) {
  if (ref == "4PL") {
    time_d_4PL(times = times, params = params)
  } else if (ref == "nc_spline") {
    time_d_nc_spline(times = times, params = params, ...)
  } else {
    stop("Unknown reference model: ", ref)
  }
}



jacobian_slowing_multiple_outcomes <- function(J,
                                               K,
                                               times,
                                               ref = "4PL",
                                               params_list = list(),
                                               slowing_only = FALSE,
                                               ...
                                               ) {
  if (!is.list(times) & !is.numeric(times)) {
    stop("times must be a numeric vector or a list of numeric vectors.")
  }
  if (is.numeric(times)) {
    times <- rep(list(times), J)
  }
  # number of outcomes
  J <- length(times)
  validate_params_list(params_list, J)

  # number of measurements
  no_of_measurements <- purrr::map_dbl(times, ~ length(.x)) %>%
    sum()

  jacobian_blocks <- purrr::map2(times, params_list, function(outcome_times, outcome_params) {
    do.call(jacobian_slowing_single_outcome_temp,
              list(
                K = length(outcome_times) - 1,
                times = outcome_times,
                ref = ref,
                params = outcome_params, 
                ...
              ))
  })
  no_params_total_vec <- purrr::map_int(jacobian_blocks, ncol)
  total_no_params <- sum(no_params_total_vec)
  
  jacobian <- matrix(0, nrow = no_of_measurements * 2, ncol = total_no_params)
  for (outcome_idx in seq_len(J)) {
    K <- length(times[[outcome_idx]]) - 1
    start_row <- if (outcome_idx == 1) {
      1
    } else {
      sum(purrr::map_dbl(times[1:(outcome_idx - 1)], ~ 2 * length(.x))) + 1
    }
    end_row <- start_row + 2 * (K + 1) - 1
    start_col <- if (outcome_idx == 1) {
      1
    } else {
      sum(no_params_total_vec[1:(outcome_idx - 1)]) + 1
    }
    end_col <- start_col + no_params_total_vec[outcome_idx] - 1

    jacobian_block <- jacobian_blocks[[outcome_idx]]

    if (ncol(jacobian_block) != no_params_total_vec[outcome_idx]) {
      stop("Jacobian block has unexpected number of columns for outcome ", outcome_idx, ".")
    }

    jacobian[start_row:end_row, start_col:end_col] <- jacobian_block
    
  }
  # Remove all columns corresponding to the reference trajectory parameters if slowing_only is TRUE. This leaves only the
  # treatment effect parameters in the Jacobian.
  if (slowing_only) {
    slowing_cols <- cumsum(no_params_total_vec)
    jacobian <- jacobian[, slowing_cols, drop = FALSE]
  }
  
  jacobian
}


jacobian_slowing_multiple_outcomes_shared <- function(J,
                                                      K,
                                                      times,
                                                      ref = "4PL",
                                                      params_list = list(),
                                                      slowing_only = FALSE,
                                                      ...) {
  if (!is.list(times) & !is.numeric(times)) {
    stop("times must be a numeric vector or a list of numeric vectors.")
  }
  if (is.numeric(times)) {
    times <- rep(list(times), J)
  }
  # number of outcomes
  J <- length(times)
  validate_params_list(params_list, J)
  
  jacobian_outcome_specific <- jacobian_slowing_multiple_outcomes(
    J = J,
    K = K,
    times = times,
    ref = ref,
    params_list = params_list,
    slowing_only = FALSE,
    ...
  )

  jacobian_blocks <- purrr::map2(times, params_list, function(outcome_times, outcome_params) {
    do.call(jacobian_slowing_single_outcome_temp, list(
      K = length(outcome_times) - 1,
      times = outcome_times,
      ref = ref,
      params = outcome_params, 
      ...
    ))
  })
  no_params_total_vec <- purrr::map_int(jacobian_blocks, ncol)
  no_params_ref_vec <- no_params_total_vec - 1
  total_no_params <- sum(no_params_total_vec)
  total_no_params_ref <- sum(no_params_ref_vec)
  
  jacobian_gamma_theta <- matrix(0,
                                 nrow = total_no_params,
                                 ncol = total_no_params_ref + 1)
  
  row_offset <- 0
  col_offset <- 0
  for (outcome_idx in seq_len(J)) {
    no_params_ref <- no_params_ref_vec[outcome_idx]
    no_params_total <- no_params_total_vec[outcome_idx]

    ref_row_idx <- (row_offset + 1):(row_offset + no_params_ref)
    shared_row_idx <- row_offset + no_params_total

    if (no_params_ref > 0) {
      ref_col_idx <- (col_offset + 2):(col_offset + no_params_ref + 1)
      jacobian_gamma_theta[ref_row_idx, ref_col_idx] <- diag(1, no_params_ref)
    }
    jacobian_gamma_theta[shared_row_idx, 1] <- 1

    row_offset <- row_offset + no_params_total
    col_offset <- col_offset + no_params_ref
  }
  
  jacobian_shared <- jacobian_outcome_specific %*% jacobian_gamma_theta
  # Only retain the shared treatment effect parameter if slowing_only is TRUE.
  if (slowing_only) {
    jacobian_shared <- jacobian_shared[, 1, drop = FALSE]
  }
  jacobian_shared
}

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

# ============================================================================
# 4PL Reference Trajectory Functions
# ============================================================================


# Derivative of the 4PL curve with respect to time, given the slope and inflection point parameters.
time_d_4PL <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  
  exp_term <- exp(-slope * (times - inflection))
  denom <- (1 + exp_term)^2
  slope * exp_term / denom
}

function_4PL <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  
  # Compute the 4PL model values given times, slope, and inflection point.
  1 / (1 + exp(-slope * (times - inflection)))
}

jacobian_4PL <- function(times, params) {
  inflection <- params[1]
  slope <- params[2]
  
  # Compute the Jacobian of the 4PL model with respect to the slope and inflection point parameters.
  exp_term <- exp(-slope * (times - inflection))
  denom <- (1 + exp_term)^2
  d_inflection <- -slope * exp_term / denom
  d_slope <- (times - inflection) * exp_term / denom
  cbind(d_inflection, d_slope)
}

# ============================================================================
# Natural Cubic Spline Interpolation Reference Trajectory Functions
# ============================================================================

build_nc_spline_basis <- function(times, knots, boundary_knots, params = NULL, ...) {
  if (!is.numeric(times)) {
    stop("times must be a numeric vector.")
  }
  if (length(times) == 0) {
    stop("times must not be empty.")
  }

  if (!is.numeric(boundary_knots) || length(boundary_knots) != 2) {
    stop("boundary knots must be a numeric vector of length 2.")
  }

  basis <- splines::ns(
    x = times,
    knots = knots,
    Boundary.knots = boundary_knots,
    intercept = TRUE
  )

  if (!is.null(params) && length(params) != ncol(basis)) {
    stop(
      "Length of params (", length(params),
      ") does not match spline basis dimension (", ncol(basis), ")."
    )
  }

  basis
}

time_d_nc_spline <- function(times, knots, boundary_knots, params) {
  if (!is.numeric(params)) {
    stop("params must be a numeric vector.")
  }
  if (!is.numeric(times)) {
    stop("times must be a numeric vector.")
  }

  # Central finite difference for d/dt f(t), with adaptive step size.
  eps <- sqrt(.Machine$double.eps) * pmax(1, abs(times))
  f_plus <- function_nc_spline(times + eps, knots = knots, boundary_knots = boundary_knots, params = params)
  f_minus <- function_nc_spline(times - eps, knots = knots, boundary_knots = boundary_knots, params = params)

  as.numeric((f_plus - f_minus) / (2 * eps))
}

function_nc_spline <- function(times, params, knots, boundary_knots) {
  if (!is.numeric(params)) {
    stop("params must be a numeric vector.")
  }
  basis <- build_nc_spline_basis(
    times = times,
    knots = knots,
    boundary_knots = boundary_knots,
    params = params
  )
  as.numeric(basis %*% params)
}

jacobian_nc_spline <- function(times, knots, boundary_knots, params) {
  if (!is.numeric(params)) {
    stop("params must be a numeric vector.")
  }
  # The model is linear in coefficients, so the Jacobian is the spline basis.
  build_nc_spline_basis(
    times = times,
    knots = knots,
    boundary_knots = boundary_knots,
    params = params
  )
}

# ============================================================================
# Covariance structure
# ============================================================================

build_kronecker_covariance <- function(K,
                                       J,
                                       rho_time = 0.8,
                                       rho_outcome = 0.75) {
  time_points <- seq(0, 3, length.out = K + 1)
  time_cov <- matrix(0, nrow = K + 1, ncol = K + 1)
  for (i in 1:(K + 1)) {
    for (j in 1:(K + 1)) {
      time_cov[i, j] <- sqrt(1 + time_points[i]) * sqrt(1 + time_points[j]) *
        rho_time^abs(time_points[i] - time_points[j])
    }
  }
  
  outcome_cov <- matrix(rho_outcome, nrow = J, ncol = J)
  diag(outcome_cov) <- 1
  
  cov <- diag(1, 2) %x% outcome_cov %x% time_cov
  
  list(Sigma = cov, time_points = time_points)
}


# ============================================================================
# Varia
# ============================================================================

# Check whether the elements of a list are of equal length.
check_equal_length <- function(lst) {
  lengths <- sapply(lst, length)
  all(lengths == lengths[1])
}

# Check whether the elements of a list are identical.
check_identical <- function(lst) {
  all(sapply(lst, function(x) identical(x, lst[[1]])))
}


reference_trajectory_f_list <- function(times, ref) {
  if (ref == "4PL") {
    no_params = 2
    list(
      eval = function(times, params) {
        check_no_params(params, no_params)
        function_4PL(times, params)
      }
        ,
      jacobian = function(times, params) {
        check_no_params(params, no_params)
        jacobian_4PL(times, params)
      },
      time_derivative = function(times, params) {
        check_no_params(params, no_params)
        time_d_4PL(times, params)
      },
      no_params = no_params
      
    )
  } else if (ref == "nc_spline") {
    knots = times[c(-1, -length(times))]
    boundary_knots = range(times)
    
    no_params = length(knots) + 2
    
    list(
      eval = function(times, params) {
        check_no_params(params, no_params)
        function_nc_spline(times = times, params = params, knots = knots, boundary_knots = boundary_knots)
      },
      jacobian = function(times, params) {
        check_no_params(params, no_params)
        jacobian_nc_spline(times = times, params = params, knots = knots, boundary_knots = boundary_knots)
      },
      time_derivative = function(times, params) {
        check_no_params(params, no_params)
        time_d_nc_spline(times = times, params = params, knots = knots, boundary_knots = boundary_knots)
      },
      no_params = no_params
    )
  }
}

check_no_params <- function(params, no_params) {
  if (length(params) != no_params) {
    stop("Length of params (", length(params), ") does not match expected number of parameters (", no_params, ").")
  }
}








