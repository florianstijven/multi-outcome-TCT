# Helper functions for local power.

#' Local power from a non-centrality parameter
#'
#' Computes local power from a non-centrality parameter, degrees of freedom,
#' and significance level, using the non-central chi-squared distribution.
#'
#' @param ncp numeric, the non-centrality parameter.
#' @param df integer, the degrees of freedom.
#' @param alpha numeric, the significance level.
#'
#' @returns numeric, the local power.
#' @export
power_from_ncp <- function(ncp, df, alpha) {
  critical_value <- qchisq(1 - alpha, df = df)
  1 - pchisq(critical_value, df = df, ncp = ncp)
}

#' Compute a local power curve
#'
#' Computes power for a range of effect sizes (`h_grid`) given the contrast
#' matrix `B`, significance level `alpha`, covariance matrix `Sigma`, and base
#' effect size `base_effect`.
#'
#' @param B numeric matrix, the plug-in contrast matrix.
#' @param h_grid numeric vector of effect-size multipliers.
#' @param alpha numeric, the significance level.
#' @param Sigma numeric matrix, the first-stage covariance estimate.
#' @param base_effect numeric vector, the base local-effect vector
#'   (see [local_shift_vector()]).
#'
#' @returns data.frame with columns `h`, `power`, `ncp`, and `df`.
#' @export
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

#' Local effect vector from a local parametric model's Jacobian
#'
#' @param jacobian numeric matrix, the Jacobian of the local parametric model.
#' @param h numeric, the local perturbation size.
#'
#' @returns numeric column vector, the local effect (base-effect) vector.
#' @export
local_shift_vector <- function(jacobian, h) {

  jacobian_gamma_theta <- matrix(0, nrow = ncol(jacobian), ncol = 1)
  jacobian_gamma_theta[1:ncol(jacobian), 1] <- h
  
  jacobian %*% jacobian_gamma_theta
}

#' Validate a per-outcome parameter list
#'
#' Checks that `params_list` is a list with exactly one element per outcome.
#'
#' @param params_list list of per-outcome parameter vectors.
#' @param J integer, the expected number of outcomes.
#'
#' @returns invisible `NULL`; called for its side effect of raising an error
#'   when validation fails.
#' @export
validate_params_list <- function(params_list, J) {
  if (!is.list(params_list)) {
    stop("params_list must be a list.")
  }
  if (length(params_list) != J) {
    stop("params_list must contain exactly one element per outcome (length J).")
  }
}



#' Local shift vector for a slowing model, one outcome
#'
#' @param ref character, the reference trajectory (`"4PL"` or `"nc_spline"`).
#' @param J integer, the number of outcomes.
#' @param K integer, the number of post-baseline time points.
#' @param times numeric vector (shared across outcomes) or list of numeric
#'   vectors, one per outcome.
#' @param h numeric, the local perturbation size.
#' @param params_list list of per-outcome reference-trajectory parameters.
#' @param ... additional arguments forwarded to the spline helpers.
#'
#' @returns numeric column vector, the local effect vector.
#' @export
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

#' Stacked mean vector for a reference trajectory model
#'
#' @param ref character, the reference trajectory (`"4PL"` or `"nc_spline"`).
#' @param J integer, the number of outcomes (inferred from `params_list`).
#' @param K integer, unused; retained for interface consistency.
#' @param times numeric vector (shared across outcomes) or list of numeric
#'   vectors, one per outcome.
#' @param params_list list of per-outcome reference-trajectory parameters.
#' @param ... additional arguments forwarded to the spline helpers.
#'
#' @returns numeric vector, the stacked mean values.
#' @export
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

#' Jacobian of a slowing model, one outcome (local-power helper)
#'
#' @param K integer, the number of post-baseline time points.
#' @param times numeric vector of measurement times.
#' @param ref character, the reference trajectory (`"4PL"` or `"nc_spline"`).
#' @param params numeric vector of reference-trajectory parameters.
#' @param ... additional arguments forwarded to the spline helpers.
#'
#' @returns numeric matrix, the stacked Jacobian for the control and
#'   experimental groups.
#' @export
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

#' Jacobian of a reference trajectory, dispatched by model type
#'
#' @inheritParams jacobian_slowing_single_outcome_temp
#'
#' @returns numeric matrix, the Jacobian of the reference trajectory.
#' @export
jacobian_ref_pm <- function(times, ref = "4PL", params, ...) {
  if (ref == "4PL") {
    jacobian_4PL(times = times, params = params)
  } else if (ref == "nc_spline") {
    jacobian_nc_spline(times = times, params = params, ...)
  } else {
    stop("Unknown reference model: ", ref)
  }
}

#' Time-derivative of a reference trajectory, dispatched by model type
#'
#' @inheritParams jacobian_slowing_single_outcome_temp
#'
#' @returns numeric vector, the time-derivative of the reference trajectory.
#' @export
ref_d <- function(times, ref = "4PL", params, ...) {
  if (ref == "4PL") {
    time_d_4PL(times = times, params = params)
  } else if (ref == "nc_spline") {
    time_d_nc_spline(times = times, params = params, ...)
  } else {
    stop("Unknown reference model: ", ref)
  }
}



#' Jacobian of a slowing model across multiple outcomes
#'
#' @param J integer, the number of outcomes.
#' @param K integer, unused; retained for interface consistency.
#' @param times numeric vector (shared across outcomes) or list of numeric
#'   vectors, one per outcome.
#' @param ref character, the reference trajectory (`"4PL"` or `"nc_spline"`).
#' @param params_list list of per-outcome reference-trajectory parameters.
#' @param slowing_only logical; if `TRUE`, retain only the columns
#'   corresponding to the treatment-effect (slowing) parameters.
#' @param ... additional arguments forwarded to the spline helpers.
#'
#' @returns numeric matrix, the stacked Jacobian across all outcomes.
#' @export
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


#' Jacobian of a slowing model with a shared treatment-effect parameter
#'
#' @inheritParams jacobian_slowing_multiple_outcomes
#'
#' @returns numeric matrix, the stacked Jacobian with treatment-effect
#'   parameters shared across outcomes.
#' @export
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
# Covariance structure
# ============================================================================

#' Kronecker-structured covariance matrix
#'
#' Builds a covariance matrix with a Kronecker-product structure across
#' treatment arms, outcomes, and time.
#'
#' @param K integer, the number of post-baseline time points.
#' @param J integer, the number of outcomes.
#' @param rho_time numeric, the AR(1)-type time correlation parameter.
#' @param rho_outcome numeric, the cross-outcome correlation parameter.
#'
#' @returns list with elements `Sigma` (the covariance matrix) and
#'   `time_points` (the time grid used to construct it).
#' @export
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

#' Check whether list elements have equal length
#'
#' @param lst a list.
#'
#' @returns logical scalar.
#' @export
check_equal_length <- function(lst) {
  lengths <- sapply(lst, length)
  all(lengths == lengths[1])
}

#' Check whether list elements are identical
#'
#' @param lst a list.
#'
#' @returns logical scalar.
#' @export
check_identical <- function(lst) {
  all(sapply(lst, function(x) identical(x, lst[[1]])))
}
