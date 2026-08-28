#' Slowing-model construction and mean-function helpers
#'
#' This file builds the working-model objects used to represent outcome-specific
#' slowing patterns and to assemble the stacked mean vectors for the null and
#' full GLS models. The functions here encode the shared reference trajectory,
#' the outcome-specific time mapping, and the concatenated mean structures used
#' later in the second-stage fitting routines.
#'
#' @keywords internal
#' @name slowing-models
NULL

#' Construct a slowing model for one outcome
#'
#' Creates a model object for a single outcome. The model includes the reference
#' trajectory shared by the control and experimental groups, along with the time
#' mapping that relates measurements under treatment to the underlying reference
#' schedule through a proportional or quadratic slowing function.
#'
#' @param times numeric vector of measurement times for the outcome.
#' @param ref reference trajectory specification passed to
#'   `reference_trajectory_f_list()`.
#' @param type character scalar. One of `"proportional"` or `"quadratic"`.
#'
#' @returns list containing the model metadata, time-mapping function, null
#'   slowing parameters, and functions for evaluating the mean vector under the
#'   null and full parameterizations.
#' @export
make_slowing_model <- function(times, ref, type) {
  if (type == "proportional") {
    time_mapping = function(t, gamma1) {
      gamma1 * t
    }
    null_gamma1 = 1
    
    time_mapping_deriv = function(t, gamma1) {
      t
    }
  } else if (type == "quadratic") {
    time_mapping = function(t, gamma1) {
      gamma1[1] * t + gamma1[2] * (t ** 2)
    }
    null_gamma1 = c(1, 0)
    
    time_mapping_deriv = function(t, gamma1) {
      cbind(t, t ** 2)
    }
  } else {
    stop("Unknown slowing model type: ", type)
  }
  
  reference_trajectory_functions_list = reference_trajectory_f_list(times, ref)
  
  # Names for model parameters
  param_names <- c(
    paste0("ref_param_", seq_len(reference_trajectory_functions_list$no_params)),
    paste0("slowing_param_", seq_len(length(null_gamma1)))
  )
  
  model(
    mean_fn = function(gamma) {
      # Number of reference-trajectory parameters.
      no_params_ref = reference_trajectory_functions_list$no_params
      
      gamma0 = gamma[1:no_params_ref]
      gamma1 = gamma[(no_params_ref + 1):length(gamma)]
      
      control_means = reference_trajectory_functions_list$eval(times, gamma0)
      exp_means = reference_trajectory_functions_list$eval(time_mapping(times, gamma1), gamma0)
      
      c(control_means, exp_means)
    },
    jacobian_fn = function(gamma) {
      # Number of reference-trajectory parameters.
      no_params_ref = reference_trajectory_functions_list$no_params
      
      gamma0 = gamma[1:no_params_ref]
      gamma1 = gamma[(no_params_ref + 1):length(gamma)]
      
      jacobian_slowing_single_outcome(
        times = times,
        ref_d = reference_trajectory_functions_list$time_derivative,
        jacobian_ref = reference_trajectory_functions_list$jacobian,
        gamma0 = gamma0,
        gamma1 = gamma1,
        time_mapping = time_mapping,
        time_mapping_deriv = time_mapping_deriv
      )
    },
    nuisance_params_position = seq_len(reference_trajectory_functions_list$no_params),
    treatment_params_null = null_gamma1,
    param_names = param_names
  )
}

#' Construct one slowing model per outcome
#'
#' Applies `make_slowing_model()` to each time vector in a list, generating a
#' model object for each outcome under the same reference-trajectory structure
#' and slowing-model type.
#'
#' @param times list of numeric measurement-time vectors, one per outcome.
#' @param ref reference trajectory specification passed to each call to
#'   `make_slowing_model()`.
#' @param type character scalar passed to `make_slowing_model()`.
#'
#' @returns list of slowing-model objects, one for each outcome.
#' @export
make_slowing_models <- function(times, ref, type, outcome_names = NULL) {
  if (is.null(outcome_names)) {
    outcome_names <- paste0("Outcome ", seq_along(times))
  }
  concatenate_models(lapply(times, function(t) make_slowing_model(t, ref, type)), submodel_names = outcome_names)
}

#' Jacobian of a slowing model, one outcome
#'
#' @param times numeric vector of measurement times.
#' @param ref_d function, the time-derivative of the reference trajectory.
#' @param jacobian_ref function, the Jacobian of the reference trajectory.
#' @param gamma0 numeric vector, the reference-trajectory parameters.
#' @param gamma1 numeric vector, the time-mapping (slowing) parameters.
#' @param time_mapping function mapping `(t, gamma1)` to mapped time points.
#' @param time_mapping_deriv function, the derivative of `time_mapping` with
#'   respect to `gamma1`.
#'
#' @returns numeric matrix, the stacked Jacobian for the control and
#'   experimental groups.
#' @export
jacobian_slowing_single_outcome <- function(times, ref_d, jacobian_ref, gamma0, gamma1, time_mapping, time_mapping_deriv){

  control_jacobian = cbind(jacobian_ref(times, gamma0),
                           matrix(0, nrow = length(times), ncol = length(gamma1)))
  exp_jacobian = cbind(
    jacobian_ref(time_mapping(times, gamma1), gamma0),
    diag(ref_d(time_mapping(times, gamma1), gamma0)) %*% time_mapping_deriv(times, gamma1)
  )
  
  rbind(control_jacobian, exp_jacobian)
}

#' Check that a parameter vector has the expected length
#'
#' @param params numeric vector of parameters.
#' @param no_params integer, the expected length.
#'
#' @returns invisible `NULL`; called for its side effect of raising an error
#'   when the lengths don't match.
#' @export
check_no_params <- function(params, no_params) {
  if (length(params) != no_params) {
    stop("Length of params (", length(params), ") does not match expected number of parameters (", no_params, ").")
  }
}


#' Reference-trajectory function list, dispatched by model type
#'
#' Builds a list of functions (`eval`, `jacobian`, `time_derivative`) and the
#' expected parameter count (`no_params`) for the requested reference
#' trajectory.
#'
#' @param times numeric vector of measurement times, used to derive knots for
#'   the `"nc_spline"` reference.
#' @param ref character, one of `"4PL"`, `"nc_spline"`, or `"4PL_asympt"`.
#'
#' @returns list with elements `eval`, `jacobian`, `time_derivative` (each a
#'   function of `times` and `params`), and `no_params` (integer).
#' @export
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
        function_nc_spline(
          times = times,
          params = params,
          knots = knots,
          boundary_knots = boundary_knots
        )
      },
      jacobian = function(times, params) {
        check_no_params(params, no_params)
        jacobian_nc_spline(
          times = times,
          params = params,
          knots = knots,
          boundary_knots = boundary_knots
        )
      },
      time_derivative = function(times, params) {
        check_no_params(params, no_params)
        time_d_nc_spline(
          times = times,
          params = params,
          knots = knots,
          boundary_knots = boundary_knots
        )
      },
      no_params = no_params
    )
  } else if (ref == "4PL_asympt") {
    no_params = 4
    
    list(
      eval = function(times, params) {
        check_no_params(params, no_params)
        function_4PL_asympt(times, params)
      }
      ,
      jacobian = function(times, params) {
        check_no_params(params, no_params)
        jacobian_4PL_asympt(times, params)
      },
      time_derivative = function(times, params) {
        check_no_params(params, no_params)
        time_d_4PL_asympt(times, params)
      },
      no_params = no_params
    )
  } else {
    stop("Unknown reference model: ", ref)
  }
}
