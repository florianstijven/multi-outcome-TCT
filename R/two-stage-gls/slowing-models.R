#' Slowing-model construction and mean-function helpers
#'
#' This file builds the working-model objects used to represent outcome-specific
#' slowing patterns and to assemble the stacked mean vectors for the null and
#' full GLS models. The functions here encode the shared reference trajectory,
#' the outcome-specific time mapping, and the concatenated mean structures used
#' later in the second-stage fitting routines.
#'
#' @keywords internal

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
    treatment_params_null = null_gamma1
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
make_slowing_models <- function(times, ref, type) {
  concatenate_models(lapply(times, function(t) make_slowing_model(t, ref, type)))
}

jacobian_slowing_single_outcome <- function(times, ref_d, jacobian_ref, gamma0, gamma1, time_mapping, time_mapping_deriv){
  
  
  
  
  control_jacobian = cbind(jacobian_ref(times, gamma0),
                           matrix(0, nrow = length(times), ncol = length(gamma1)))
  exp_jacobian = cbind(
    jacobian_ref(time_mapping(times, gamma1), gamma0),
    diag(ref_d(time_mapping(times, gamma1), gamma0)) %*% time_mapping_deriv(times, gamma1)
  )
  
  rbind(control_jacobian, exp_jacobian)
}

#' #' Create a full-model mean function from a list of slowing models
#' #'
#' #' Builds a function that maps a concatenated parameter vector to the vector of
#' #' control and experimental mean responses across all outcomes. The parameter
#' #' vector is arranged outcome by outcome, with reference-trajectory parameters
#' #' followed by slowing parameters.
#' #'
#' #' @param slowing_models list returned by `make_slowing_models()`.
#' #'
#' #' @returns function taking the full parameter vector and returning the stacked
#' #'   mean vector for all outcomes.
#' mu_from_gamma_f_factory <- function(slowing_models) {
#'   function(gamma) {
#'     # Number of outcomes
#'     J <- length(slowing_models)
#'     
#'     # Number of reference-trajectory parameters for each outcome.
#'     no_params_ref_vec <- sapply(slowing_models, function(model)
#'       model$reference_trajectory_functions_list$no_params)
#'     
#'     # Number of slowing parameters for each outcome.
#'     no_params_slowing_vec <- sapply(slowing_models, function(model)
#'       length(model$null_gamma1))
#'     
#'     check_no_params(gamma, sum(no_params_ref_vec) + sum(no_params_slowing_vec))
#'     
#'     mu_blocks <- vector("list", J)
#'     for (outcome_idx in seq_len(J)) {
#'       # K <- no_of_measurements_vec[outcome_idx] - 1
#'       if (outcome_idx == 1) {
#'         start_col <- 1
#'       } else {
#'         start_col <- sum(no_params_ref_vec[1:(outcome_idx - 1)]) + sum(no_params_slowing_vec[1:(outcome_idx - 1)]) + 1
#'       }
#'       t1 <- start_col + no_params_ref_vec[outcome_idx] - 1
#'       gamma0 <- gamma[start_col:t1]
#'       gamma1 <- gamma[t1:(t1 + no_params_slowing_vec[outcome_idx] - 1)]
#'       
#'       mu_blocks[[outcome_idx]] <- slowing_models[[outcome_idx]]$mu_from_gamma(c(gamma0, gamma1))
#'     }
#'     do.call(c, mu_blocks)
#'   }
#' }
#' 
#' #' Create the null-model mean function for a list of slowing models
#' #'
#' #' Under the null model, only the reference-trajectory parameters are estimated,
#' #' while each outcome uses its own model-specific null slowing values. This
#' #' factory returns a function that maps the concatenated reference parameters to
#' #' the stacked null mean vector.
#' #'
#' #' @param slowing_models list returned by `make_slowing_models()`.
#' #'
#' #' @returns function taking the concatenated reference-trajectory parameters and
#' #'   returning the corresponding null mean vector.
#' mu_from_gamma_null_f_factory <- function(slowing_models) {
#'   function(gamma) {
#'     # Number of outcomes
#'     J <- length(slowing_models)
#'     
#'     # Number of reference-trajectory parameters for each outcome.
#'     no_params_ref_vec <- sapply(slowing_models, function(model)
#'       model$reference_trajectory_functions_list$no_params)
#'     
#'     check_no_params(gamma, sum(no_params_ref_vec))
#'     
#'     mu_blocks <- vector("list", J)
#'     for (outcome_idx in seq_len(J)) {
#'       # K <- no_of_measurements_vec[outcome_idx] - 1
#'       if (outcome_idx == 1) {
#'         start_col <- 1
#'       } else {
#'         start_col <- sum(no_params_ref_vec[1:(outcome_idx - 1)]) + 1
#'       }
#'       t1 <- start_col + no_params_ref_vec[outcome_idx] - 1
#'       gamma0 <- gamma[start_col:t1]
#'       
#'       mu_blocks[[outcome_idx]] <- slowing_models[[outcome_idx]]$mu_from_gamma_null(gamma0)
#'     }
#'     do.call(c, mu_blocks)
#'   }
#' }
#' 
#' #' Create the Jacobian of the null-model mean function
#' #'
#' #' Computes the block-diagonal Jacobian for the null mean function. Each block
#' #' corresponds to one outcome and contains the derivatives of the stacked control
#' #' and experimental means with respect to that outcome's reference-trajectory
#' #' parameters.
#' #'
#' #' @param slowing_models list returned by `make_slowing_models()`.
#' #'
#' #' @returns function taking the null parameter vector and returning the block-
#' #'   diagonal Jacobian matrix.
#' jac_mu_from_gamma_null_f_factory <- function(slowing_models) {
#'   function(gamma) {
#'     # Number of outcomes
#'     J <- length(slowing_models)
#'     
#'     # Number of reference-trajectory parameters for each outcome.
#'     no_params_ref_vec <- sapply(slowing_models, function(model)
#'       model$reference_trajectory_functions_list$no_params)
#'     
#'     check_no_params(gamma, sum(no_params_ref_vec))
#'     
#'     jac_blocks <- vector("list", J)
#'     for (outcome_idx in seq_len(J)) {
#'       # K <- no_of_measurements_vec[outcome_idx] - 1
#'       if (outcome_idx == 1) {
#'         start_col <- 1
#'       } else {
#'         start_col <- sum(no_params_ref_vec[1:(outcome_idx - 1)]) + 1
#'       }
#'       t1 <- start_col + no_params_ref_vec[outcome_idx] - 1
#'       gamma0 <- gamma[start_col:t1]
#'       
#'       times = rep(slowing_models[[outcome_idx]]$times, 2)
#'       jac_blocks[[outcome_idx]] <-
#'         slowing_models[[outcome_idx]]$reference_trajectory_functions_list$jacobian(times = times, params = gamma0)[, ]
#'       
#'     }
#'     Matrix::bdiag(jac_blocks)
#'   }
#' }