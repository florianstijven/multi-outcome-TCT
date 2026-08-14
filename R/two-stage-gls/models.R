make_slowing_model <- function(times, ref, type) {
  if (type == "proportional") {
    time_mapping = function(t, gamma1) {
      gamma1 * t
    }
    null_gamma1 = 1
  } else if (type == "quadratic") {
    time_mapping = function(t, gamma1) {
      gamma1[1] * t + gamma1[2] * (t ** 2)
    }
    null_gamma1 = c(1, 0)
  } else {
    stop("Unknown slowing model type: ", type)
  }
  
  reference_trajectory_functions_list = reference_trajectory_f_list(times, ref)
  
  list(
    type        = "proportional",
    ref         = ref,
    times = times,
    reference_trajectory_functions_list = reference_trajectory_functions_list,
    null_gamma1 = null_gamma1,
    time_mapping = time_mapping,
    mu_from_gamma_null = function(gamma0) {
      control_means = reference_trajectory_functions_list$eval(times, gamma0)
      exp_means = reference_trajectory_functions_list$eval(time_mapping(times, null_gamma1), gamma0)
      c(control_means, exp_means)
    },
    mu_from_gamma = function(gamma) {
      # Number of reference-trajectory parameters.
      no_params_ref = reference_trajectory_functions_list$no_params
      
      gamma0 = gamma[1:no_params_ref]
      gamma1 = gamma[(no_params_ref + 1):length(gamma)]
      
      control_means = reference_trajectory_functions_list$eval(times, gamma0)
      exp_means = reference_trajectory_functions_list$eval(time_mapping(times, gamma1), gamma0)
      
      c(control_means, exp_means)
    }
  )
}

make_slowing_models <- function(times, ref, type) {
  lapply(times, make_slowing_model, ref = ref, type = type)
}

# Function factory to create the outcome-regression function for a list of
# slowing models. For a single outcome, this would be equivalent to the
# mu_from_gamma function in the slowing model.
mu_from_gamma_f_factory <- function(slowing_models) {
  function(gamma) {
    # Number of outcomes
    J <- length(slowing_models)
    
    # Number of reference-trajectory parameters for each outcome.
    no_params_ref_vec <- sapply(slowing_models, function(model)
      model$reference_trajectory_functions_list$no_params)
    
    # Number of slowing parameters for each outcome.
    no_params_slowing_vec <- sapply(slowing_models, function(model)
      length(model$null_gamma1))
    
    check_no_params(gamma, sum(no_params_ref_vec) + sum(no_params_slowing_vec))
    
    mu_blocks <- vector("list", J)
    for (outcome_idx in seq_len(J)) {
      # K <- no_of_measurements_vec[outcome_idx] - 1
      if (outcome_idx == 1) {
        start_col <- 1
      } else {
        start_col <- sum(no_params_ref_vec[1:(outcome_idx - 1)]) + sum(no_params_slowing_vec[1:(outcome_idx - 1)]) + 1
      }
      t1 <- start_col + no_params_ref_vec[outcome_idx] - 1
      gamma0 <- gamma[start_col:t1]
      gamma1 <- gamma[t1:(t1 + no_params_slowing_vec[outcome_idx] - 1)]
      
      mu_blocks[[outcome_idx]] <- slowing_models[[outcome_idx]]$mu_from_gamma(c(gamma0, gamma1))
    }
    do.call(c, mu_blocks)
  }
}

mu_from_gamma_null_f_factory <- function(slowing_models) {
  function(gamma) {
    # Number of outcomes
    J <- length(slowing_models)
    
    # Number of reference-trajectory parameters for each outcome.
    no_params_ref_vec <- sapply(slowing_models, function(model)
      model$reference_trajectory_functions_list$no_params)
    
    check_no_params(gamma, sum(no_params_ref_vec))
    
    mu_blocks <- vector("list", J)
    for (outcome_idx in seq_len(J)) {
      # K <- no_of_measurements_vec[outcome_idx] - 1
      if (outcome_idx == 1) {
        start_col <- 1
      } else {
        start_col <- sum(no_params_ref_vec[1:(outcome_idx - 1)]) + 1
      }
      t1 <- start_col + no_params_ref_vec[outcome_idx] - 1
      gamma0 <- gamma[start_col:t1]
      
      mu_blocks[[outcome_idx]] <- slowing_models[[outcome_idx]]$mu_from_gamma_null(gamma0)
    }
    do.call(c, mu_blocks)
  }
}

jac_mu_from_gamma_null_f_factory <- function(slowing_models) {
  function(gamma) {
    # Number of outcomes
    J <- length(slowing_models)

    # Number of reference-trajectory parameters for each outcome.
    no_params_ref_vec <- sapply(slowing_models, function(model)
      model$reference_trajectory_functions_list$no_params)
    
    check_no_params(gamma, sum(no_params_ref_vec))
    
    jac_blocks <- vector("list", J)
    for (outcome_idx in seq_len(J)) {
      # K <- no_of_measurements_vec[outcome_idx] - 1
      if (outcome_idx == 1) {
        start_col <- 1
      } else {
        start_col <- sum(no_params_ref_vec[1:(outcome_idx - 1)]) + 1
      }
      t1 <- start_col + no_params_ref_vec[outcome_idx] - 1
      gamma0 <- gamma[start_col:t1]
      
      times = rep(slowing_models[[outcome_idx]]$times, 2)
      jac_blocks[[outcome_idx]] <- 
        slowing_models[[outcome_idx]]$reference_trajectory_functions_list$jacobian(times = times, params = gamma0)[, ]
      
    }
    Matrix::bdiag(jac_blocks)
  }
}
