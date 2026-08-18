new_model <- function(x = list()) {
  stopifnot(is.list(x))
  
  # Add null functions for mean and Jacobian to the model object.
  mean_fn_null <- mean_fn_null_constructor(
    mean_fn = x$mean_fn,
    nuisance_params_position = x$nuisance_params_position,
    treatment_params_null = x$treatment_params_null
  )
  
  jacobian_fn_null <- jacobian_fn_null_constructor(
    jacobian_fn = x$jacobian_fn,
    nuisance_params_position = x$nuisance_params_position,
    treatment_params_null = x$treatment_params_null
  )
  
  # Jacobian of the mean function with respect to treatment-effect parameters at
  # the null, given the (estimated) nuisance parameters.
  jacobian_fn_treatment_null <- jacobian_fn_treatment_null_constructor(
    jacobian_fn = x$jacobian_fn,
    nuisance_params_position = x$nuisance_params_position,
    treatment_params_null = x$treatment_params_null
  )
  
  x$mean_fn_null <- mean_fn_null
  x$jacobian_fn_null <- jacobian_fn_null
  x$jacobian_fn_treatment_null <- jacobian_fn_treatment_null
  
  
  structure(
    x,
    class = "model"
  )
}

mean_fn_null_constructor <- function(mean_fn, nuisance_params_position, treatment_params_null) {
  function(gamma0) {
    gamma <- numeric(length = length(nuisance_params_position) + length(treatment_params_null))
    gamma[nuisance_params_position] <- gamma0
    gamma[-nuisance_params_position] <- treatment_params_null
    
    mean_fn(gamma)
  }
}

jacobian_fn_null_constructor <- function(jacobian_fn, nuisance_params_position, treatment_params_null) {
  function(gamma0) {
    gamma <- numeric(length = length(nuisance_params_position) + length(treatment_params_null))
    gamma[nuisance_params_position] <- gamma0
    gamma[-nuisance_params_position] <- treatment_params_null

    jacobian <- jacobian_fn(gamma)

    jacobian[, nuisance_params_position]
  }
}

jacobian_fn_treatment_null_constructor <- function(jacobian_fn, nuisance_params_position, treatment_params_null) {
  function(gamma0) {
    gamma <- numeric(length = length(nuisance_params_position) + length(treatment_params_null))
    gamma[nuisance_params_position] <- gamma0
    gamma[-nuisance_params_position] <- treatment_params_null
    
    jacobian <- jacobian_fn(gamma)
    
    jacobian[, -nuisance_params_position, drop = FALSE]
  }
}

model <- function(mean_fn, jacobian_fn, nuisance_params_position, treatment_params_null) {
  new_model(
    list(
      mean_fn = mean_fn,
      jacobian_fn = jacobian_fn,
      nuisance_params_position = nuisance_params_position,
      treatment_params_null = treatment_params_null
    )
  )
}

concatenate_models <- function(models) {
  if (!all(sapply(models, inherits, "model"))) {
    stop("All elements of 'models' must be of class 'model'.")
  } else {
    lapply(models, validate_model)
  }
  
  # Vector with the number of nuisance parameters for each model.
  nuisance_params_counts <- sapply(models, function(model) length(model$nuisance_params_position))
  # Vector with the number of treatment parameters for each model.
  treatment_params_counts <- sapply(models, function(model) length(model$treatment_params_null))
  # Vector with the total number of parameters for each model.
  total_params_counts <- nuisance_params_counts + treatment_params_counts
  # Vector with the starting index of each model's parameters in the combined parameter vector.
  start_indices_full <- cumsum(c(1, total_params_counts[-length(total_params_counts)]))
  start_indices_null <- cumsum(c(1, nuisance_params_counts[-length(nuisance_params_counts)]))
  
  mean_fn_combined <- function(gamma) {
    do.call(c,
            purrr::map2(
              .x = models,
              .y = start_indices_full,
              .f =  function(model, start_index) {
                end_index <- start_index + length(model$nuisance_params_position) + length(model$treatment_params_null) - 1
                model$mean_fn(gamma[start_index:end_index])
              }
            ))
  }
  
  jacobian_fn_combined <- function(gamma) {
    do.call(Matrix::bdiag, purrr::map2(.x = models, .y = start_indices_full, .f = function(model, start_index) {
      end_index <- start_index + length(model$nuisance_params_position) + length(model$treatment_params_null) - 1
      model$jacobian_fn(gamma[start_index:end_index])
    }))
  }
  
  nuisance_params_position_combined <- unlist(purrr::map2(.x = models, .y = start_indices_full, function(model, start_index) model$nuisance_params_position + start_index - 1))
  treatment_params_null_combined <- unlist(lapply(models, function(model) model$treatment_params_null))
  
  new_model(
    list(
      mean_fn = mean_fn_combined,
      jacobian_fn = jacobian_fn_combined,
      nuisance_params_position = nuisance_params_position_combined,
      treatment_params_null = treatment_params_null_combined
    )
  )
}

shared_parameter_model <- function(model, shared_param_positions) {
  # `shared_param_positions` must be a numeric vector of positive integers
  # indicating the positions of the shared parameters in the full parameter
  # vector, or it should be a list of such vectors, one for each shared
  # parameter.
  if (!(is.list(shared_param_positions) || is.numeric(shared_param_positions))) {
    stop("'shared_param_positions' must be a numeric vector of positive integers, or a list of such vectors.")
  }
  if (!is.list(shared_param_positions)) {
    # We further work with a list of shared parameter positions, even if there
    # is only one shared parameter.
    shared_param_positions <- list(shared_param_positions)
  }
  
  # Construct a model for the parameters of `model`, taking into account the
  # shared parameters. The new parameter vector will have one entry for each
  # shared parameter, and one entry for each non-shared parameter. The mean
  # function of the new model will map the new parameter vector to the full
  # parameter vector of `model`.
  no_of_non_shared_params <- length(model$nuisance_params_position) + length(model$treatment_params_null) - length(unlist(shared_param_positions))
  no_of_params <- length(model$nuisance_params_position) + length(model$treatment_params_null)
  no_of_non_shared_treatment_params <- length(model$treatment_params_null) - length(unlist(shared_param_positions))
  
  mean_fn_params <- function(gamma_shared) {
    # Construct the full parameter vector for `model` from the parameters of
    # `gamma_shared`, taking into account the shared parameters.
    gamma_full <- numeric(length = no_of_params)
    gamma_full[-unlist(shared_param_positions)] <- gamma_shared[seq_along(no_of_non_shared_params)]
    
    # Fill in the shared parameters.
    for (i in seq_along(shared_param_positions)) {
      shared_param_positions_i <- shared_param_positions[[i]]
      gamma_full[shared_param_positions_i] <- gamma_shared[no_of_non_shared_params + i]
    }
    gamma_full
  }
  
  jacobian_fn_params <- function(gamma_shared) {
    jacobian_full <- matrix(0, nrow = no_of_params, ncol = length(gamma_shared))
    # Fill in the Jacobian for the non-shared parameters.
    jacobian_full[-unlist(shared_param_positions), seq_len(no_of_non_shared_params)] <- diag(1, nrow = no_of_non_shared_params)
    
    # Fill in the Jacobian for the shared parameters.
    for (i in seq_along(shared_param_positions)) {
      shared_param_positions_i <- shared_param_positions[[i]]
      jacobian_full[shared_param_positions_i, no_of_non_shared_params + i] <- 1
    }
    jacobian_full
  }
  
  # The function assumes that only treatment-effect parameters can be shared.
  # The treatment-effect parameters are grouped together at the end of the
  # parameter vector, so the nuisance parameters the first
  # `length(model$nuisance_params_position)` entries of the parameter vector.
  nuisance_params_position = seq_len(length(model$nuisance_params_position))
  
  # Determine the null values for the non-shared treatment-effect parameters.
  non_shared_treatment_params_position_complement <- c(model$nuisance_params_position, unlist(shared_param_positions))
  non_shared_treatment_params_position <- seq_len(no_of_params)[-non_shared_treatment_params_position_complement]
  treatment_params_position <-  seq_len(no_of_params)[-model$nuisance_params_position]
  # Relative position of the non-shared treatment-effect parameters in the
  # original treatment-effect parameter vector.
  non_shared_treatment_params_position_subset <- sapply(
    non_shared_treatment_params_position,
    function(pos) {
      which.max(treatment_params_position == pos)
    }
  )
  shared_treatment_params_position_subset <- 
    sapply(
      shared_param_positions,
      function(positions) {
        which.max(treatment_params_position %in% positions)
      }
    )
  
  if (length(non_shared_treatment_params_position_subset) == 0) {
    non_shared_treatment_params_position_subset <- integer(0)
  }

  treatment_params_null <- numeric(length = no_of_non_shared_treatment_params + length(shared_param_positions))
  
  treatment_params_null[seq_len(no_of_non_shared_treatment_params)] <- model$treatment_params_null[non_shared_treatment_params_position_subset]
  treatment_params_null[1:length(shared_param_positions) + no_of_non_shared_treatment_params] <- model$treatment_params_null[shared_treatment_params_position_subset]
  
  model_for_params <- model(
    mean_fn = mean_fn_params,
    jacobian_fn = jacobian_fn_params,
    nuisance_params_position = nuisance_params_position,
    treatment_params_null = treatment_params_null
  )
  model_for_model_params(model, model_for_params)
}

model_for_model_params <- function(model, model_for_params) {
  mean_fn_model_params <- function(gamma_model_params) {
    gamma_full <- model_for_params$mean_fn(gamma_model_params)
    model$mean_fn(gamma_full)
  }
  
  jacobian_fn_model_params <- function(gamma_model_params) {
    gamma_full <- model_for_params$mean_fn(gamma_model_params)
    jacobian_full <- model$jacobian_fn(gamma_full)
    
    jacobian_model_params <- model_for_params$jacobian_fn(gamma_model_params)
    
    jacobian_full %*% jacobian_model_params
  }
  
  new_model(
    list(
      mean_fn = mean_fn_model_params,
      jacobian_fn = jacobian_fn_model_params,
      nuisance_params_position = model_for_params$nuisance_params_position,
      treatment_params_null = model_for_params$treatment_params_null
    )
  )
}

validate_model <- function(model) {
  required_fields <- c("mean_fn", "jacobian_fn", "nuisance_params_position", "treatment_params_null")
  
  missing_fields <- setdiff(required_fields, names(model))
  if (length(missing_fields) > 0) {
    stop("Model is missing required fields: ", paste(missing_fields, collapse = ", "))
  }
  
  if (!is.function(model$mean_fn)) {
    stop("'mean_fn' must be a function.")
  }
  
  if (!is.function(model$jacobian_fn)) {
    stop("'jacobian_fn' must be a function.")
  }
  
  if (!is.numeric(model$nuisance_params_position) || any(model$nuisance_params_position <= 0)) {
    stop("'nuisance_params_position' must be a numeric vector of positive integers.")
  }
  
  if (!is.numeric(model$treatment_params_null)) {
    stop("'treatment_params_null' must be a numeric vector.")
  }
  
  TRUE
}

