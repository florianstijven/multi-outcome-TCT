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
  function(gamma) {
    gamma <- numeric(length = length(nuisance_params_position) + length(treatment_params_null))
    gamma[nuisance_params_position] <- gamma0
    gamma[-nuisance_params_position] <- treatment_params_null
    
    jacobian <- jacobian_fn(gamma)
    
    jacobian[, -nuisance_params_position]
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
  
  # mean_fn_null_combined <- function(gamma0) {
  #   do.call(
  #     c,
  #     purrr::map2(.x = models, .y = start_indices_null),
  #     .f =  function(model, start_index) {
  #       end_index <- start_index + length(model$nuisance_params_position) - 1
  #       model$mean_fn_null(gamma0[start_index:end_index])
  #     }
  #   )
  # }
  
  jacobian_fn_combined <- function(gamma) {
    do.call(Matrix::bdiag, purrr::map2(.x = models, .y = start_indices_full, .f = function(model, start_index) {
      end_index <- start_index + length(model$nuisance_params_position) + length(model$treatment_params_null) - 1
      model$jacobian_fn(gamma[start_index:end_index])
    }))
  }
  
  # jacobian_fn_null_combined <- function(gamma0) {
  #   do.call(Matrix::bdiag, purrr::map2(.x = models, .y = start_indices_null, .f = function(model, start_index) {
  #     end_index <- start_index + length(model$nuisance_params_position) - 1
  #     model$jacobian_fn_null(gamma0[start_index:end_index])
  #   }))
  # }
  # 
  # jacobian_fn_treatment_null_combined <- function(gamma0) {
  #   do.call(Matrix::bdiag, purrr::map2(.x = models, .y = start_indices_null, .f = function(model, start_index) {
  #     end_index <- start_index + length(model$nuisance_params_position) - 1
  #     model$jacobian_fn_treatment_null(gamma0[start_index:end_index])
  #   }))
  # }
  
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

