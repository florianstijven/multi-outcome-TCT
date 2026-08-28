#' Finalize a working-model object
#'
#' Given a list with a mean function, Jacobian function, and information about
#' which parameters are nuisance vs. treatment-effect parameters, computes the
#' derived null/treatment mean and Jacobian functions and assigns the `model`
#' class.
#'
#' @param x list with elements `mean_fn`, `jacobian_fn`,
#'   `nuisance_params_position`, and `treatment_params_null`.
#'
#' @returns an object of class `model`.
#' @export
new_model <- function(x = list()) {
  stopifnot(is.list(x))
  
  # Compute location of treatment-effect parameters in the full parameter
  # vector.
  x$no_of_params <- length(x$nuisance_params_position) + length(x$treatment_params_null)
  x$treatment_params_position <- seq_len(x$no_of_params)[-x$nuisance_params_position]
  
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

#' Build the null mean function for a model
#'
#' @param mean_fn function, the full mean function.
#' @param nuisance_params_position integer vector, positions of the nuisance
#'   parameters in the full parameter vector.
#' @param treatment_params_null numeric vector, the null values of the
#'   treatment-effect parameters.
#'
#' @returns a function of `gamma0` (the nuisance parameters) returning the
#'   mean vector at the null.
#' @keywords internal
mean_fn_null_constructor <- function(mean_fn, nuisance_params_position, treatment_params_null) {
  function(gamma0) {
    gamma <- numeric(length = length(nuisance_params_position) + length(treatment_params_null))
    gamma[nuisance_params_position] <- gamma0
    gamma[-nuisance_params_position] <- treatment_params_null
    
    mean_fn(gamma)
  }
}

#' Build the null Jacobian function for a model
#'
#' @inheritParams mean_fn_null_constructor
#' @param jacobian_fn function, the full Jacobian function.
#'
#' @returns a function of `gamma0` returning the Jacobian with respect to the
#'   nuisance parameters, evaluated at the null.
#' @keywords internal
jacobian_fn_null_constructor <- function(jacobian_fn, nuisance_params_position, treatment_params_null) {
  function(gamma0) {
    gamma <- numeric(length = length(nuisance_params_position) + length(treatment_params_null))
    gamma[nuisance_params_position] <- gamma0
    gamma[-nuisance_params_position] <- treatment_params_null

    jacobian <- jacobian_fn(gamma)

    jacobian[, nuisance_params_position]
  }
}

#' Build the treatment-effect Jacobian function at the null
#'
#' @inheritParams jacobian_fn_null_constructor
#'
#' @returns a function of `gamma0` returning the Jacobian with respect to the
#'   treatment-effect parameters, evaluated at the null.
#' @keywords internal
jacobian_fn_treatment_null_constructor <- function(jacobian_fn, nuisance_params_position, treatment_params_null) {
  function(gamma0) {
    gamma <- numeric(length = length(nuisance_params_position) + length(treatment_params_null))
    gamma[nuisance_params_position] <- gamma0
    gamma[-nuisance_params_position] <- treatment_params_null
    
    jacobian <- jacobian_fn(gamma)
    
    jacobian[, -nuisance_params_position, drop = FALSE]
  }
}

#' Construct a working-model object
#'
#' Constructs an object of class `model` from a mean function and, optionally,
#' its Jacobian. If `jacobian_fn` is not supplied, it is approximated via
#' numerical differentiation with [numDeriv::jacobian()].
#'
#' @param mean_fn function, the mean function.
#' @param jacobian_fn optional function, the Jacobian of `mean_fn`. If
#'   `NULL`, computed numerically.
#' @param nuisance_params_position integer vector, positions of the nuisance
#'   parameters in the full parameter vector.
#' @param treatment_params_null numeric vector, the null values of the
#'   treatment-effect parameters.
#'
#' @returns an object of class `model`.
#' @export
model <- function(mean_fn, jacobian_fn = NULL, nuisance_params_position, treatment_params_null, param_names = NULL) {
  # If the Jacobian function is not provided, we approximate it using numerical
  # differentiation.
  if (is.null(jacobian_fn)) {
    jacobian_fn <- function(gamma) {
      numDeriv::jacobian(func = mean_fn, x = gamma)
    }
  }
  
  new_model(
    list(
      mean_fn = mean_fn,
      jacobian_fn = jacobian_fn,
      nuisance_params_position = nuisance_params_position,
      treatment_params_null = treatment_params_null,
      param_names = param_names
    )
  )
}

#' Concatenate several models into one
#'
#' Stacks the mean functions and (block-diagonal) Jacobians of several
#' `model` objects into a single combined `model`, used to build one
#' multi-outcome working model out of per-outcome models.
#'
#' @param models list of objects of class `model`.
#'
#' @returns an object of class `model` combining all elements of `models`.
#' @export
concatenate_models <- function(models, submodel_names = NULL) {
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
              .f =  function(working_model, start_index) {
                end_index <- start_index + length(working_model$nuisance_params_position) + length(working_model$treatment_params_null) - 1
                working_model$mean_fn(gamma[start_index:end_index])
              }
            ))
  }
  
  jacobian_fn_combined <- function(gamma) {
    do.call(Matrix::bdiag, purrr::map2(.x = models, .y = start_indices_full, .f = function(working_model, start_index) {
      end_index <- start_index + length(working_model$nuisance_params_position) + length(working_model$treatment_params_null) - 1
      working_model$jacobian_fn(gamma[start_index:end_index])
    }))
  }
  
  nuisance_params_position_combined <- unlist(purrr::map2(.x = models, .y = start_indices_full, function(working_model, start_index) working_model$nuisance_params_position + start_index - 1))
  treatment_params_null_combined <- unlist(lapply(models, function(working_model) working_model$treatment_params_null))
  
  if (is.null(submodel_names)) {
    submodel_names <- rep(list(NULL), length(models))
  }
  
  param_names_combined <- unlist(purrr::map2(
    .x = models,
    .y = submodel_names,
    .f = function(working_model, name) { 
      extra_name <- if (!is.null(name)) paste0(" (", name, ")") else ""
      paste0(working_model$param_names, extra_name)
      }
  ))
  
  new_model(
    list(
      mean_fn = mean_fn_combined,
      jacobian_fn = jacobian_fn_combined,
      nuisance_params_position = nuisance_params_position_combined,
      treatment_params_null = treatment_params_null_combined,
      param_names = param_names_combined
    )
  )
}

#' Split a working model by parameter/mean-vector stratum
#'
#' Splits `working_model` into one sub-model per unique value of
#' `split_indices_params`, restricting the mean and Jacobian functions to the
#' corresponding rows/columns.
#'
#' @param working_model an object of class `model`.
#' @param split_indices_mu integer vector, the stratum of each element of the
#'   mean vector.
#' @param split_indices_params list of integers, the stratum of each
#'   parameter (elements may have length greater than one for shared-parameter
#'   models).
#'
#' @returns list of objects of class `model`, one per unique stratum.
#' @export
split_model <- function(working_model, split_indices_mu, split_indices_params) {
  if (!inherits(working_model, "model")) {
    stop("'working_model' must be of class 'model'.")
  }
  # 
  # if (!is.numeric(split_indices) || any(split_indices <= 0)) {
  #   stop("'split_indices' must be a numeric vector of positive integers.")
  # }
  
  # Split the model into sub-models based on the provided indices.
  unique_split_indices <- unique(unlist(split_indices_params))
  n_models <- length(unique_split_indices)
  
  models <- rep(list(NULL), n_models)
  start_index <- 1
  
  for (i in seq_along(unique_split_indices)) {
    params_match <- sapply(
      split_indices_params,
      function(x) unique_split_indices[i] %in% x
    )

    mean_fn_sub <- function(gamma_sub) {
      gamma <- rep(0, working_model$no_of_params)
      gamma[params_match] <- gamma_sub
      working_model$mean_fn(gamma)[split_indices_mu == unique_split_indices[i]]
    }
    jacobian_fn_sub <- function(gamma_sub) {
      gamma <- rep(0, working_model$no_of_params)
      gamma[params_match] <- gamma_sub
      working_model$jacobian_fn(gamma)[split_indices_mu == unique_split_indices[i], params_match, drop = FALSE]
    }
    
    params_position_sub <- which(params_match)
    nuisance_params_position_sub <- which(params_position_sub %in% working_model$nuisance_params_position)
    treatment_params_null_sub_position <- which(
      working_model$treatment_params_position %in% which(
        params_match &
          seq_len(working_model$no_of_params) %in% working_model$treatment_params_position
      )
    )
    treatment_params_null_sub <- working_model$treatment_params_null[treatment_params_null_sub_position]
    
    param_names_sub <- working_model$param_names[params_position_sub]
    
    models[[i]] <- model(
      mean_fn = mean_fn_sub,
      jacobian_fn = jacobian_fn_sub,
      nuisance_params_position = nuisance_params_position_sub,
      treatment_params_null = treatment_params_null_sub,
      param_names = param_names_sub
    )
  }
  
  models
}

#' Reparameterize a model to share treatment-effect parameters
#'
#' Wraps `working_model` in a reduced parameterization in which the
#' treatment-effect parameters at the positions listed in
#' `shared_param_positions` are constrained to be equal.
#'
#' @param working_model an object of class `model`.
#' @param shared_param_positions numeric vector of positive integers (or a
#'   list of such vectors, one per shared parameter), indicating positions in
#'   the full parameter vector that should share a single value.
#'
#' @returns an object of class `model` with a reduced parameter vector.
#' @export
shared_parameter_model <- function(working_model, shared_param_positions) {
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
  
  # Construct a model for the parameters of `working_model`, taking into account the
  # shared parameters. The new parameter vector will have one entry for each
  # shared parameter, and one entry for each non-shared parameter. The mean
  # function of the new model will map the new parameter vector to the full
  # parameter vector of `working_model`.
  no_of_non_shared_params <- length(working_model$nuisance_params_position) + length(working_model$treatment_params_null) - length(unlist(shared_param_positions))
  no_of_non_shared_treatment_params <- length(working_model$treatment_params_null) - length(unlist(shared_param_positions))
  
  mean_fn_params <- function(gamma_shared) {
    # Construct the full parameter vector for `working_model` from the parameters of
    # `gamma_shared`, taking into account the shared parameters.
    gamma_full <- numeric(length = working_model$no_of_params)
    gamma_full[-unlist(shared_param_positions)] <- gamma_shared[seq_len(no_of_non_shared_params)]
    
    # Fill in the shared parameters.
    for (i in seq_along(shared_param_positions)) {
      shared_param_positions_i <- shared_param_positions[[i]]
      gamma_full[shared_param_positions_i] <- gamma_shared[no_of_non_shared_params + i]
    }
    gamma_full
  }
  
  jacobian_fn_params <- function(gamma_shared) {
    jacobian_full <- matrix(0, nrow = working_model$no_of_params, ncol = length(gamma_shared))
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
  # `length(working_model$nuisance_params_position)` entries of the parameter vector.
  nuisance_params_position = seq_len(length(working_model$nuisance_params_position))
  
  # Determine the null values for the non-shared treatment-effect parameters.
  non_shared_treatment_params_position_complement <- c(working_model$nuisance_params_position, unlist(shared_param_positions))
  non_shared_treatment_params_position <- seq_len(working_model$no_of_params)[-non_shared_treatment_params_position_complement]
  # Relative position of the non-shared treatment-effect parameters in the
  # original treatment-effect parameter vector.
  non_shared_treatment_params_position_subset <- sapply(
    non_shared_treatment_params_position,
    function(pos) {
      which.max(working_model$treatment_params_position == pos)
    }
  )
  shared_treatment_params_position_subset <- 
    sapply(
      shared_param_positions,
      function(positions) {
        which.max(working_model$treatment_params_position %in% positions)
      }
    )
  
  if (length(non_shared_treatment_params_position_subset) == 0) {
    non_shared_treatment_params_position_subset <- integer(0)
  }

  treatment_params_null <- numeric(length = no_of_non_shared_treatment_params + length(shared_param_positions))
  
  treatment_params_null[seq_len(no_of_non_shared_treatment_params)] <- working_model$treatment_params_null[non_shared_treatment_params_position_subset]
  treatment_params_null[1:length(shared_param_positions) + no_of_non_shared_treatment_params] <- working_model$treatment_params_null[shared_treatment_params_position_subset]
  
  model_for_params <- new_model(
    list(
      mean_fn = mean_fn_params,
      jacobian_fn = jacobian_fn_params,
      nuisance_params_position = nuisance_params_position,
      treatment_params_null = treatment_params_null,
      param_names = c(
        working_model$param_names[-unlist(shared_param_positions)],
        sapply(seq_along(shared_param_positions), function(i) working_model$param_names[shared_param_positions[[i]][1]])
      )
    )
  )
  model_for_model_params(working_model, model_for_params)
}

#' Compose a model with its reduced-parameter model
#'
#' @param working_model an object of class `model`.
#' @param model_for_params an object of class `model` mapping the reduced
#'   (e.g. shared) parameter vector onto `working_model`'s full parameter
#'   vector.
#'
#' @returns an object of class `model` parameterized by
#'   `model_for_params`'s parameters.
#' @keywords internal
model_for_model_params <- function(working_model, model_for_params) {
  mean_fn_model_params <- function(gamma_model_params) {
    gamma_full <- model_for_params$mean_fn(gamma_model_params)
    working_model$mean_fn(gamma_full)
  }
  
  jacobian_fn_model_params <- function(gamma_model_params) {
    gamma_full <- model_for_params$mean_fn(gamma_model_params)
    jacobian_full <- working_model$jacobian_fn(gamma_full)
    
    jacobian_model_params <- model_for_params$jacobian_fn(gamma_model_params)
    
    jacobian_full %*% jacobian_model_params
  }
  
  new_model(
    list(
      mean_fn = mean_fn_model_params,
      jacobian_fn = jacobian_fn_model_params,
      nuisance_params_position = model_for_params$nuisance_params_position,
      treatment_params_null = model_for_params$treatment_params_null,
      param_names = model_for_params$param_names
    )
  )
}

#' Validate a working-model object
#'
#' Checks that `working_model` has the required fields (`mean_fn`,
#' `jacobian_fn`, `nuisance_params_position`, `treatment_params_null`) with
#' the expected types.
#'
#' @param working_model an object of class `model`.
#'
#' @returns `TRUE`, invisibly on success; raises an error otherwise.
#' @export
validate_model <- function(working_model) {
  required_fields <- c("mean_fn", "jacobian_fn", "nuisance_params_position", "treatment_params_null")
  
  missing_fields <- setdiff(required_fields, names(working_model))
  if (length(missing_fields) > 0) {
    stop("Model is missing required fields: ", paste(missing_fields, collapse = ", "))
  }
  
  if (!is.function(working_model$mean_fn)) {
    stop("'mean_fn' must be a function.")
  }
  
  if (!is.function(working_model$jacobian_fn)) {
    stop("'jacobian_fn' must be a function.")
  }
  
  if (!is.numeric(working_model$nuisance_params_position) || any(working_model$nuisance_params_position <= 0)) {
    stop("'nuisance_params_position' must be a numeric vector of positive integers.")
  }
  
  if (!is.numeric(working_model$treatment_params_null)) {
    stop("'treatment_params_null' must be a numeric vector.")
  }
  
  TRUE
}

#' Print a working-model object
#'
#' @param x an object of class `model`.
#' @param ... unused; retained for S3 method consistency.
#'
#' @returns `x`, invisibly.
#' @export
print.model <- function(x, ...) {
  cat("Model object:\n")
  cat("Number of parameters:", x$no_of_params, "\n")
  cat("Nuisance parameters positions:", paste(x$nuisance_params_position, collapse = ", "), "\n")
  cat("Treatment parameters null values:", paste(x$treatment_params_null, collapse = ", "), "\n")
}
