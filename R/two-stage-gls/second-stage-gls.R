# Second-stage GLS: fit working-model parameters to first-stage estimates, build
# the plug-in contrast matrix B_n, and run the targeted test.
#
# All functions operate on the stacked mean vector / covariance convention
# Evaluate the GLS criterion function 
#' 
#' `gls_criterion()` evaluates the GLS criterion function 
#' Q(gamma) = (m_tilde - mu(gamma))' Sigma^{-1} (m_tilde - mu(gamma)).
#' 
#' The model is specified by the mean function `mean_fn`, which takes a vector
#' of parameters `gamma` and returns the stacked mean vector. The first-stage 
#' estimates are given by `m_tilde` and `Sigma_inv` is the inverse of the 
#' first-stage covariance matrix.
#' 
#' @param gamma (numeric) Vector of parameters for the mean function.
#' @param mean_fn (function) that takes a vector of parameters (`gamma`) and 
#' returns the stacked mean vector.
#' @param m_tilde (numeric) Stacked mean vector from the first-stage estimates.
#' @param Sigma_inv (matrix) Inverse of the first-stage covariance matrix.
#'
#' @returns (numeric) The value of the GLS criterion function evaluated at `gamma`.
gls_criterion <- function(gamma, mean_fn, m_tilde, Sigma_inv) {
  diff <- m_tilde - mean_fn(gamma)
  as.numeric(crossprod(diff, Sigma_inv %*% diff))
}

#' Fit a GLS model to first-stage estimates
#'
#' `fit_gls()` fits a generalized least squares (GLS) model to the first-stage
#' estimates `m_tilde` and `Sigma`. The model is specified by the mean function
#' `mean_fn` and its Jacobian `jac_mean_fn`. The optimization is performed using
#' the `optim()` function, starting from the initial parameter values provided
#' in `start`. The optimization method can be specified using the `method`
#' argument.
#'
#' @param m_tilde (numeric) Stacked mean vector from the first-stage estimates.
#' @param Sigma (matrix) Covariance matrix from the first-stage estimates.
#' @param mean_fn (function) that takes a vector of parameters (`gamma`) and
#'   returns the stacked mean vector. This argument defines the model.
#' @param jac_mean_fn (function) that takes a vector of parameters (`gamma`) and
#'   returns the Jacobian of the mean function. This is used to compute the
#'   gradient of the GLS criterion function.
#' @param start (numeric) Starting values for the parameters to be estimated.
#' @param method (character) Optimization method to be used in `optim()`.
#'   Default is `"BFGS"`.
#'
#' @returns (list) A list containing:
#'  - `gamma_hat`: Estimated parameters (vector).
#'  - `criterion`: Value of the GLS criterion function at the optimum.
#'  - `optim`: Raw output from the `optim()` function.
fit_gls <- function(m_tilde, Sigma, mean_fn, jac_mean_fn, start, method = "BFGS") {
  Sigma_inv <- solve(Sigma)
  gls_gradient <- function(gamma, ...) {
    mu_hat <- mean_fn(gamma)
    jacobian <- as.matrix(jac_mean_fn(gamma))
    diff <- m_tilde - mu_hat
    -2 * t(jacobian) %*% Sigma_inv %*% diff
  }
  
  opt <- stats::optim(
    par     = start,
    fn      = gls_criterion,
    gr      = gls_gradient,
    mean_fn = mean_fn,
    m_tilde = m_tilde,
    Sigma_inv = Sigma_inv,
    method  = method,
    control = list(maxit = 1e3)
  )
  
  # Raise a warning if the optimization did not converge.
  if (opt$convergence != 0) {
    warning("GLS optimization did not converge. Convergence code: ", opt$convergence)
  }
  
  list(
    gamma_hat = opt$par,
    criterion  = opt$value,
    optim      = opt
  )
}

#' Process GLS fit estimates and assign names based on the model structure
#'
#' @param gamma_hat (numeric) Vector of estimated parameters from the GLS fit.
#' @param null_model (boolean) Indicates whether the model is a null model
#'   (TRUE) where the treatment effect is zero or a full model (FALSE).
#' @param slowing_models (list) A list of slowing model specifications for each
#'   outcome, as returned by `build_slowing_models()`.
#'
#' @returns (named numeric) A named vector of estimated parameters, where the
#'   names are constructed based on the model structure and outcome indices.
fit_gls_est_names <- function(gamma_hat, null_model, slowing_models) {
  # Names vector for gamma_hat
  names_vec <- c()
  
  # Number of outcomes
  J <- length(slowing_models)
  
  gamma0_hat <- c()
  
  for (outcome_idx in seq_len(J)) {
    no_params_ref <- slowing_models[[outcome_idx]]$reference_trajectory_functions_list$no_params
    ref_model <- slowing_models[[outcome_idx]]$reference_trajectory_functions_list$ref
    
    names_gamma0_hat <- paste0("gamma0_",
                               ref_model,
                               "_outcome",
                               outcome_idx,
                               "_param",
                               seq_len(no_params_ref))
    names_vec <- c(names_vec, names_gamma0_hat)
    
    if (!null_model) {
      no_params_slowing <- length(slowing_models[[outcome_idx]]$null_gamma1)
      type <- slowing_models[[outcome_idx]]$type
      
      names_gamma1_hat <- paste0(
        "gamma1_",
        type,
        "_outcome",
        outcome_idx,
        "_param",
        seq_len(no_params_slowing)
      )
      names_vec <- c(names_vec, names_gamma1_hat)
    }
  }
  return(setNames(gamma_hat, names_vec))
}


#' Split Estimated Slowing Parameters
#'
#' This function takes a vector of estimated parameters (`gamma_hat`) from a GLS
#' fit and splits it into two lists: one for the estimated reference-trajectory
#' parameters (`gamma0_hat`) and another for the estimated slowing parameters
#' (`gamma1_hat`). The splitting is based on the structure of the provided
#' `slowing_models` and whether the model is a null model or not.
#'
#' @inheritParams fit_gls_est_names
#'
#' @returns (list) A list containing two elements:
#' - `gamma0_hat`: A list of estimated reference-trajectory parameters for each outcome.
#' - `gamma1_hat`: A list of estimated slowing parameters for each outcome.
fit_gls_split_gamma_hat <- function(gamma_hat, null_model, working_model) {

  gamma0_hat_list <- vector("list", J)
  gamma1_hat_list <- vector("list", J)
  
  start_idx <- 1
  for (outcome_idx in seq_len(J)) {
    no_params_ref <- slowing_models[[outcome_idx]]$reference_trajectory_functions_list$no_params
    end_idx <- start_idx + no_params_ref - 1
    gamma0_hat_list[[outcome_idx]] <- gamma_hat[start_idx:end_idx]
    
    no_params_slowing <- length(slowing_models[[outcome_idx]]$null_gamma1)
    if (!null_model) {
      start_idx <- end_idx + 1
      end_idx <- start_idx + no_params_slowing - 1
      gamma1_hat_list[[outcome_idx]] <- gamma_hat[start_idx:end_idx]
    } else {
      gamma1_hat_list[[outcome_idx]] <- slowing_models[[outcome_idx]]$null_gamma1
    }
    
    start_idx <- end_idx + 1
  }
  
  list(
    gamma0_hat = gamma0_hat_list,
    gamma1_hat = gamma1_hat_list
  )
}



#' Fit a slowing model using GLS
#'
#' [two_stage_gls_full()] fits a generalized least squares (GLS) model to the
#' first-stage estimates. [two_stage_gls_null()] fits the same model under the
#' null of no treatment effect. It constructs the mean function based on
#' the provided slowing models, and then optimizes the GLS criterion function to
#' estimate the reference-trajectory parameters (and slowing parameters for [two_stage_gls_full()]).
#'
#' @inheritParams fit_gls
#' @inheritParams fit_gls_split_gamma_hat
#'
#' @returns (list) A list containing:
#' - `gamma_hat`: Named vector of estimated parameters under the null.
#' - `gamma_hat_list`: List of estimated parameters split into reference-trajectory and slowing
#'   parameters for each outcome.
#'  - `criterion`: Value of the GLS criterion function at the optimum.
#'  - `optim`: Raw output from the `optim()` function.
two_stage_gls_null <- function(m_tilde, Sigma, working_model, start) {
  if (!inherits(working_model, "model")) {
    stop("Object is not of class 'model'.")
  } else {
    validate_model(working_model)
  }
  
  # Starting values for gamma0
  if (is.null(start)) {
    start <- rep(1, sum(
      sapply(slowing_models, function(model)
        model$reference_trajectory_functions_list$no_params)
    ))
  }
  
  # Fit the GLS model under the null
  gls_fitted <- fit_gls(
    m_tilde = m_tilde,
    Sigma = Sigma,
    mean_fn = working_model$mean_fn_null,
    jac_mean_fn = working_model$jacobian_fn_null,
    start = start
  )
  
  list(
    gamma_hat = gls_fitted$gamma_hat,
    criterion  = gls_fitted$criterion,
    optim      = gls_fitted$optim
  )
}


#' @rdname two_stage_gls_null 
two_stage_gls_full <- function(m_tilde, Sigma, working_model, start) {
  
  # Fit the GLS model under the full model
  gls_fitted <- fit_gls(
    m_tilde = m_tilde,
    Sigma = Sigma,
    mean_fn = working_model$mean_fn,
    start = start
  )
  
  list(
    gamma_hat = gls_fitted$gamma_hat,
    criterion  = gls_fitted$criterion,
    optim      = gls_fitted$optim
  )
}





