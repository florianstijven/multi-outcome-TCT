# Second-stage GLS: fit working-model parameters to first-stage estimates,
# build the plug-in contrast matrix B_n, and run the targeted test.
#
# All functions operate on the stacked mean vector / covariance convention
# established in helper-functions-local-power.R.  The stacked mean vector
# m_tilde is ordered as:
#   (ctrl_j1, trt_j1, ctrl_j2, trt_j2, ..., ctrl_jJ, trt_jJ)
# where each block has K+1 entries (baseline + K follow-up visits).


# ============================================================================
# GLS objective
# ============================================================================

# Q(gamma0; m_tilde, Sigma_inv) for the null-constrained optimisation.
# gamma0_vec  : concatenation of reference-trajectory parameters for all J
#               outcomes (nuisance parameters; treatment-effect parameters are
#               held at their null values).
# mean_fn     : function(gamma0_vec) -> stacked mean vector mu(gamma0, gamma1_null)
# m_tilde     : first-stage stacked mean estimate
# Sigma_inv   : pre-computed solve(Sigma_n)
gls_criterion <- function(gamma0_vec, mean_fn, m_tilde, Sigma_inv) {
  diff <- m_tilde - mean_fn(gamma0_vec)
  as.numeric(crossprod(diff, Sigma_inv %*% diff))
}

# Evaluate the reference trajectory at given times.
ref_eval <- function(times, ref, params, ...) {
  if (ref == "4PL") {
    function_4PL(times, params)
  } else if (ref == "nc_spline") {
    do.call(function_nc_spline, c(list(times = times, params = params), list(...)))
  } else {
    stop("Unknown reference model: ", ref)
  }
}


# ============================================================================
# Null-constrained GLS estimator
# ============================================================================

# Minimise Q over reference-trajectory parameters gamma0 while holding
# gamma1 = gamma1_null.
#
# Returns a list with:
#   $gamma0_hat  : estimated reference-trajectory parameters (vector)
#   $criterion   : value of Q at the optimum
#   $optim       : raw optim() output
#
# m_tilde     : first-stage stacked mean estimate
# Sigma       : first-stage covariance estimate
# mean_fn     : null mean function from make_null_mean_fn_*
# start       : starting values for gamma0
# method      : optimisation method passed to optim(); default "BFGS"

fit_gls <- function(m_tilde, Sigma, mean_fn, jac_mean_fn, start, method = "BFGS") {
  Sigma_inv <- solve(Sigma)
  gls_gradient <- function(gamma0_vec, ...) {
    mu_hat <- mean_fn(gamma0_vec)
    jacobian <- as.matrix(jac_mean_fn(gamma0_vec))
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

# Function to process the results from the GLS fit. Mainly to annotate the
# estimates.
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

# Split the estimated parameters into two lists of estimated parameters. The
# first list is a list where each element are the estimated reference-trajectory
# parameters for each outcome. The second list is a list where each element are
# the estimated slowing parameters for each outcome.
fit_gls_split_gamma_hat <- function(gamma_hat, null_model, slowing_models) {
  J <- length(slowing_models)
  
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



two_stage_gls_null <- function(m_tilde, Sigma, slowing_models, start) {
  # Build the null mean function
  mean_fn_null <- mu_from_gamma_null_f_factory(slowing_models)

  jac_mean_fn_null <- jac_mu_from_gamma_null_f_factory(slowing_models)
  
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
    mean_fn = mean_fn_null,
    jac_mean_fn = jac_mean_fn_null,
    start = start
  )
  
  # Annotate the estimates with names
  gamma_hat_named <- fit_gls_est_names(gls_fitted$gamma_hat,
                                       null_model = TRUE,
                                       slowing_models = slowing_models)
  
  gamma_hat_list <- fit_gls_split_gamma_hat(gls_fitted$gamma_hat,
                                            null_model = TRUE,
                                            slowing_models = slowing_models)
  
  list(
    gamma_hat = gamma_hat_named,
    gamma_hat_list = gamma_hat_list,
    criterion  = gls_fitted$criterion,
    optim      = gls_fitted$optim
  )
}

two_stage_gls_full <- function(m_tilde, Sigma, slowing_models) {
  # Build the full mean function
  mean_fn_full <- mu_from_gamma_f_factory(slowing_models)
  
  # Starting values for gamma0 and gamma1
  start <- unlist(lapply(slowing_models, function(model) {
    c(
      seq_len(model$reference_trajectory_functions_list$no_params),
      model$null_gamma1
    )
  }))
  # Fit the GLS model under the full model
  gls_fitted <- fit_gls(
    m_tilde = m_tilde,
    Sigma = Sigma,
    mean_fn = mean_fn_full,
    start = start
  )
  
  # Annotate the estimates with names
  gamma_hat_named <- fit_gls_est_names(gls_fitted$gamma_hat,
                                       null_model = FALSE,
                                       slowing_models = slowing_models)
  
  gamma_hat_list <- fit_gls_split_gamma_hat(gls_fitted$gamma_hat,
                                            null_model = FALSE,
                                            slowing_models = slowing_models)
  
  list(
    gamma_hat = gamma_hat_named,
    gamma_hat_list = gamma_hat_list,
    criterion  = gls_fitted$criterion,
    optim      = gls_fitted$optim
  )
}





