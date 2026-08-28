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
#' @export
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
#' @export
fit_gls <- function(m_tilde, Sigma, mean_fn, jac_mean_fn, start, method = "L-BFGS-B") {
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
#' @export
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
#' @export
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
#' `two_stage_gls_full()` fits a generalized least squares (GLS) model to the
#' first-stage estimates. `two_stage_gls_null()` fits the same model under the
#' null of no treatment effect. It constructs the mean function based on the
#' provided slowing models, and then optimizes the GLS criterion function to
#' estimate the reference-trajectory parameters (and slowing parameters for
#' `two_stage_gls_full()`).
#'
#' @inheritParams fit_gls
#' @inheritParams fit_gls_split_gamma_hat
#' @param split_indices_mu integer vector where each element denotes the
#'   stratum of the the corresponding element in `m_tilde`.
#' @param split_indices_params list of integers where each element denotes the stratum of the
#'   the corresponding parameter in parameter vector for `working_model`. The
#'   list elements can have length greater than one, which may be required for
#'   shared-parameter models.
#'
#' @returns an object of class `model_fit` (see [new_model_fit()]).
#' @export
two_stage_gls_null <- function(m_tilde,
                               Sigma,
                               working_model,
                               start,
                               ols = FALSE,
                               split_indices_mu = NULL,
                               split_indices_params = NULL) {
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
  
  # OLS is equivalent to GLS with an identity covariance matrix. If ols = TRUE,
  # we replace Sigma with an identity matrix.
  if (ols) {
    Sigma_new <- diag(1, nrow = length(m_tilde))
  } else{
    Sigma_new <- Sigma
  }
  
  if (is.null(split_indices_mu)) {
    # Fit the GLS model under the null
    split_indices_mu <- rep(1, length(m_tilde))
    n_params <- length(working_model$nuisance_params_position) + length(working_model$treatment_params_null)
    split_indices_params <- rep(1, n_params)
  }
  split_indices_params_null <- split_indices_params[working_model$nuisance_params_position]
  
  # Fit the GLS model by stratum.
  unique_strata <- unique(split_indices_mu)
  n_strata <- length(unique_strata)
  
  m_tilde_list <- rep(list(NULL), n_strata)
  Sigma_list <- rep(list(NULL), n_strata)
  Sigma_true_list <- rep(list(NULL), n_strata)
  start_list <- rep(list(NULL), n_strata)
  for (i in seq_along(unique_strata)) {
    stratum_idx <- split_indices_mu == unique_strata[i]
    m_tilde_list[[i]] <- m_tilde[stratum_idx]
    Sigma_list[[i]] <- Sigma_new[stratum_idx, stratum_idx]
    Sigma_true_list[[i]] <- Sigma[stratum_idx, stratum_idx]
    start_list[[i]] <- start[split_indices_params_null == unique_strata[i]]
  }
  models_list <- split_model(
    working_model = working_model,
    split_indices_mu = split_indices_mu,
    split_indices_params = split_indices_params
  )
  
  gls_fitted_list <- lapply(seq_along(m_tilde_list), function(i) {
    fit_gls(
      m_tilde = m_tilde_list[[i]],
      Sigma = Sigma_list[[i]],
      mean_fn = models_list[[i]]$mean_fn_null,
      jac_mean_fn = models_list[[i]]$jacobian_fn_null,
      start = start_list[[i]]
    )
  })
  
  # The GLS weight used for fitting (Sigma_list) may differ from the true
  # first-stage covariance (Sigma_true_list) when ols = TRUE; the sandwich
  # formula in gls_vcov() accounts for that.
  vcov_list <- lapply(seq_along(gls_fitted_list), function(i) {
    jacobian_i <- models_list[[i]]$jacobian_fn_null(gls_fitted_list[[i]]$gamma_hat)
    gls_vcov(jacobian = jacobian_i, Sigma = Sigma_true_list[[i]], Sigma_working = Sigma_list[[i]])
  })
  
  coefficients <- purrr::map(gls_fitted_list, "gamma_hat") %>% unlist()
  names(coefficients) <- working_model$param_names[working_model$nuisance_params_position]
  
  new_model_fit(
    coefficients  = coefficients,
    vcov          = as.matrix(Matrix::bdiag(vcov_list)),
    criterion     = sum(purrr::map_dbl(gls_fitted_list, "criterion")),
    optim         = purrr::map(gls_fitted_list, "optim"),
    working_model = working_model,
    null_model    = TRUE,
    data          = list(m_tilde = m_tilde, Sigma = Sigma)
  )
}


#' @rdname two_stage_gls_null 
#' @export
two_stage_gls_full <- function(m_tilde,
                               Sigma,
                               working_model,
                               start,
                               ols = FALSE,
                               split_indices_mu = NULL,
                               split_indices_params = NULL) {
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
  
  # OLS is equivalent to GLS with an identity covariance matrix. If ols = TRUE,
  # we replace Sigma with an identity matrix.
  if (ols) {
    Sigma_new <- diag(1, nrow = length(m_tilde))
  } else{
    Sigma_new <- Sigma
  }
  
  if (is.null(split_indices_mu)) {
    # Fit the GLS model under the null
    split_indices_mu <- rep(1, length(m_tilde))
    n_params <- length(working_model$nuisance_params_position) + length(working_model$treatment_params_null)
    split_indices_params <- rep(1, n_params)
  }
  
  # Fit the GLS model by stratum.
  unique_strata <- unique(split_indices_mu)
  n_strata <- length(unique_strata)
  
  m_tilde_list <- rep(list(NULL), n_strata)
  Sigma_list <- rep(list(NULL), n_strata)
  Sigma_true_list <- rep(list(NULL), n_strata)
  start_list <- rep(list(NULL), n_strata)
  for (i in seq_along(unique_strata)) {
    stratum_idx <- split_indices_mu == unique_strata[i]
    m_tilde_list[[i]] <- m_tilde[stratum_idx]
    Sigma_list[[i]] <- Sigma_new[stratum_idx, stratum_idx]
    Sigma_true_list[[i]] <- Sigma[stratum_idx, stratum_idx]
    params_match <- sapply(
      split_indices_params,
      function(x) unique_strata[i] %in% x
    )
    start_list[[i]] <- start[params_match]
  }
  models_list <- split_model(
    working_model = working_model,
    split_indices_mu = split_indices_mu,
    split_indices_params = split_indices_params
  )
  
  gls_fitted_list <- lapply(seq_along(m_tilde_list), function(i) {
    fit_gls(
      m_tilde = m_tilde_list[[i]],
      Sigma = Sigma_list[[i]],
      mean_fn = models_list[[i]]$mean_fn,
      jac_mean_fn = models_list[[i]]$jacobian_fn,
      start = start_list[[i]]
    )
  })
  
  vcov_list <- lapply(seq_along(gls_fitted_list), function(i) {
    jacobian_i <- models_list[[i]]$jacobian_fn(gls_fitted_list[[i]]$gamma_hat)
    gls_vcov(jacobian = jacobian_i, Sigma = Sigma_true_list[[i]], Sigma_working = Sigma_list[[i]])
  })
  
  coefficients <- purrr::map(gls_fitted_list, "gamma_hat") %>% unlist()
  names(coefficients) <- working_model$param_names
  
  new_model_fit(
    coefficients  = coefficients,
    vcov          = as.matrix(Matrix::bdiag(vcov_list)),
    criterion     = sum(purrr::map_dbl(gls_fitted_list, "criterion")),
    optim         = purrr::map(gls_fitted_list, "optim"),
    working_model = working_model,
    null_model    = FALSE,
    data          = list(m_tilde = m_tilde, Sigma = Sigma)
  )
}


#' Plot GLS-fitted values against first-stage estimates
#'
#' Requires the `ggplot2` package.
#'
#' @param gls_fitted an object of class `model_fit` (see [two_stage_gls_null()]
#'   or [two_stage_gls_full()]).
#' @param treatment_strata optional character/factor vector labeling each
#'   element of `m_tilde` by treatment arm.
#' @param outcome_strata optional character/factor vector labeling each
#'   element of `m_tilde` by outcome.
#' @param times optional numeric vector of time points for each element of
#'   `m_tilde`.
#'
#' @returns a `ggplot` object.
#' @export
plot_gls_fitted <- function(gls_fitted, treatment_strata = NULL, outcome_strata = NULL, times = NULL) {
  # Create a data frame for plotting.
  df_gls_fitted_predictions <- df_gls_fitted(gls_fitted, treatment_strata, outcome_strata, times)
  
  # Check whether ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required for plotting. Please install it using install.packages('ggplot2').")
  }
  
  ggplot2::ggplot(df_gls_fitted_predictions, ggplot2::aes(x = time_points, y = fitted_value)) +
    ggplot2::geom_line(ggplot2::aes(color = as.factor(treatment_strata)), size = 1) +
    ggplot2::geom_point(ggplot2::aes(y = first_stage_estimate, color = treatment_strata), size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = first_stage_estimate - 1.96 * SE_first_stage_estimate,
                                        ymax = first_stage_estimate + 1.96 * SE_first_stage_estimate,
                                        color = treatment_strata), width = 0.1) +
    ggplot2::facet_wrap(~ outcome_strata, scales = "free_y") +
    ggplot2::labs(x = "Time Points", y = "Fitted Value / First Stage Estimate",
                  title = "GLS Fitted Values and First Stage Estimates") +
    ggplot2::theme_minimal()
}

#' Build a data frame of GLS-fitted values and first-stage estimates
#'
#' @inheritParams plot_gls_fitted
#'
#' @returns data.frame with columns `fitted_value`, `first_stage_estimate`,
#'   `treatment_strata`, `outcome_strata`, `time_points`, and
#'   `SE_first_stage_estimate`.
#' @export
df_gls_fitted <- function(gls_fitted, treatment_strata = NULL, outcome_strata = NULL, times = NULL) {
  # Create a data frame for plotting with the following columns: fitted value,
  # first stage estimate, treatment strata, outcome strata, time points, and SEs
  # for the first stage estimates.
  if (gls_fitted$null_model) {
    prediction_f <- gls_fitted$working_model$mean_fn_null
  } else {
    prediction_f <- gls_fitted$working_model$mean_fn
  }
  
  data.frame(
    fitted_value = prediction_f(stats::coef(gls_fitted)),
    first_stage_estimate = gls_fitted$data$m_tilde,
    treatment_strata = treatment_strata,
    outcome_strata = outcome_strata,
    time_points = times,
    SE_first_stage_estimate = sqrt(diag(gls_fitted$data$Sigma))
  )
  
  
}
