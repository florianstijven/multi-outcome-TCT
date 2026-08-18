# ============================================================================
# All-in-one targeted test for proportional slowing
# ============================================================================

# Run the full Algorithm 1 from the manuscript (Section sec:score-test-estimated)
# for the proportional slowing working model.
#
# Returns a list with:
#   $B_n        : plug-in contrast matrix
#   $gamma0_hat : estimated reference-trajectory parameters
#   $statistic  : T_hat_n
#   $df         : degrees of freedom
#   $p_value    : chi-squared p-value
#
# m_tilde     : first-stage stacked mean estimate (length 2*sum_j(K_j+1))
# Sigma       : first-stage covariance estimate
# times       : list of time vectors, one per outcome
# start       : list of starting values for the reference-trajectory parameters (gamma0)
# shared      : if TRUE, share the time-mapping working parameters across outcomes
# working_model : List (to be specified)
# ...         : extra args forwarded to spline helpers

#' Data-Adaptive Targeted Test Based on Slowing Model
#'
#' [targeted_test_proportional_slowing()] estimates the slowing model parameters
#' under the null hypothesis and computes a targeted test statistic for the
#' omnibus null hypothesis of no treatment effect across multiple outcomes.
#'
#' @inheritParams two_stage_gls_null
#' @param shared (boolean) If TRUE, the slowing model parameters are shared
#'   across outcomes; if FALSE, they are outcome-specific.
#'
#' @returns (list) A list containing the following elements:
#'  - `B_n`: The plug-in contrast matrix.
#'  - `gls_fitted_null`: The fitted model under the null hypothesis.
#'  - `statistic`: The computed test statistic.
#'  - `df`: The degrees of freedom for the test statistic.
#'  - `p_value`: The p-value associated with the test statistic.
targeted_test <- function(m_tilde, Sigma, working_model, A, shared_matrix, start) {
  # Step 1: Fit working model under the null.
  gls_fitted_null <- two_stage_gls_null(m_tilde = m_tilde,
                                        Sigma   = Sigma,
                                        working_model = working_model, 
                                        start = start)
  
  # Step 2: build the Jacobian at the null estimate
  jacobian <- working_model$jacobian_fn_treatment_null(
    gamma0 = gls_fitted_null$gamma_hat
  )
  
  # If parameters are shared across outcomes, we need to adjust the Jacobian to
  # account for the shared structure.
  if (is.null(shared_matrix)) {
    shared_matrix <- diag(1, nrow = ncol(jacobian))
  }
  jacobian <- jacobian %*% shared_matrix
  
  B_n <- build_contrast_matrix(jacobian = jacobian,
                               A = A,
                               Sigma = Sigma)
  
  # Step 4: evaluate the test statistic
  test <- targeted_test_statistic(B = B_n,
                                  m_tilde = m_tilde,
                                  Sigma = Sigma)
  
  list(
    B_n        = B_n,
    gls_fitted_null = gls_fitted_null,
    statistic  = test$statistic,
    df         = test$df,
    p_value    = test$p_value
  )
}


#' #' Build Jacobian of Slowing Model
#' #'
#' #' @param gamma0 (list) A list of reference-trajectory parameters under the null
#' #'   hypothesis. Each element of the list corresponds to an outcome.
#' #' @inheritParams two_stage_gls_null
#' #'
#' #' @returns (matrix) The Jacobian matrix of the slowing model evaluated at the
#' #'   null estimate.
#' build_jacobian_at_null <- function(gamma0, slowing_models, shared) {
#'   if (slowing_models[[1]]$type != "proportional") {
#'     stop("Only proportional slowing is currently supported for the targeted test.")
#'   }
#'   if (slowing_models[[1]]$ref == "nc_spline") {
#'     times = slowing_models[[1]]$times
#'     knots = times[-c(1, length(times))]
#'     boundary_knots = c(times[1], times[length(times)])
#'   }
#'   if (shared) {
#'     jacobian_slowing_multiple_outcomes_shared(
#'       times = slowing_models %>% lapply(function(model) model$times),
#'       ref = slowing_models[[1]]$ref,
#'       params_list = gamma0,
#'       slowing_only = TRUE,
#'       knots = knots,
#'       boundary_knots = boundary_knots
#'     )
#'   } else {
#'     jacobian_slowing_multiple_outcomes(
#'       times = slowing_models %>% lapply(function(model) model$times),
#'       ref = slowing_models[[1]]$ref,
#'       params_list = gamma0,
#'       slowing_only = TRUE,
#'       knots = knots,
#'       boundary_knots = boundary_knots
#'     )
#'   }
#' }



#' Compute Chi-Squared Test Statistic
#'
#' @param B (matrix) The contrast matrix.
#' @inheritParams two_stage_gls_null
#' 
#' @returns (list) A list containing the following elements:
#' - `statistic`: The computed test statistic.
#' - `df`: The degrees of freedom for the test statistic.
#' - `p_value`: The p-value associated with the test statistic.
targeted_test_statistic <- function(B, m_tilde, Sigma) {
  BSB_inv <- solve(B %*% Sigma %*% Matrix::t(B))
  T_hat <- as.numeric(t(m_tilde) %*% Matrix::t(B) %*% BSB_inv %*% B %*% m_tilde)
  df <- nrow(B)
  p_value <- stats::pchisq(T_hat, df = df, lower.tail = FALSE)
  list(statistic = T_hat, df = df, p_value = p_value)
}