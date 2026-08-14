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

targeted_test_proportional_slowing <- function(m_tilde, Sigma, shared, slowing_models, start) {
  # Step 1: Fit working model under the null.
  gls_fitted_null <- two_stage_gls_null(m_tilde = m_tilde,
                                        Sigma   = Sigma,
                                        slowing_models = slowing_models, 
                                        start = start)
  
  # Step 2: build the Jacobian at the null estimate
  
  jac <- build_jacobian_at_null(
    gamma0  = gls_fitted_null$gamma_hat_list$gamma0_hat,
    slowing_models = slowing_models,
    shared      = shared
  )
  
  
  # Step 3: build the omnibus contrast matrix and then B_n
  A <- build_omnibus_contrast_multi_outcome(
    times = slowing_models %>% lapply(function(model) model$times)
  )
  B_n <- build_contrast_matrix(jacobian = jac,
                               A = A,
                               Sigma = Sigma)
  
  # Step 4: evaluate the test statistic
  test <- targeted_test_statistic(B_n = B_n,
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


# ============================================================================
# Plug-in contrast matrix B_n
# ============================================================================

# Compute B_n = dot_mu(gamma_hat)^T A^T (A Sigma A^T)^{-1} A
# using only the treatment-effect columns of the Jacobian (Remark in manuscript:
# reference-trajectory columns are orthogonal to A and can be dropped).
#
# jacobian    : full Jacobian matrix from build_jacobian_at_null()
# A           : omnibus contrast matrix
# Sigma       : first-stage covariance estimate
# shared      : if TRUE, the last column of the Jacobian is the shared slowing
#               parameter; if FALSE, the last J columns are outcome-specific
# J           : number of outcomes (only needed when shared = FALSE)

build_Bn <- function(jacobian, A, Sigma, shared = TRUE, J = NULL) {
  # Extract only the treatment-effect columns: last 1 (shared) or last J
  # (outcome-specific) columns. These are the D_time columns from Remark
  # in the manuscript; reference-trajectory columns satisfy A D_1 = 0.
  if (shared) {
    jac_time <- jacobian[, ncol(jacobian), drop = FALSE]
  } else {
    if (is.null(J)) stop("J must be supplied when shared = FALSE.")
    jac_time <- jacobian[, (ncol(jacobian) - J + 1):ncol(jacobian), drop = FALSE]
  }
  
  build_contrast_matrix(jacobian = jac_time, A = A, Sigma = Sigma)
}

build_jacobian_at_null <- function(gamma0, slowing_models, shared) {
  if (slowing_models[[1]]$type != "proportional") {
    stop("Only proportional slowing is currently supported for the targeted test.")
  }
  if (slowing_models[[1]]$ref == "nc_spline") {
    times = slowing_models[[1]]$times
    knots = times[-c(1, length(times))]
    boundary_knots = c(times[1], times[length(times)])
  }
  if (shared) {
    jacobian_slowing_multiple_outcomes_shared(
      times = slowing_models %>% lapply(function(model) model$times),
      ref = slowing_models[[1]]$ref,
      params_list = gamma0,
      slowing_only = TRUE,
      knots = knots,
      boundary_knots = boundary_knots
    )
  } else {
    jacobian_slowing_multiple_outcomes(
      times = slowing_models %>% lapply(function(model) model$times),
      ref = slowing_models[[1]]$ref,
      params_list = gamma0,
      slowing_only = TRUE,
      knots = knots,
      boundary_knots = boundary_knots
    )
  }
}


# ============================================================================
# Targeted test statistic
# ============================================================================

# Compute T_hat_n = n * m_tilde^T B_n^T (B_n Sigma B_n^T)^{-1} B_n m_tilde
# and return the chi-squared p-value with df = nrow(B_n).
#
# B_n     : contrast matrix from build_Bn()
# m_tilde : first-stage stacked mean estimate
# Sigma   : first-stage covariance estimate
# n       : total sample size (used to scale the statistic)

targeted_test_statistic <- function(B_n, m_tilde, Sigma) {
  BSB_inv <- solve(B_n %*% Sigma %*% t(B_n))
  T_hat <- as.numeric(t(m_tilde) %*% t(B_n) %*% BSB_inv %*% B_n %*% m_tilde)
  df <- nrow(B_n)
  p_value <- stats::pchisq(T_hat, df = df, lower.tail = FALSE)
  list(statistic = T_hat, df = df, p_value = p_value)
}