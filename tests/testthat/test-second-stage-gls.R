if (!file.exists(helper_file2)) {
  helper_file2 <- file.path("R", "two-stage-gls", "second-stage-gls.R")
  helper_file3 <- file.path("R", "local-power", "helper-functions-local-power.R")
  helper_file4 <- file.path("R", "two-stage-gls", "models.R")
}
source(helper_file2)
source(helper_file3)
source(helper_file4)
# ============================================================================
# Shared fixtures used across multiple tests
# ============================================================================

K  <- 4
J  <- 3
times_equal <- rep(list(0:K), J)

m_tilde <- rep(0:K, J * 2) 
Sigma <- build_kronecker_covariance(K = K, J = J)$Sigma

prop_slow_models_4PL <- make_slowing_models(ref = "4PL", times = times_equal, type = "proportional")
prop_slow_models_NC <- make_slowing_models(ref = "nc_spline", times = times_equal, type = "proportional")


# ============================================================================
# GLS Estimation
# ============================================================================

testthat::test_that("two_stage_gls_null() works for 4PL and NC spline", {
  gls_fitted_4PL <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    slowing_models = prop_slow_models,
    start = rep(1:2, J)
  )
  
  gls_fitted_NC <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    slowing_models = prop_slow_models_NC,
    start = rep(0:K, J)
  )
  
  # Check values for 4PL model
  
  value <- gls_fitted_4PL$optim$value
  param1 <- gls_fitted_4PL$gamma_hat[1]
  
  testthat::expect_equal(
    value, 7.411704234
  )
  testthat::expect_equal(
    param1, c("gamma0__outcome1_param1" = 0.9685776232)
  )
  
  # Check values for NC spline model
  
  value <- gls_fitted_NC$optim$value
  param1 <- gls_fitted_NC$gamma_hat[1]
  
  testthat::expect_equal(
    value, 0
  )
  testthat::expect_equal(
    param1, c("gamma0__outcome1_param1" = 1)
  )
  
})


