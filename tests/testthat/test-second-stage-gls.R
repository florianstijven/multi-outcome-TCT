helper_file2 <- file.path("..", "..", "R", "two-stage-gls", "second-stage-gls.R")
helper_file3 <- file.path("..", "..","R", "local-power", "helper-functions-local-power.R")
helper_file4 <- file.path("..", "..","R", "two-stage-gls", "models.R")
helper_file5 <- file.path("..", "..","R", "two-stage-gls", "slowing-models.R")


source(helper_file2)
source(helper_file3)
source(helper_file4)
source(helper_file5)
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

shared_quadratic_slow_model_4PL <- shared_parameter_model(
  make_slowing_models(ref = "4PL", times = times_equal, type = "quadratic"),
  list(c(3, 7, 11), c(4, 8, 12))
)

treatment_effect_params_positions <- list(
  c(K + 2, 2 * (K + 1) + 2 + 1, 3 * (K + 1) + 4 + 1),
  c(K + 3, 2 * (K + 1) + 2 + 2, 3 * (K + 1) + 4 + 2)
)
shared_quadratic_slow_model_NC <- shared_parameter_model(
  make_slowing_models(ref = "nc_spline", times = times_equal, type = "quadratic"),
  treatment_effect_params_positions
)


# ============================================================================
# GLS Estimation
# ============================================================================

testthat::test_that("two_stage_gls_null() works for 4PL and NC spline", {
  gls_fitted_4PL <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_4PL,
    start = rep(1:2, J)
  )
  
  gls_fitted_NC <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_NC,
    start = rep(0:K, J)
  )
  
  # Check values for 4PL model
  
  value <- gls_fitted_4PL$optim[[1]]$value
  param1 <- gls_fitted_4PL$gamma_hat[1]
  
  testthat::expect_equal(
    value, 7.411634625
  )
  testthat::expect_equal(
    param1, c(0.9798581652)
  )
  
  # Check values for NC spline model
  
  value <- gls_fitted_NC$optim[[1]]$value
  param1 <- gls_fitted_NC$gamma_hat[1]
  
  testthat::expect_equal(
    value, 0
  )
  testthat::expect_equal(
    param1, 1.000018456
  )
  
})


testthat::test_that("two_stage_gls_null() works for 4PL and NC spline with by outcome fitting", {
  gls_fitted_4PL <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_4PL,
    start = rep(1:2, J),
    split_indices_params = rep(1:J, each = 3),
    split_indices_mu = rep(1:J, each = 2 * (K + 1))
  )
  
  gls_fitted_NC <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_NC,
    start = rep(0:K, J),
    split_indices_params = rep(1:J, each = K + 2),
    split_indices_mu = rep(1:J, each = 2 * (K + 1))
  )
  
  gls_fitted_4PL_shared_quadratic <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = shared_quadratic_slow_model_4PL,
    start = rep(1:2, J),
    split_indices_params = c(as.list(rep(1:J, each = 2)), list(1:J), list(1:J)),
    split_indices_mu = rep(1:J, each = 2 * (K + 1))
  )
  
  gls_fitted_NC_shared_quadratic <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = shared_quadratic_slow_model_NC,
    start = rep(0:K, J),
    split_indices_params = c(as.list(rep(1:J, each = K + 1)), list(1:J), list(1:J)),
    split_indices_mu = rep(1:J, each = 2 * (K + 1))
  )
  
  # For the NC spline model, the minimized values should be equal to zero
  # because m_tilde satisfies the null exactly.
  testthat::expect_equal(
    purrr::map_dbl(gls_fitted_NC$optim, "value"),
    rep(0, J)
  )
  testthat::expect_equal(
    purrr::map_dbl(gls_fitted_NC_shared_quadratic$optim, "value"),
    rep(0, J)
  )
  
  # For the 4PL model, the minimized values should be equal some pre-computed
  # values.
  testthat::expect_equal(
    purrr::map_dbl(gls_fitted_4PL$optim, "value"),
    c(3.529349807, 6.176362163, 3.529349807)
  )
  testthat::expect_equal(
    purrr::map_dbl(gls_fitted_4PL_shared_quadratic$optim, "value"),
    c(3.529349807, 6.176362163, 3.529349807)
  )
})


testthat::test_that("two_stage_gls_null() works for 4PL and NC spline with OLS fitting", {
  gls_fitted_4PL <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_4PL,
    start = rep(1:2, J),
    ols = TRUE
  )
  
  gls_fitted_NC <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_NC,
    start = rep(0:K, J),
    ols = TRUE
  )
  
  # For the NC spline model, the minimized values should be equal to zero
  # because m_tilde satisfies the null exactly.
  testthat::expect_equal(
    purrr::map_dbl(gls_fitted_NC$optim, "value"),
    0
  )
  # For the 4PL model, the minimized values should be equal some pre-computed
  # values.
  testthat::expect_equal(
    purrr::map_dbl(gls_fitted_4PL$optim, "value"),
    c(84.00000203)
  )
})


