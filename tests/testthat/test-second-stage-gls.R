# ============================================================================
# Shared fixtures used across multiple tests
# ============================================================================

K  <- 4
J  <- 3
times_equal <- rep(list(0:K), J)

m_tilde <- rep(0:K, J * 2) 
Sigma <- build_kronecker_covariance(K = K, J = J)$Sigma

prop_slow_models_4PL <- make_slowing_models(ref = "4PL", times = times_equal, type = "proportional")
prop_slow_models_4PL_asympt <- make_slowing_models(ref = "4PL_asympt", times = times_equal, type = "proportional")
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
# GLS Estimation under the null
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
  param1 <- coef(gls_fitted_4PL)[1]
  
  testthat::expect_equal(
    value, 7.411634625
  )
  testthat::expect_equal(
    param1, c("ref_param_1 (Outcome 1)" = 0.9798581652)
  )
  
  # Check values for NC spline model
  
  value <- gls_fitted_NC$optim[[1]]$value
  param1 <- coef(gls_fitted_NC)[1]
  
  testthat::expect_equal(
    value, 0
  )
  testthat::expect_equal(
    param1, c("ref_param_1 (Outcome 1)" = 1.000018456)
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

testthat::test_that("two_stage_gls_null() works for 4PL with asymptotes as parameters", {
  gls_fitted_4PL_asympt <- two_stage_gls_null(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_4PL_asympt,
    start = rep(c(1:2, 5, 0), J),
    ols = TRUE
  )
  
  # For the 4PL model, the minimized values should be equal some pre-computed
  # values.
  testthat::expect_equal(
    purrr::map_dbl(gls_fitted_4PL_asympt$optim, "value"),
    c(0.8909656715)
  )
})


# ============================================================================
# GLS Estimation of the full model
# ============================================================================

testthat::test_that("two_stage_gls_full() works for 4PL and NC spline", {
  gls_fitted_4PL <- two_stage_gls_full(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_4PL,
    start = rep(c(1:2, 1), J)
  )
  
  gls_fitted_NC <- two_stage_gls_full(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_NC,
    start = rep(c(0:K, 1), J)
  )
  
  gls_fitted_4PL_shared_quadratic <- two_stage_gls_full(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = shared_quadratic_slow_model_4PL,
    start = c(rep(1:2, J), 1, 0)
  )
  
  gls_fitted_NC_shared_quadratic <- two_stage_gls_full(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = shared_quadratic_slow_model_NC,
    start = c(rep(0:K, J), 1, 0)
  )
  
  # Check values for 4PL model
  
  value <- gls_fitted_4PL$optim[[1]]$value
  param1 <- coef(gls_fitted_4PL)[1]
  
  testthat::expect_equal(
    value, 7.411634625
  )
  testthat::expect_equal(
    param1, c("ref_param_1 (Outcome 1)" = 0.9796054359)
  )
  
  # Check values for NC spline model
  
  value <- gls_fitted_NC$optim[[1]]$value
  param1 <- coef(gls_fitted_NC)[1]
  
  testthat::expect_equal(
    value, 0
  )
  testthat::expect_equal(
    param1, c("ref_param_1 (Outcome 1)" = 1.000005415)
  )
  
  # Check values for the quadratic models
  value3 <- gls_fitted_4PL_shared_quadratic$optim[[1]]$value
  value4 <- gls_fitted_NC_shared_quadratic$optim[[1]]$value
  param3 <- coef(gls_fitted_4PL_shared_quadratic)[1]
  param4 <- coef(gls_fitted_NC_shared_quadratic)[1]
  
  testthat::expect_equal(
    value3, 7.411634582
  )
  testthat::expect_equal(
    value4, 1.955419344e-07
  )
  testthat::expect_equal(
    param3, c("ref_param_1 (Outcome 1)" = 0.9796109976)
  )
  testthat::expect_equal(
    param4, c("ref_param_1 (Outcome 1)" = 0.9999877145)
  )
})


testthat::test_that("two_stage_gls_full() works for 4PL and NC spline with by outcome fitting", {
  gls_fitted_4PL <- two_stage_gls_full(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_4PL,
    start = rep(c(1:2, 1), J),
    split_indices_params = rep(1:J, each = 3),
    split_indices_mu = rep(1:J, each = 2 * (K + 1))
  )
  
  gls_fitted_NC <- two_stage_gls_full(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_NC,
    start = rep(c(0:K, 1), J),
    split_indices_params = rep(1:J, each = K + 2),
    split_indices_mu = rep(1:J, each = 2 * (K + 1))
  )
  

  
  # For the NC spline model, the minimized values should be equal to zero
  # because m_tilde satisfies the null exactly.
  testthat::expect_equal(
    purrr::map_dbl(gls_fitted_NC$optim, "value"),
    rep(0, J)
  )
  
  # For the 4PL model, the minimized values should be equal some pre-computed
  # values.
  testthat::expect_equal(
    purrr::map_dbl(gls_fitted_4PL$optim, "value"),
    c(3.529349807, 6.176362163, 3.529349807)
  )
})







