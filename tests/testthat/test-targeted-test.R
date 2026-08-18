helper_file2 <- file.path("..", "..", "R", "two-stage-gls", "second-stage-gls.R")
helper_file3 <- file.path("..", "..","R", "local-power", "helper-functions-local-power.R")
helper_file4 <- file.path("..", "..","R", "two-stage-gls", "models.R")
helper_file5 <- file.path("..", "..","R", "two-stage-gls", "slowing-models.R")
helper_file6 <- file.path("..", "..","R", "targeted-test", "targeted-test.R")


source(helper_file2)
source(helper_file3)
source(helper_file4)
source(helper_file5)
source(helper_file6)
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

quadratic_slow_models_4PL <- make_slowing_models(ref = "4PL", times = times_equal, type = "quadratic")
quadratic_slow_models_NC <- make_slowing_models(ref = "nc_spline", times = times_equal, type = "quadratic")

shared_slowing_model_4PL <- shared_parameter_model(prop_slow_models_4PL, c(3, 6, 9))
shared_quadratic_slow_model_4PL <- shared_parameter_model(quadratic_slow_models_4PL, list(c(3, 7, 11), c(4, 8, 12)))

# ============================================================================
# GLS Estimation
# ============================================================================

testthat::test_that("targeted_test() returns expected value when the null holds exactly for proportional slowing models", {
  A <- build_omnibus_contrast_multi_outcome(J, K)
  targeted_test_4PL <- targeted_test(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_4PL,
    A = A,
    start = rep(1:2, J),
    shared_matrix = NULL
  )
  targeted_test_NC <- targeted_test(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = prop_slow_models_NC,
    A = A,
    start = rep(0:K, J),
    shared_matrix = NULL
  )
  targeted_test_shared_4PL <- targeted_test(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = shared_slowing_model_4PL,
    A = A,
    start = rep(1:2, J),
    shared_matrix = NULL
  )
  
  # Check values for 4PL model
  
  p_value_4PL <- targeted_test_4PL$p_value
  testthat::expect_equal(p_value_4PL, 1)
  testthat::expect_equal(targeted_test_4PL$df, J)
  
  p_value_NC <- targeted_test_NC$p_value
  testthat::expect_equal(p_value_NC, 1)
  testthat::expect_equal(targeted_test_NC$df, J)
  
  p_value_shared_4PL <- targeted_test_shared_4PL$p_value
  testthat::expect_equal(p_value_shared_4PL, 1)
  testthat::expect_equal(targeted_test_shared_4PL$df, 1)
})

testthat::test_that("targeted_test() returns expected value when the null holds exactly for quadratic slowing models", {
  A <- build_omnibus_contrast_multi_outcome(J, K)
  targeted_test_4PL <- targeted_test(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = quadratic_slow_models_4PL,
    A = A,
    start = rep(1:2, J),
    shared_matrix = NULL
  )
  targeted_test_NC <- targeted_test(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = quadratic_slow_models_NC,
    A = A,
    start = rep(0:K, J),
    shared_matrix = NULL
  )
  targeted_test_shared_4PL <- targeted_test(
    m_tilde = m_tilde,
    Sigma   = Sigma,
    working_model = shared_quadratic_slow_model_4PL,
    A = A,
    start = rep(1:2, J),
    shared_matrix = NULL
  )
  
  # Check values for 4PL model
  
  p_value_4PL <- targeted_test_4PL$p_value
  testthat::expect_equal(p_value_4PL, 1)
  testthat::expect_equal(targeted_test_4PL$df,  J * 2)
  
  p_value_NC <- targeted_test_NC$p_value
  testthat::expect_equal(p_value_NC, 1)
  testthat::expect_equal(targeted_test_NC$df,  J * 2)
  
  p_value_shared_4PL <- targeted_test_shared_4PL$p_value
  testthat::expect_equal(p_value_shared_4PL, 1)
  testthat::expect_equal(targeted_test_shared_4PL$df,  2)
})

testthat::test_that("targeted_test() returns expected value when the null does not hold for m_tilde for proportional slowing", {
  m_tilde_alt <- m_tilde + rep(c(rep(0, 6), rep(1, 4)), J)
  A <- build_omnibus_contrast_multi_outcome(J, K)
  targeted_test_4PL <- targeted_test(
    m_tilde = m_tilde_alt,
    Sigma   = Sigma,
    working_model = prop_slow_models_4PL,
    A = A,
    start = rep(1:2, J),
    shared_matrix = NULL
  )
  targeted_test_NC <- targeted_test(
    m_tilde = m_tilde_alt,
    Sigma   = Sigma,
    working_model = prop_slow_models_NC,
    A = A,
    start = rep(0:K, J),
    shared_matrix = NULL
  )
  targeted_test_shared_4PL <- targeted_test(
    m_tilde = m_tilde_alt,
    Sigma   = Sigma,
    working_model = shared_slowing_model_4PL,
    A = A,
    start = rep(1:2, J),
    shared_matrix = NULL
  )
  
  # Check values for 4PL model
  
  p_value_4PL <- targeted_test_4PL$p_value
  testthat::expect_equal(p_value_4PL, 0.8457513252)
  
  p_value_NC <- targeted_test_NC$p_value
  testthat::expect_equal(p_value_NC, 0.8767681933)
  
  p_value_shared_4PL <- targeted_test_shared_4PL$p_value
  testthat::expect_equal(p_value_shared_4PL, 0.3583682035)
})

testthat::test_that("targeted_test() returns expected value when the null does not hold for m_tilde for quadratic slowing", {
  m_tilde_alt <- m_tilde + rep(c(rep(0, 5), 0:4), J)
  A <- build_omnibus_contrast_multi_outcome(J, K)
  targeted_test_4PL <- targeted_test(
    m_tilde = m_tilde_alt,
    Sigma   = Sigma,
    working_model = quadratic_slow_models_4PL,
    A = A,
    start = rep(1:2, J),
    shared_matrix = NULL
  )
  targeted_test_NC <- targeted_test(
    m_tilde = m_tilde_alt,
    Sigma   = Sigma,
    working_model = quadratic_slow_models_NC,
    A = A,
    start = rep(0:K, J),
    shared_matrix = NULL
  )
  targeted_test_shared_4PL <- targeted_test(
    m_tilde = m_tilde_alt,
    Sigma   = Sigma,
    working_model = shared_quadratic_slow_model_4PL,
    A = A,
    start = rep(1:2, J),
    shared_matrix = NULL
  )
  
  # Check values for 4PL model
  
  p_value_4PL <- targeted_test_4PL$p_value
  testthat::expect_equal(p_value_4PL, 0.9309020851)
  
  p_value_NC <- targeted_test_NC$p_value
  testthat::expect_equal(p_value_NC, 0.002891371836)
  
  p_value_shared_4PL <- targeted_test_shared_4PL$p_value
  testthat::expect_equal(p_value_shared_4PL, 0.1318756374)
})

