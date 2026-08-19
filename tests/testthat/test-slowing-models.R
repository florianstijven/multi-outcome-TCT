K  <- 4
J  <- 3
times_equal <- rep(list(0:K), J)

m_tilde <- rep(0:K, J * 2)
Sigma <- build_kronecker_covariance(K = K, J = J)$Sigma

prop_slow_models_4PL <- make_slowing_models(ref = "4PL", times = times_equal, type = "proportional")
prop_slow_models_NC <- make_slowing_models(ref = "nc_spline", times = times_equal, type = "proportional")


ref = "4PL"
type = "proportional"

# Extract from make_slowing_models() to test numeric Jacobian.
time_mapping = function(t, gamma1) {
  gamma1 * t
}
null_gamma1 = 1

time_mapping_deriv = function(t, gamma1) {
  t
}

reference_trajectory_functions_list = reference_trajectory_f_list(times, ref)

times = times_equal[[1]]

prop_slow_model_4PL_numerical <- model(
  mean_fn = function(gamma) {
    # Number of reference-trajectory parameters.
    no_params_ref = reference_trajectory_functions_list$no_params
    
    gamma0 = gamma[1:no_params_ref]
    gamma1 = gamma[(no_params_ref + 1):length(gamma)]
    
    control_means = reference_trajectory_functions_list$eval(times, gamma0)
    exp_means = reference_trajectory_functions_list$eval(time_mapping(times, gamma1), gamma0)
    
    c(control_means, exp_means)
  },
  nuisance_params_position = seq_len(reference_trajectory_functions_list$no_params),
  treatment_params_null = null_gamma1
)

prop_slow_models_4PL_numerical <- concatenate_models(rep(list(prop_slow_model_4PL_numerical), 3))

testthat::test_that("Analytical and numerical Jacobians match for selected slowing models",
                    {
                      # Outcome-Specific Slowing
                      jac_num1 <- prop_slow_models_4PL_numerical$jacobian_fn(rep(1, 9))
                      jac_an1 <- prop_slow_models_4PL$jacobian_fn(rep(1, 9))
                      
                      testthat::expect_equal(jac_num, jac_an, tolerance = 1e-10)
                      
                      # Shared Slowing
                      jac_num2 <- shared_parameter_model(prop_slow_models_4PL_numerical, c(3, 6, 9))$jacobian_fn(rep(1, 7))
                      jac_an2 <- shared_parameter_model(prop_slow_models_4PL, c(3, 6, 9))$jacobian_fn(rep(1, 7))
                      
                      testthat::expect_equal(jac_num2, jac_an2, tolerance = 1e-10)
                    })
