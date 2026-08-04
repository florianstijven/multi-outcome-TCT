testthat::test_that("power_from_ncp returns type-I level at ncp=0 and increases with ncp", {
  alpha <- 0.05
  p0 <- power_from_ncp(ncp = 0, df = 3, alpha = alpha)
  p1 <- power_from_ncp(ncp = 2, df = 3, alpha = alpha)
  p2 <- power_from_ncp(ncp = 8, df = 3, alpha = alpha)

  testthat::expect_equal(p0, alpha, tolerance = 1e-10)
  testthat::expect_true(p0 < p1)
  testthat::expect_true(p1 < p2)
})


testthat::test_that("4PL reference helpers return expected shapes and ranges", {
  times <- 0:4
  slope <- 1.2
  inflection <- 2.0
  params <- c(inflection, slope)

  curve <- function_4PL(times, params)
  jac <- jacobian_4PL(times, params)
  d_time <- time_d_4PL(times, params)

  testthat::expect_length(curve, length(times))
  testthat::expect_true(all(curve > 0 & curve < 1))
  testthat::expect_equal(dim(jac), c(length(times), 2))
  testthat::expect_equal(d_time, jac[, 1], tolerance = 1e-12)
})

testthat::test_that("nc spline reference helpers return expected shapes and ranges", {
  times <- 0:8
  params <- 1:5
  boundary_knots <- range(times)

  t <- build_nc_spline_basis(times = times, knots = 1:3, boundary_knots = boundary_knots, params = params)
  curve <- function_nc_spline(times, knots = 1:3, boundary_knots = boundary_knots, params = params)
  jac <- jacobian_nc_spline(times, knots = 1:3, boundary_knots = boundary_knots, params = params)
  d_time <- time_d_nc_spline(times, knots = 1:3, boundary_knots = boundary_knots, params = params)

  testthat::expect_length(curve, length(times))
  testthat::expect_equal(dim(jac), c(length(times), length(params)))
})

testthat::test_that("nc spline helpers enforce boundary_knots input", {
  times <- 0:8
  params <- 1:5
  internal_knots <- 1:3

  testthat::expect_error(
    build_nc_spline_basis(times = times, knots = internal_knots, params = params),
    "argument \"boundary_knots\" is missing"
  )
  testthat::expect_error(
    function_nc_spline(times = times, knots = internal_knots, params = params),
    "argument \"boundary_knots\" is missing"
  )
  testthat::expect_error(
    jacobian_nc_spline(times = times, knots = internal_knots, params = params),
    "argument \"boundary_knots\" is missing"
  )
  testthat::expect_error(
    time_d_nc_spline(times = times, knots = internal_knots, params = params),
    "argument \"boundary_knots\" is missing"
  )

  bad_boundary_knots <- c(0, 4, 8)
  testthat::expect_error(
    build_nc_spline_basis(
      times = times,
      knots = internal_knots,
      boundary_knots = bad_boundary_knots,
      params = params
    ),
    "boundary knots must be a numeric vector of length 2"
  )
  testthat::expect_error(
    function_nc_spline(
      times = times,
      knots = internal_knots,
      boundary_knots = bad_boundary_knots,
      params = params
    ),
    "boundary knots must be a numeric vector of length 2"
  )
  testthat::expect_error(
    jacobian_nc_spline(
      times = times,
      knots = internal_knots,
      boundary_knots = bad_boundary_knots,
      params = params
    ),
    "boundary knots must be a numeric vector of length 2"
  )
  testthat::expect_error(
    time_d_nc_spline(
      times = times,
      knots = internal_knots,
      boundary_knots = bad_boundary_knots,
      params = params
    ),
    "boundary knots must be a numeric vector of length 2"
  )
})




testthat::test_that("omnibus and summing contrast matrices have expected dimensions and entries", {
  k <- 3
  single <- build_omnibus_contrast_single_outcome(k)

  expected_single <- matrix(c(
    0, -1, 0, 0, 0, 1, 0, 0,
    0, 0, -1, 0, 0, 0, 1, 0,
    0, 0, 0, -1, 0, 0, 0, 1
  ), nrow = k, byrow = TRUE)

  testthat::expect_equal(single, expected_single)

  j <- 2
  multi <- build_omnibus_contrast_multi_outcome(j, k)
  testthat::expect_equal(dim(multi), c(j * k, j * 2 * (k + 1)))

  sum_contrast <- build_summing_contrast_multi_outcome(j, k)
  testthat::expect_equal(dim(sum_contrast), c(k, j * 2 * (k + 1)))

  expected_sum <- cbind(expected_single, expected_single)
  testthat::expect_equal(sum_contrast, expected_sum)
})

testthat::test_that("build_omnibus_contrast_multi_outcome() and build_summing_contrast_multi_outcome() work with list-valued argument", {
  K <- 4
  J <- 3
  times = list(0:K, 0:K, 0:K)
  
  
  testthat::expect_equal(
    build_omnibus_contrast_multi_outcome(J, K),
    build_omnibus_contrast_multi_outcome(times = times)
  )
  testthat::expect_equal(
    build_summing_contrast_multi_outcome(J, K),
    build_summing_contrast_multi_outcome(times = times)
  )
})

testthat::test_that("build_linear_contrast_multi_outcome() and build_linear_contrast_common() return correct values", {
  times <- list(0:3, 0:3, 0:5)
  Sigma <- diag(1, sum(sapply(times, length)) * 2)

  testthat::expect_equal(
    build_linear_contrast_multi_outcome(times = times, Sigma = Sigma)[c(2, 3), c(16, 27)],
    matrix(c(-1.5, 0, 0, -2), nrow = 2, byrow = TRUE)
  )
  
  testthat::expect_equal(
    build_linear_contrast_common(times = times, Sigma = Sigma)[1, c(16, 27)],
    c(-1.5, -2), nrow = 1
  )
})

testthat::test_that("jacobian constructors return correct dimensions and slowing-only selections", {
  j <- 2
  k <- 3
  times <- 0:k
  params_list <- list(
    c(1.5, 1.1),
    c(2.2, 1.4)
  )

  jac_single <- jacobian_slowing_single_outcome(
    K = k,
    times = times,
    ref = "4PL",
    params = params_list[[1]]
  )

  testthat::expect_equal(dim(jac_single), c(2 * (k + 1), 3))
  testthat::expect_true(all(jac_single[1:(k + 1), 3] == 0))

  jac_multi <- jacobian_slowing_multiple_outcomes(
    J = j,
    K = k,
    times = times,
    ref = "4PL",
    params_list = params_list
  )
  testthat::expect_equal(dim(jac_multi), c(j * 2 * (k + 1), 3 * j))

  jac_multi_slow <- jacobian_slowing_multiple_outcomes(
    J = j,
    K = k,
    times = times,
    ref = "4PL",
    params_list = params_list,
    slowing_only = TRUE
  )
  testthat::expect_equal(dim(jac_multi_slow), c(j * 2 * (k + 1), j))
  testthat::expect_equal(jac_multi_slow, jac_multi[, seq(3, 3 * j, by = 3), drop = FALSE])

  jac_shared_slow <- jacobian_slowing_multiple_outcomes_shared(
    J = j,
    K = k,
    times = times,
    ref = "4PL",
    params_list = params_list,
    slowing_only = TRUE
  )
  testthat::expect_equal(dim(jac_shared_slow), c(j * 2 * (k + 1), 1))
})


testthat::test_that("mean and shift vectors are well-formed", {
  j <- 3
  k <- 2
  times <- rep(0:k, 2)
  params_list <- list(
    c(1, 1),
    c(1.5, 1.2),
    c(2, 1.5)
  )

  mu <- mean_vector(
    ref = "4PL",
    J = j,
    K = k,
    times = times,
    params_list = params_list
  )
  testthat::expect_length(mu, j * length(times))
  testthat::expect_true(all(mu > 0 & mu < 1))

  shift <- local_shift_vector_slowing_outcome(
    ref = "4PL",
    J = j,
    K = k,
    times = 0:k,
    h = c(1, 2, 3),
    params_list = params_list
  )
  testthat::expect_equal(dim(shift), c(j * 2 * (k + 1), 1))
})


testthat::test_that("covariance builder returns symmetric covariance with expected dimensions", {
  K <- 4
  J <- 3
  out <- build_kronecker_covariance(K = K, J = J, rho_time = 0.7, rho_outcome = 0.6)

  testthat::expect_named(out, c("Sigma", "time_points"))
  testthat::expect_equal(length(out$time_points), K + 1)
  testthat::expect_equal(dim(out$Sigma), c(2 * J * (K + 1), 2 * J * (K + 1)))
  testthat::expect_equal(out$Sigma, t(out$Sigma), tolerance = 1e-12)
  testthat::expect_true(all(diag(out$Sigma) > 0))
})


testthat::test_that("compute_power_curve returns one row per h and monotone power for simple identity setup", {
  B <- matrix(c(1, 0, 0, 1), nrow = 2)
  Sigma <- diag(2)
  base_effect <- c(1, 2)
  h_grid <- c(0, 0.5, 1.0)

  out <- compute_power_curve(
    B = B,
    h_grid = h_grid,
    alpha = 0.05,
    Sigma = Sigma,
    base_effect = base_effect
  )

  testthat::expect_equal(nrow(out), length(h_grid))
  testthat::expect_equal(out$h, h_grid)
  testthat::expect_equal(out$df, rep(nrow(B), length(h_grid)))
  testthat::expect_true(all(diff(out$ncp) > 0))
  testthat::expect_true(all(diff(out$power) > 0))
})


testthat::test_that("local_shift_vector_slowing_outcome() and mean_vector() return correct values", {
  j <- 3
  k <- 3
  times <- 0:k
  ref <- "4PL"
  params_list <- list(
    c(1, 1),
    c(1.5, 1.2),
    c(3, 1.5)
  )
  
  result1 <- mean_vector(
    ref = ref,
    J = j,
    K = k,
    times = rep(times, 2),
    params_list = params_list
  )
  result2 <- local_shift_vector_slowing_outcome(
    ref = ref,
    J = j,
    K = k,
    times = times,
    h = c(1, 2, 3),
    params_list = params_list
  )
  
  expected1 <- c(
    0.26894142137,
    0.50000000000,
    0.73105857863,
    0.88079707798,
    0.26894142137,
    0.50000000000,
    0.73105857863,
    0.88079707798,
    0.14185106490,
    0.35434369377,
    0.64565630623,
    0.85814893510,
    0.14185106490,
    0.35434369377,
    0.64565630623,
    0.85814893510,
    0.01098694263,
    0.04742587318,
    0.18242552381,
    0.50000000000,
    0.01098694263,
    0.04742587318,
    0.18242552381,
    0.50000000000
  )
  expected2 <- c(
    0,
    0,
    0,
    0,
    0,
    -0.25,
    -0.3932238665,
    -0.3149807562,
    0,
    0,
    0,
    0,
    0,
    -0.5490821771,
    -1.0981643542,
    -0.8764512501,
    0,
    0,
    0,
    0,
    0,
    -0.2032949688,
    -1.3423180686,
    -3.3750000000
  )

  testthat::expect_equal(result1, expected1, tolerance = 1e-10)
  testthat::expect_equal(as.numeric(result2), expected2, tolerance = 1e-10)

  ref <- "nc_spline"
  params_list <- list(
    c(0.5, 1.0, 1.5, 2),
    c(0.5, 1.0, 1.5, 2),
    c(0.5, 1.0, 1.5, 2)
  )

  result3 <- mean_vector(
    ref = ref,
    J = j,
    K = k,
    times = rep(times, 2),
    params_list = params_list,
    knots = c(1.5, 2.5),
    boundary_knots = range(times)
  )
  result4 <- local_shift_vector_slowing_outcome(
    ref = ref,
    J = j,
    K = k,
    times = times,
    h = c(1, 2, 3),
    params_list = params_list,
    knots = c(1.5, 2.5),
    boundary_knots = range(times)
  )
})
# 
# testthat::test_that("jacobian constructors support different parameter counts per outcome", {
#   times <- list(0:3, 0:4)
#   j <- length(times)
#   params <- list(
#     list(knots = c(1.5), coeffs = c(0.2, 0.4, 0.6)),
#     list(knots = c(1.5, 2.5), coeffs = c(0.2, 0.4, 0.6, 0.8))
#   )
# 
#   jac_multi <- jacobian_slowing_multiple_outcomes(
#     J = j,
#     K = NULL,
#     times = times,
#     ref = "nc_spline",
#     params_list = params
#   )
# 
#   # outcome 1: 3 spline coefficients + 1 slowing; outcome 2: 4 coefficients + 1 slowing
#   testthat::expect_equal(ncol(jac_multi), 9)
#   testthat::expect_equal(nrow(jac_multi), 2 * (length(times[[1]]) + length(times[[2]])))
# 
#   jac_multi_slow <- jacobian_slowing_multiple_outcomes(
#     J = j,
#     K = NULL,
#     times = times,
#     ref = "nc_spline",
#     params_list = params,
#     slowing_only = TRUE
#   )
#   testthat::expect_equal(dim(jac_multi_slow), c(nrow(jac_multi), j))
# 
#   jac_shared <- jacobian_slowing_multiple_outcomes_shared(
#     J = j,
#     K = NULL,
#     times = times,
#     ref = "nc_spline",
#     params_list = params,
#     slowing_only = TRUE
#   )
#   testthat::expect_equal(dim(jac_shared), c(nrow(jac_multi), 1))
# })


