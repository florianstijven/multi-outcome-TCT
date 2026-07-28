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
  time_index <- 0:4
  slope <- 1.2
  inflection <- 2.0

  curve <- function_4PL(time_index, slope = slope, inflection = inflection)
  jac <- jacobian_4PL(time_index, slope = slope, inflection = inflection)
  d_time <- time_d_4PL(time_index, slope = slope, inflection = inflection)

  testthat::expect_length(curve, length(time_index))
  testthat::expect_true(all(curve > 0 & curve < 1))
  testthat::expect_equal(dim(jac), c(length(time_index), 2))
  testthat::expect_equal(d_time, jac[, 1], tolerance = 1e-12)
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


testthat::test_that("jacobian constructors return correct dimensions and slowing-only selections", {
  j <- 2
  k <- 3
  time_index <- 0:k
  params <- list(slope = c(1.1, 1.4), inflection = c(1.5, 2.2))

  jac_single <- jacobian_slowing_single_outcome(
    k = k,
    time_index = time_index,
    ref = "4PL",
    slope = params$slope[1],
    inflection = params$inflection[1]
  )

  testthat::expect_equal(dim(jac_single), c(2 * (k + 1), 3))
  testthat::expect_true(all(jac_single[1:(k + 1), 3] == 0))

  jac_multi <- jacobian_slowing_multiple_outcomes(
    j = j,
    k = k,
    time_index = time_index,
    ref = "4PL",
    params_list = params
  )
  testthat::expect_equal(dim(jac_multi), c(j * 2 * (k + 1), 3 * j))

  jac_multi_slow <- jacobian_slowing_multiple_outcomes(
    j = j,
    k = k,
    time_index = time_index,
    ref = "4PL",
    params_list = params,
    slowing_only = TRUE
  )
  testthat::expect_equal(dim(jac_multi_slow), c(j * 2 * (k + 1), j))
  testthat::expect_equal(jac_multi_slow, jac_multi[, seq(3, 3 * j, by = 3), drop = FALSE])

  jac_shared_slow <- jacobian_slowing_multiple_outcomes_shared(
    j = j,
    k = k,
    time_index = time_index,
    ref = "4PL",
    params_list = params,
    slowing_only = TRUE
  )
  testthat::expect_equal(dim(jac_shared_slow), c(j * 2 * (k + 1), 1))
})


testthat::test_that("mean and shift vectors are well-formed", {
  j <- 3
  k <- 2
  time_index <- rep(0:k, 2)
  params <- list(slope = c(1.0, 1.2, 1.5), inflection = c(1.0, 1.5, 2.0))

  mu <- mean_vector(
    ref = "4PL",
    j = j,
    k = k,
    time_index = time_index,
    params_list = params
  )
  testthat::expect_length(mu, j * length(time_index))
  testthat::expect_true(all(mu > 0 & mu < 1))

  shift <- local_shift_vector(
    ref = "4PL",
    j = j,
    k = k,
    time_index = 0:k,
    h = c(1, 2, 3),
    params_list = params
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


testthat::test_that("local_shift_vector() and mean_vector() return correct values", {
  j <- 3
  k <- 3
  time_index <- 0:k
  ref <- "4PL"
  params_list = list(inflection = c(1.0, 1.5, 3), slope = c(1.0, 1.2, 1.5))
  
  result1 <- mean_vector(ref = ref, j = j, k = k, time_index = rep(time_index, 2), params_list = params_list)
  result2 <- local_shift_vector(ref = ref, j = j, k = k, time_index = time_index, h = c(1, 2, 3), params_list = params_list)
  
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
})
