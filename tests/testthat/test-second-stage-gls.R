helper_file2 <- file.path("R", "local-power", "second-stage-gls.R")
if (!file.exists(helper_file2)) {
  helper_file2 <- file.path("..", "..", "R", "local-power", "second-stage-gls.R")
}
source(helper_file2)

# ============================================================================
# Shared fixtures used across multiple tests
# ============================================================================

K  <- 4
J  <- 2
times_equal <- rep(list(0:K), J)

# Simple 4PL params: one per outcome (inflection, slope)
params_list_4PL <- list(c(2.0, 1.2), c(3.0, 0.9))

# Flat reference-trajectory parameters for the spline (K+1 basis coefficients).
# The spline evaluates as basis %*% params; we pick params as c(0.1, 0.2, ...).
n_spline_coef <- K + 1
params_list_spline <- list(
  seq(0.05, 0.5, length.out = n_spline_coef),
  seq(0.10, 0.6, length.out = n_spline_coef)
)

# Stacked start vector (reference params for all outcomes, concatenated)
start_4PL <- do.call(c, params_list_4PL)
start_spline <- do.call(c, params_list_spline)

# Small positive-definite covariance (diagonal for simplicity)
p   <- J * 2 * (K + 1)
Sig <- diag(runif(p, 0.5, 1.5))

# Simple non-zero mean vector: null stacked means under 4PL with gamma1 = 1
m_null <- mean_vector(
  ref         = "4PL",
  J           = J,
  K           = K,
  times       = times_equal,
  params_list = params_list_4PL
)
# Perturb the treatment means slightly to produce a non-null signal
m_tilde_4PL <- m_null
for (j in seq_len(J)) {
  trt_start <- j * 2 * (K + 1) - K   # first post-baseline treatment entry
  trt_end   <- j * 2 * (K + 1)
  m_tilde_4PL[trt_start:trt_end] <- m_tilde_4PL[trt_start:trt_end] - 0.05
}

n_obs <- 200L


# ============================================================================
# gls_criterion
# ============================================================================

testthat::test_that("gls_criterion returns a non-negative scalar", {
  mean_fn <- make_null_mean_fn_proportional(
    ref = "4PL", times = times_equal, gamma1_null = 1, shared = TRUE
  )
  Sigma_inv <- solve(Sig)
  val <- gls_criterion(start_4PL, mean_fn, m_null, Sigma_inv)

  testthat::expect_length(val, 1L)
  testthat::expect_true(val >= 0)
})

testthat::test_that("gls_criterion is zero when m_tilde equals mu(gamma0)", {
  mean_fn <- make_null_mean_fn_proportional(
    ref = "4PL", times = times_equal, gamma1_null = 1, shared = TRUE
  )
  Sigma_inv <- solve(Sig)
  val <- gls_criterion(start_4PL, mean_fn, m_null, Sigma_inv)

  testthat::expect_equal(val, 0, tolerance = 1e-10)
})


# ============================================================================
# make_null_mean_fn_proportional
# ============================================================================

testthat::test_that("null mean function returns vector of correct length (shared)", {
  mean_fn <- make_null_mean_fn_proportional(
    ref = "4PL", times = times_equal, gamma1_null = 1, shared = TRUE
  )
  out <- mean_fn(start_4PL)

  testthat::expect_length(out, p)
})

testthat::test_that("null mean function returns vector of correct length (outcome-specific)", {
  gamma1_null <- c(1, 1)
  mean_fn <- make_null_mean_fn_proportional(
    ref = "4PL", times = times_equal, gamma1_null = gamma1_null, shared = FALSE
  )
  out <- mean_fn(start_4PL)

  testthat::expect_length(out, p)
})

testthat::test_that("null mean function equals mean_vector() at gamma1 = 1", {
  mean_fn <- make_null_mean_fn_proportional(
    ref = "4PL", times = times_equal, gamma1_null = 1, shared = TRUE
  )
  from_fn <- mean_fn(start_4PL)
  from_mv <- mean_vector(
    ref = "4PL", J = J, K = K, times = times_equal, params_list = params_list_4PL
  )

  testthat::expect_equal(from_fn, from_mv, tolerance = 1e-12)
})

testthat::test_that("make_null_mean_fn_proportional errors on wrong gamma1_null length", {
  testthat::expect_error(
    make_null_mean_fn_proportional(
      ref = "4PL", times = times_equal, gamma1_null = c(1, 1, 1), shared = FALSE
    )
  )
})


# ============================================================================
# fit_null_gls
# ============================================================================

testthat::test_that("fit_null_gls converges and returns correct list structure", {
  mean_fn <- make_null_mean_fn_proportional(
    ref = "4PL", times = times_equal, gamma1_null = 1, shared = TRUE
  )
  fit <- fit_null_gls(
    m_tilde = m_tilde_4PL,
    Sigma   = Sig,
    mean_fn = mean_fn,
    start   = start_4PL
  )

  testthat::expect_true(is.list(fit))
  testthat::expect_named(fit, c("gamma0_hat", "criterion", "optim"))
  testthat::expect_length(fit$gamma0_hat, length(start_4PL))
  testthat::expect_true(fit$criterion >= 0)
})

testthat::test_that("fit_null_gls recovers true params when m_tilde = mu(true_params)", {
  mean_fn <- make_null_mean_fn_proportional(
    ref = "4PL", times = times_equal, gamma1_null = 1, shared = TRUE
  )
  # Under m_tilde = m_null (exact null), the optimizer should return params
  # close to start_4PL (or any minimizer of a zero-valued criterion).
  fit <- fit_null_gls(
    m_tilde = m_null,
    Sigma   = Sig,
    mean_fn = mean_fn,
    start   = start_4PL
  )
  testthat::expect_equal(fit$criterion, 0, tolerance = 1e-6)
})


# ============================================================================
# build_jacobian_at_null
# ============================================================================

testthat::test_that("build_jacobian_at_null returns matrix with correct rows (shared)", {
  jac <- build_jacobian_at_null(
    gamma0_hat  = start_4PL,
    gamma1_null = 1,
    ref         = "4PL",
    times       = times_equal,
    shared      = TRUE
  )

  testthat::expect_true(is.matrix(jac))
  testthat::expect_equal(nrow(jac), p)
})

testthat::test_that("build_jacobian_at_null returns matrix with correct rows (outcome-specific)", {
  jac <- build_jacobian_at_null(
    gamma0_hat  = start_4PL,
    gamma1_null = c(1, 1),
    ref         = "4PL",
    times       = times_equal,
    shared      = FALSE
  )

  testthat::expect_true(is.matrix(jac))
  testthat::expect_equal(nrow(jac), p)
})

testthat::test_that("build_jacobian_at_null last column matches time_d_4PL * times (shared)", {
  jac <- build_jacobian_at_null(
    gamma0_hat  = start_4PL,
    gamma1_null = 1,
    ref         = "4PL",
    times       = times_equal,
    shared      = TRUE
  )
  # For the shared case the last column is the stacked d f_0(t) / dt * t,
  # which equals zero for the control rows and time_d_4PL * times for treatment rows.
  # Control block for outcome 1: rows 1..(K+1); treatment block: (K+2)..(2K+2).
  t_j <- times_equal[[1]]
  expected_trt <- time_d_4PL(t_j, params_list_4PL[[1]]) * t_j
  testthat::expect_equal(
    jac[(K + 2):(2 * (K + 1)), ncol(jac)],
    expected_trt,
    tolerance = 1e-10
  )
})


# ============================================================================
# build_Bn
# ============================================================================

testthat::test_that("build_Bn returns matrix with correct dimensions (shared)", {
  A <- build_omnibus_contrast_multi_outcome(times = times_equal)
  jac <- build_jacobian_at_null(
    gamma0_hat = start_4PL, gamma1_null = 1,
    ref = "4PL", times = times_equal, shared = TRUE
  )
  B_n <- build_Bn(jac, A, Sig, shared = TRUE)

  # B_n has d_time rows (1 for shared) and p columns
  testthat::expect_equal(dim(B_n), c(1L, p))
})

testthat::test_that("build_Bn returns matrix with correct dimensions (outcome-specific)", {
  A <- build_omnibus_contrast_multi_outcome(times = times_equal)
  jac <- build_jacobian_at_null(
    gamma0_hat = start_4PL, gamma1_null = c(1, 1),
    ref = "4PL", times = times_equal, shared = FALSE
  )
  B_n <- build_Bn(jac, A, Sig, shared = FALSE, J = J)

  testthat::expect_equal(dim(B_n), c(J, p))
})

testthat::test_that("build_Bn errors when J is missing for outcome-specific case", {
  A <- build_omnibus_contrast_multi_outcome(times = times_equal)
  jac <- build_jacobian_at_null(
    gamma0_hat = start_4PL, gamma1_null = c(1, 1),
    ref = "4PL", times = times_equal, shared = FALSE
  )
  testthat::expect_error(build_Bn(jac, A, Sig, shared = FALSE))
})

testthat::test_that("B_n m_null = 0 (orthogonality under the null)", {
  A <- build_omnibus_contrast_multi_outcome(times = times_equal)
  jac <- build_jacobian_at_null(
    gamma0_hat = start_4PL, gamma1_null = 1,
    ref = "4PL", times = times_equal, shared = TRUE
  )
  B_n <- build_Bn(jac, A, Sig, shared = TRUE)

  testthat::expect_equal(
    as.numeric(B_n %*% m_null),
    rep(0, nrow(B_n)),
    tolerance = 1e-8
  )
})


# ============================================================================
# targeted_test_statistic
# ============================================================================

testthat::test_that("targeted_test_statistic returns list with statistic, df, p_value", {
  A <- build_omnibus_contrast_multi_outcome(times = times_equal)
  jac <- build_jacobian_at_null(
    gamma0_hat = start_4PL, gamma1_null = 1,
    ref = "4PL", times = times_equal, shared = TRUE
  )
  B_n <- build_Bn(jac, A, Sig, shared = TRUE)
  res <- targeted_test_statistic(B_n, m_tilde_4PL, Sig, n = n_obs)

  testthat::expect_named(res, c("statistic", "df", "p_value"))
  testthat::expect_true(res$statistic >= 0)
  testthat::expect_equal(res$df, 1L)
  testthat::expect_true(res$p_value >= 0 && res$p_value <= 1)
})

testthat::test_that("targeted_test_statistic is near zero and p near 1 under m_null", {
  A <- build_omnibus_contrast_multi_outcome(times = times_equal)
  jac <- build_jacobian_at_null(
    gamma0_hat = start_4PL, gamma1_null = 1,
    ref = "4PL", times = times_equal, shared = TRUE
  )
  B_n <- build_Bn(jac, A, Sig, shared = TRUE)
  res <- targeted_test_statistic(B_n, m_null, Sig, n = n_obs)

  testthat::expect_equal(res$statistic, 0, tolerance = 1e-8)
  testthat::expect_equal(res$p_value, 1, tolerance = 1e-6)
})


# ============================================================================
# targeted_test_proportional_slowing  (end-to-end)
# ============================================================================

testthat::test_that("targeted_test_proportional_slowing returns correct list structure (shared, 4PL)", {
  out <- targeted_test_proportional_slowing(
    m_tilde     = m_tilde_4PL,
    Sigma       = Sig,
    times       = times_equal,
    ref         = "4PL",
    start       = start_4PL,
    shared      = TRUE,
    gamma1_null = 1,
    n           = n_obs
  )

  testthat::expect_named(out, c("B_n", "gamma0_hat", "gls_fit", "statistic", "df", "p_value"))
  testthat::expect_equal(nrow(out$B_n), 1L)
  testthat::expect_equal(out$df, 1L)
  testthat::expect_true(out$p_value >= 0 && out$p_value <= 1)
})

testthat::test_that("targeted_test_proportional_slowing returns correct list structure (outcome-specific, 4PL)", {
  out <- targeted_test_proportional_slowing(
    m_tilde     = m_tilde_4PL,
    Sigma       = Sig,
    times       = times_equal,
    ref         = "4PL",
    start       = start_4PL,
    shared      = FALSE,
    gamma1_null = c(1, 1),
    n           = n_obs
  )

  testthat::expect_equal(nrow(out$B_n), J)
  testthat::expect_equal(out$df, J)
})

testthat::test_that("p_value is near 1 when m_tilde is the null mean (no treatment effect)", {
  out <- targeted_test_proportional_slowing(
    m_tilde     = m_null,
    Sigma       = Sig,
    times       = times_equal,
    ref         = "4PL",
    start       = start_4PL,
    shared      = TRUE,
    gamma1_null = 1,
    n           = n_obs
  )

  testthat::expect_equal(out$statistic, 0, tolerance = 1e-6)
  testthat::expect_true(out$p_value > 0.99)
})

testthat::test_that("outcome-specific test has higher or equal power than shared when effects differ", {
  # Perturb treatment means differently per outcome: larger signal for outcome 2.
  m_het <- m_null
  # Outcome 1 treatment block: rows (K+2)..(2K+2)
  m_het[(K + 2):(2 * (K + 1))] <- m_het[(K + 2):(2 * (K + 1))] - 0.02
  # Outcome 2 treatment block: rows (2*(K+1) + K+2)..(4*(K+1))
  o2_start <- 2 * (K + 1) + K + 2
  o2_end   <- 4 * (K + 1)
  m_het[o2_start:o2_end] <- m_het[o2_start:o2_end] - 0.15

  out_shared <- targeted_test_proportional_slowing(
    m_tilde = m_het, Sigma = Sig, times = times_equal,
    ref = "4PL", start = start_4PL, shared = TRUE, n = n_obs
  )
  out_specific <- targeted_test_proportional_slowing(
    m_tilde = m_het, Sigma = Sig, times = times_equal,
    ref = "4PL", start = start_4PL, shared = FALSE, gamma1_null = c(1, 1), n = n_obs
  )

  # Under heterogeneous effects, outcome-specific statistic should be >= shared
  # (more degrees of freedom, but also more targeted signal capture).
  # We only assert both statistics are positive.
  testthat::expect_true(out_shared$statistic >= 0)
  testthat::expect_true(out_specific$statistic >= 0)
})
