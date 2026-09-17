# Packages and functions ---------

library(tidyverse)
library(furrr)
library(tctHelpers)

source("R/simulations/simulation-A4LEARN-setup.R")

figures_dir <- "results/local-power/figures/"

# Simulations options. This should later be specified through command line
# arguments or a configuration file.

# Number of MC replications.
n_MC <- 500

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}

# Simulations ---------

## Setup ------------

# Define all simulation scenarios in a tibble. Each row corresponds to a unique
# scenario.

scenarios_setup_tbl <- tibble(
  outcome = c("MMSE", "CDR-SB"),
  true_analysis_list = list(
    analyze_A4LEARN(MMSE_tbl_complete_cases),
    analyze_A4LEARN(CDRSB_tbl_complete_cases)
  )
)

scenarios_setup_tbl <- scenarios_setup_tbl %>%
  mutate(
    times = purrr::map(true_analysis_list, "times"),
    m_0 = purrr::map(true_analysis_list, "m_tilde"),
    Sigma_0 = purrr::map(true_analysis_list, "Sigma_n"),
    outcome_names = purrr::map(
      m_0,
      ~ stringr::str_split_i(
        string = rownames(.x),
        pattern = "_",
        i = 2
      ) %>% unique()
    ),
    J = purrr::map_int(outcome_names, length),
    K = purrr::map_int(times, ~ length(.x) - 1),
    times_list = purrr::map2(times, J, ~ rep(list(.x / 240), .y))
  )


models_setup_tbl <- expand_grid(
  ref = c("4PL", "nc_spline"),
  slowing_shared = c("common", "outcome-specific")
)

construct_working_model_helper <- function(times, ref, slowing_shared, type = "proportional") {
  working_model <- make_slowing_models(ref = ref,
                                       times = times,
                                       type = type)
  
  if (slowing_shared == "common") {
    working_model <- shared_parameter_model(
      working_model, 
      shared_param_positions = working_model$treatment_params_position
      )
  }
  
  return(working_model)
}

scenarios_setup_tbl <- scenarios_setup_tbl %>%
  cross_join(models_setup_tbl) %>%
  rowwise(everything()) %>%
  summarise(
    working_model = construct_working_model_helper(
      times = times_list,
      ref = ref,
      slowing_shared = slowing_shared
    ) %>% list()
  ) %>% ungroup()

scenarios_setup_tbl = scenarios_setup_tbl %>%
  rowwise(everything()) %>%
  summarise(
    gls_fitted = two_stage_gls_null(
      m_tilde = m_0,
      Sigma   = Sigma_0,
      working_model = working_model,
      ols = TRUE,
      start = {
        if (ref == "4PL")
          rep(c(2, 1), J)
        else
          rep(0:K, J)
      }
    ) %>% list()
  ) %>%
  ungroup()


alternatives_tbl <- expand_grid(
  slowing_factor = c(0.25, 0.5, 0.75, 1),
  alternative_type = c("shared", "outcome-specific")
)


scenarios_dgm_tbl <- cross_join(scenarios_setup_tbl, alternatives_tbl) %>%
  mutate(
    treatment_effect_parameter = purrr::pmap(
      .l = list(slowing_factor, slowing_shared, J),
      .f = function(slowing_factor, slowing_shared, J){
        if (slowing_shared == "common") {
          slowing_factor
        } else {
          1 - slowing_factor * ((0:(J - 1)) / (J - 1))
        }
      }),
    params = purrr::pmap(
      .l = list(gls_fitted, treatment_effect_parameter, slowing_shared),
      .f = function(gls_fitted, treatment_effect_parameter, slowing_shared) {
        if (slowing_shared == "common") {
          return(c(coef(gls_fitted), treatment_effect_parameter))
        } else {
          params <- rep(NA, gls_fitted$working_model$no_of_params)
          params[gls_fitted$working_model$treatment_params_position] <- treatment_effect_parameter
          params[gls_fitted$working_model$nuisance_params_position] <- coef(gls_fitted)
          return(params)
        }
      }
    ),
    Delta = purrr::pmap(
      .l = list(working_model, params, times_list),
      .f = function(working_model, params, times_list) {
        compute_treatment_shift(
          model = working_model,
          params = params,
          times = times_list
        )
      }
    )
  )

## Helper Functions ---------

compute_p_value <- function(m_tilde, Sigma, working_model, J, K, ref) {
  if (ref == "4PL") {
    start = rep(c(2, 1), J)
  } else {
    start = rep(0:K, J)
  }
  
  p_value <- NA
  
  try({
    p_value <- targeted_test(
      m_tilde = m_tilde,
      Sigma   = Sigma,
      working_model = working_model,
      A = build_omnibus_contrast_multi_outcome(J, K),
      ols = TRUE,
      start = start
    )$p_value
  }, silent = TRUE)
  
  return(p_value)
}

simulate_p_values <- function(i, outcome, Delta) {
  data_set <- sample_A4LEARN_null(outcome)
  analysis_results <- analyze_A4LEARN(data_set)
  
  times <- analysis_results$times
  
  m_tilde_new <- analysis_results$m_tilde + Delta
  Sigma_new   <- analysis_results$Sigma_n
  
  outcome_names = stringr::str_split_i(string = rownames(m_tilde_new),
                                       pattern = "_",
                                       i = 2) %>% unique()
  J = length(outcome_names)
  K = length(times) - 1
  times_list = rep(list(times / 240), J)

  p_values_tbl <- expand_grid(
    test_ref = c("4PL", "nc_spline"),
    test_slowing_shared = c("common", "outcome-specific")
  ) %>%
    rowwise(everything()) %>%
    summarise(
      working_model = construct_working_model_helper(
        times = times_list,
        ref = test_ref,
        slowing_shared = test_slowing_shared
      ) %>% list(),
      p_value = compute_p_value(
        m_tilde = m_tilde_new,
        Sigma   = Sigma_new,
        working_model = working_model,
        J = J,
        K = K,
        ref = test_ref
      )
    ) %>%
    ungroup() %>%
    select(-working_model)
  
  p_values_tbl <- p_values_tbl %>%
    bind_rows(
      tibble(
        test_ref = "summing",
        test_slowing_shared = "summing",
        p_value = targeted_test_statistic(
          B = build_summing_contrast_multi_outcome(J, K),
          m_tilde = m_tilde_new,
          Sigma = Sigma_new
        )$p_value
      )
    )
  
  return(
    p_values_tbl
  )
}


## Running Simulations ---------


# Helper functions now live in the `tctHelpers` package instead of
# `.GlobalEnv`, so their closures resolve via a stable namespace and furrr
# only needs to know to attach the package on each worker.
scenarios_dgm_tbl %>%
  slice_sample(n = 10) %>%
  rowwise(everything()) %>%
  summarise(
    p_values = future_map(
      .x = 1:n_MC,
      .f = simulate_p_values,
      outcome = outcome,
      Delta = Delta,
      .options = furrr_options(
        seed = TRUE,
        stdout = FALSE,
        conditions = character(),
        packages = "tctHelpers"
      )
    ) %>% list()
  ) %>%
  ungroup() %>%
  saveRDS(file = "results/simulations/intermediate-objects/p_values_tbl.rds")


## MMSE -------------

K  <- MMSE_summary_tbl$weeks_since_randomization %>% unique() %>% length() - 1
J  <- MMSE_summary_tbl$item %>% unique() %>% length()

times_list <- rep(list(time_points_MMSE / 240), J)

analysis_results_temp <- analyze_A4LEARN(MMSE_tbl_complete_cases)

prop_slow_models_4PL <- make_slowing_models(
  ref = "4PL",
  times = times_list,
  type = "proportional",
  outcome_names = stringr::str_split_i(
    string = rownames(analysis_results_temp$m_tilde),
    pattern = "_",
    i = 2
  ) %>% unique()
)
prop_slow_models_NC <- make_slowing_models(
  ref = "nc_spline",
  times = times_list,
  type = "proportional",
  outcome_names = stringr::str_split_i(
    string = rownames(analysis_results_temp$m_tilde),
    pattern = "_",
    i = 2
  ) %>% unique()
)

shared_prop_slow_models_4PL <- shared_parameter_model(
  working_model = prop_slow_models_4PL,
  shared_param_positions = list((1:J) * 3)
  )

n_MC <- 500
p_values <- numeric(n_MC)
p_values_summing <- numeric(n_MC)

MMSE_tbl_temp <- bind_rows(
  MMSE_tbl_complete_cases_control %>%
    mutate(TX = "Placebo"),
  MMSE_tbl_complete_cases_control %>%
    mutate(TX = "Experimental", BID = paste0(BID, "A"))
)


gls_fitted <- two_stage_gls_null(
  m_tilde = analysis_results_temp$m_tilde,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = shared_prop_slow_models_4PL,
  ols = TRUE,
  start = rep(c(2, 1), J)
)

gls_fitted_full <- two_stage_gls_full(
  m_tilde = analysis_results_temp$m_tilde,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = shared_prop_slow_models_4PL,
  ols = TRUE,
  start = c(rep(c(2, 1), J), 1)
)

gls_fitted_full_non_shared <- two_stage_gls_full(
  m_tilde = analysis_results_temp$m_tilde,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = prop_slow_models_4PL,
  ols = TRUE,
  start = rep(c(2, 1, 1), J)
)

summary(gls_fitted_full)
summary(gls_fitted_full_non_shared)

outcome_names <- sapply(rownames(analysis_results_temp$m_tilde), function(x) {
  stringr::str_split_i(string = x,
                       pattern = "_",
                       i = 2)
})

treatment_names <- sapply(rownames(analysis_results_temp$m_tilde), function(x) {
  stringr::str_split_i(string = x,
                       pattern = "_",
                       i = 5)
})

plot_gls_fitted(
  gls_fitted = gls_fitted,
  treatment_strata = treatment_names,
  outcome_strata = outcome_names,
  times = rep(unlist(times_list), 2)
)

slowing_factor <- 0.75

Delta <- compute_treatment_shift(
  model = shared_prop_slow_models_4PL,
  params = c(coef(gls_fitted), slowing_factor),
  times = times_list
)

simulate_p_values <- function(i, Delta) {
  data_set <- sample_A4LEARN_null("MMSE")
  analysis_results <- analyze_A4LEARN(data_set)
  
  p_value_targeted <- targeted_test(
    m_tilde = analysis_results$m_tilde + Delta,
    Sigma   = analysis_results$Sigma_n,
    working_model = shared_prop_slow_models_4PL,
    A = build_omnibus_contrast_multi_outcome(J, K),
    ols = TRUE,
    start = rep(c(2, 1), J)
  )$p_value
  p_value_summing <- targeted_test_statistic(
    B = build_summing_contrast_multi_outcome(J, K),
    m_tilde = analysis_results$m_tilde + Delta,
    Sigma = analysis_results$Sigma_n
  )$p_value
  
  return(list(p_value_targeted = p_value_targeted, p_value_summing = p_value_summing))
}

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}

# Helper functions now live in the `tctHelpers` package instead of
# `.GlobalEnv`, so their closures resolve via a stable namespace and furrr
# only needs to know to attach the package on each worker.
p_values_list <- future_map(
  .x = 1:n_MC,
  .f = simulate_p_values,
  Delta = Delta,
  .options = furrr_options(
    seed = TRUE,
    stdout = FALSE,
    conditions = character(),
    packages = "tctHelpers"
  )
)

plan(sequential)


hist(
  purrr::map_dbl(p_values_list, "p_value_targeted"),
  main = "Histogram of p-values from targeted test",
  xlab = "p-value"
)
hist(
  purrr::map_dbl(p_values_list, "p_value_summing"),
  main = "Histogram of p-values from sum test",
  xlab = "p-value"
)
mean(purrr::map_dbl(p_values_list, "p_value_targeted") <= 0.05)
mean(purrr::map_dbl(p_values_list, "p_value_summing") <= 0.05)

## CDR-SB -------------


K  <- CDRSB_summary_tbl$weeks_since_randomization %>% unique() %>% length() - 1
J  <- CDRSB_summary_tbl$item %>% unique() %>% length()

times_list <- rep(list(time_points_CDRSB / 240), J)

analysis_results_temp <- analyze_A4LEARN(CDRSB_tbl_temp)

prop_slow_models_4PL <- make_slowing_models(
  ref = "4PL",
  times = times_list,
  type = "proportional",
  outcome_names = stringr::str_split_i(
    string = rownames(analysis_results_temp$m_tilde),
    pattern = "_",
    i = 2
  ) %>% unique()
)
prop_slow_models_NC <- make_slowing_models(
  ref = "nc_spline",
  times = times_list,
  type = "proportional",
  outcome_names = stringr::str_split_i(
    string = rownames(analysis_results_temp$m_tilde),
    pattern = "_",
    i = 2
  ) %>% unique()
)

shared_prop_slow_models_4PL <- shared_parameter_model(working_model = prop_slow_models_4PL,
                                                      shared_param_positions = list((1:J) * 3))

n_MC <- 500
p_values <- numeric(n_MC)
p_values_summing <- numeric(n_MC)

CDRSB_tbl_temp <- bind_rows(
  CDRSB_tbl_complete_cases_control %>%
    mutate(TX = "Placebo"),
  CDRSB_tbl_complete_cases_control %>%
    mutate(TX = "Experimental", BID = paste0(BID, "A"))
)


gls_fitted <- two_stage_gls_null(
  m_tilde = analysis_results_temp$m_tilde,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = shared_prop_slow_models_4PL,
  ols = TRUE,
  start = rep(c(2, 1), J)
)

gls_fitted_full <- two_stage_gls_full(
  m_tilde = analysis_results_temp$m_tilde + Delta,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = shared_prop_slow_models_4PL,
  ols = TRUE,
  start = c(rep(c(2, 1), J), 1)
)

gls_fitted_full_non_shared <- two_stage_gls_full(
  m_tilde = analysis_results_temp$m_tilde + Delta,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = prop_slow_models_4PL,
  ols = TRUE,
  start = rep(c(2, 1, 1), J)
)

summary(gls_fitted_full)
summary(gls_fitted_full_non_shared)

outcome_names <- sapply(rownames(analysis_results_temp$m_tilde), function(x) {
  stringr::str_split_i(string = x,
                       pattern = "_",
                       i = 2)
})

treatment_names <- sapply(rownames(analysis_results_temp$m_tilde), function(x) {
  stringr::str_split_i(string = x,
                       pattern = "_",
                       i = 5)
})

plot_gls_fitted(
  gls_fitted = gls_fitted,
  treatment_strata = treatment_names,
  outcome_strata = outcome_names,
  times = rep(unlist(times_list), 2)
)

slowing_factor <- 0.25

Delta <- compute_treatment_shift(
  model = shared_prop_slow_models_4PL,
  params = c(coef(gls_fitted), slowing_factor),
  times = times_list
)

simulate_p_values <- function(i, Delta) {
  data_set <- sample_A4LEARN_null("CDR-SB")
  analysis_results <- analyze_A4LEARN(data_set)
  
  p_value_targeted <- targeted_test(
    m_tilde = analysis_results$m_tilde + Delta,
    Sigma   = analysis_results$Sigma_n,
    working_model = shared_prop_slow_models_4PL,
    A = build_omnibus_contrast_multi_outcome(J, K),
    ols = TRUE,
    start = rep(c(2, 1), J)
  )$p_value
  p_value_summing <- targeted_test_statistic(
    B = build_summing_contrast_multi_outcome(J, K),
    m_tilde = analysis_results$m_tilde + Delta,
    Sigma = analysis_results$Sigma_n
  )$p_value
  
  return(list(p_value_targeted = p_value_targeted, p_value_summing = p_value_summing))
}

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}

# Helper functions now live in the `tctHelpers` package instead of
# `.GlobalEnv`, so their closures resolve via a stable namespace and furrr
# only needs to know to attach the package on each worker.
p_values_list <- future_map(
  .x = 1:n_MC,
  .f = simulate_p_values,
  Delta = Delta,
  .options = furrr_options(
    seed = TRUE,
    stdout = FALSE,
    conditions = character(),
    packages = "tctHelpers"
  )
)

plan(sequential)


hist(
  purrr::map_dbl(p_values_list, "p_value_targeted"),
  main = "Histogram of p-values from targeted test",
  xlab = "p-value"
)
hist(
  purrr::map_dbl(p_values_list, "p_value_summing"),
  main = "Histogram of p-values from sum test",
  xlab = "p-value"
)
mean(purrr::map_dbl(p_values_list, "p_value_targeted") <= 0.05)
mean(purrr::map_dbl(p_values_list, "p_value_summing") <= 0.05)

# ddde ---------

analysis_results_temp <- analyze_A4LEARN(CDRSB_tbl_complete_cases)


gls_fitted <- two_stage_gls_null(
  m_tilde = analysis_results_temp$m_tilde + Delta,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = shared_prop_slow_models_4PL,
  ols = TRUE,
  start = rep(c(2, 1), J)
)

gls_fitted_full <- two_stage_gls_full(
  m_tilde = analysis_results_temp$m_tilde + Delta,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = shared_prop_slow_models_4PL,
  ols = TRUE,
  start = c(rep(c(2, 1), J), 1)
)

gls_fitted_full_non_shared <- two_stage_gls_full(
  m_tilde = analysis_results_temp$m_tilde + Delta,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = prop_slow_models_4PL,
  ols = TRUE,
  start = rep(c(2, 1, 1), J)
)

summary(gls_fitted_full)
summary(gls_fitted_full_non_shared)

outcome_names <- sapply(rownames(analysis_results_temp$m_tilde), function(x) {
  stringr::str_split_i(string = x,
                       pattern = "_",
                       i = 2)
})

treatment_names <- sapply(rownames(analysis_results_temp$m_tilde), function(x) {
  stringr::str_split_i(string = x,
                       pattern = "_",
                       i = 5)
})

plot_gls_fitted(
  gls_fitted = gls_fitted,
  treatment_strata = treatment_names,
  outcome_strata = outcome_names,
  times = rep(unlist(times_list), 2)
)

plot_gls_fitted(
  gls_fitted = gls_fitted_full,
  treatment_strata = treatment_names,
  outcome_strata = outcome_names,
  times = rep(unlist(times_list), 2)
)

plot_gls_fitted(
  gls_fitted = gls_fitted_full_non_shared,
  treatment_strata = treatment_names,
  outcome_strata = outcome_names,
  times = rep(unlist(times_list), 2)
)

targeted_test(
  m_tilde = analysis_results_temp$m_tilde + Delta,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = shared_prop_slow_models_4PL,
  A = build_omnibus_contrast_multi_outcome(J, K),
  ols = TRUE,
  start = rep(c(2, 1), J)
)

targeted_test_statistic(
  B = build_summing_contrast_multi_outcome(J, K),
  m_tilde = analysis_results_temp$m_tilde + Delta,
  Sigma = analysis_results_temp$Sigma_n
)

targeted_test(
  m_tilde = analysis_results_temp$m_tilde + Delta,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = prop_slow_models_4PL,
  A = build_omnibus_contrast_multi_outcome(J, K),
  ols = TRUE,
  start = rep(c(2, 1), J)
)
