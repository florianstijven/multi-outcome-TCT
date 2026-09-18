# Packages and functions ---------

library(tidyverse)
library(furrr)
library(tctHelpers)

source("R/simulations/simulation-A4LEARN-setup.R")

# Simulations options. This should later be specified through command line
# arguments or a configuration file.

# Number of MC replications.
n_MC <- 10

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}

# Simulations ---------

## Setup ------------

scenarios_dgm_tbl <- readRDS("results/simulations/intermediate-objects/scenarios_dgm_tbl.rds")

## Helper Functions ---------

# The function below is copied from `R/simulations/data-generating-mechanism.R`.
# Later updates could remove this redundancy and define this function only once.
construct_working_model_helper <- function(times, ref, slowing_shared, type = "proportional") {
  working_model <- make_slowing_models(ref = ref,
                                       times = times,
                                       type = type)
  
  if (slowing_shared == "common") {
    working_model <- shared_parameter_model(working_model,
                                            shared_param_positions = working_model$treatment_params_position)
  }
  
  return(working_model)
}

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
  select(outcome, ref, slowing_shared, slowing_factor, alternative_type, p_values) %>%
  saveRDS(file = "results/simulations/intermediate-objects/p_values_tbl.rds")
