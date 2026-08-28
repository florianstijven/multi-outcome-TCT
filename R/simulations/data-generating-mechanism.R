# Packages and functions ---------

library(tidyverse)
library(furrr)

list.files(path = file.path("R", "helper-functions"), pattern = "\\.R$", full.names = TRUE) %>%
  lapply(source, echo = FALSE)

source("R/simulations/simulation-A4LEARN-setup.R")


# Simulations ---------

## MMSE -------------
# 
# K  <- MMSE_summary_tbl$weeks_since_randomization %>% unique() %>% length() - 1
# J  <- MMSE_summary_tbl$item %>% unique() %>% length()
# 
# times_list <- rep(times_MMSE / 240, J)
# 
# prop_slow_models_4PL <- make_slowing_models(ref = "4PL", times = times_list, type = "proportional")
# prop_slow_models_NC <- make_slowing_models(ref = "nc_spline", times = times_list, type = "proportional")
# 
# shared_prop_slow_models_4PL <- shared_parameter_model(
#   model = prop_slow_models_4PL,
#   shared_param_positions = list((1:J) * 3)
# )
# 
# compute_treatment_shift(model = prop_slow_models_4PL, params = c(1, 0.1, 1e3, 1, 1, 1, 1, 1, 1), times = times_equal)
# 
# n_MC <- 30
# p_values <- numeric(n_MC)
# p_values_summing <- numeric(n_MC)
# 
# for (i in 1:n_MC) {
#   data_set <- sample_A4LEARN_null("MMSE")
#   analysis_results <- analyze_A4LEARN(data_set)
#   p_values[i] <- targeted_test(
#     m_tilde = analysis_results$m_tilde,
#     Sigma   = analysis_results$Sigma_n,
#     working_model = shared_prop_slow_models_4PL,
#     A = build_omnibus_contrast_multi_outcome(J, K),
#     ols = TRUE,
#     start = rep(c(5, 0.1), J)
#   )$p_value
#   p_values_summing[i] <- targeted_test_statistic(
#     B = build_summing_contrast_multi_outcome(J, K),
#     m_tilde = analysis_results$m_tilde,
#     Sigma = analysis_results$Sigma_n
#   )$p_value
# }
# 
# hist(p_values, main = "Histogram of p-values from targeted test", xlab = "p-value")
# hist(unlist(p_values_summing), main = "Histogram of p-values from omnibus test", xlab = "p-value")


## CDR-SB -------------


K  <- CDRSB_summary_tbl$weeks_since_randomization %>% unique() %>% length() - 1
J  <- CDRSB_summary_tbl$item %>% unique() %>% length()

times_list <- rep(list(time_points_CDRSB / 240), J)

prop_slow_models_4PL <- make_slowing_models(ref = "4PL", times = times_list, type = "proportional")
prop_slow_models_NC <- make_slowing_models(ref = "nc_spline", times = times_list, type = "proportional")

shared_prop_slow_models_4PL <- shared_parameter_model(model = prop_slow_models_4PL, shared_param_positions = list((1:J) * 3))

n_MC <- 50
p_values <- numeric(n_MC)
p_values_summing <- numeric(n_MC)

CDRSB_tbl_temp <- bind_rows(
  CDRSB_tbl_complete_cases_control %>%
    mutate(TX = "Placebo"),
  CDRSB_tbl_complete_cases_control %>% 
    mutate(TX = "Experimental", BID = paste0(BID, "A"))
)
analysis_results_temp <- analyze_A4LEARN(CDRSB_tbl_temp)

gls_fitted <- two_stage_gls_null(
  m_tilde = analysis_results_temp$m_tilde,
  Sigma   = analysis_results_temp$Sigma_n,
  working_model = shared_prop_slow_models_4PL,
  ols = TRUE,
  start = rep(c(2, 1), J)
)

outcome_names <- sapply(
  rownames(analysis_results_temp$m_tilde),
  function(x) {
    stringr::str_split_i(string = x, pattern = "_", i = 2)
  }
)

treatment_names <- sapply(
  rownames(analysis_results_temp$m_tilde),
  function(x) {
    stringr::str_split_i(string = x, pattern = "_", i = 5)
  }
)

plot_gls_fitted(
  gls_fitted = gls_fitted, 
  treatment_strata = treatment_names, 
  outcome_strata = outcome_names,
  times = rep(unlist(times_list), 2)
)

slowing_factor <- 1

Delta <- compute_treatment_shift(
  model = shared_prop_slow_models_4PL,
  params = c(gls_fitted$gamma_hat, slowing_factor),
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

list.files(path = file.path("R", "helper-functions"), pattern = "\\.R$", full.names = TRUE) %>%
  lapply(source, echo = FALSE)

model

# Set up parallel computing
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}

p_values_list <- future_map(
  .x = 1:n_MC,
  .f = simulate_p_values,
  Delta = Delta,
  .options = furrr_options(
    seed = TRUE,
    stdout = FALSE,
    conditions = character()
  )
)

plan(sequential)
# Errsequential()# Error in `model()`:
#   ! could not find function "model"
# Run `rlang::last_trace()` to see where the error occurred.
# Warning messages:
#   1: In serializedSize(x) :
#   'package:forcats' may not be available when loading
# 2: In serializedSize(x) :
#   'package:lubridate' may not be available when loading


hist(purrr::map_dbl(p_values_list, "p_value_targeted"), main = "Histogram of p-values from targeted test", xlab = "p-value")
hist(purrr::map_dbl(p_values_list, "p_value_summing"), main = "Histogram of p-values from sum test", xlab = "p-value")
mean(purrr::map_dbl(p_values_list, "p_value_targeted") <= 0.05)
mean(purrr::map_dbl(p_values_list, "p_value_summing") <= 0.05)
