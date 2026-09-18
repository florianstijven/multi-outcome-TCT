# Packages and functions ---------


library(tidyverse)
library(furrr)
library(tctHelpers)

source("R/simulations/simulation-A4LEARN-setup.R")

figures_dir <- "results/local-power/figures/"

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
    working_model <- shared_parameter_model(working_model,
                                            shared_param_positions = working_model$treatment_params_position)
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


alternatives_tbl <- expand_grid(slowing_factor = c(0, 0.25, 0.5, 0.75, 1))


scenarios_dgm_tbl <- cross_join(scenarios_setup_tbl, alternatives_tbl) %>%
  mutate(
    treatment_effect_parameter = purrr::pmap(
      .l = list(slowing_factor, slowing_shared, J),
      .f = function(slowing_factor, slowing_shared, J) {
        if (slowing_shared == "common") {
          slowing_factor
        } else {
          1 - (1 - slowing_factor) * ((0:(J - 1)) / (J - 1))
        }
      }
    ),
    params = purrr::pmap(
      .l = list(gls_fitted, treatment_effect_parameter, slowing_shared),
      .f = function(gls_fitted,
                    treatment_effect_parameter,
                    slowing_shared) {
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
        compute_treatment_shift(model = working_model,
                                params = params,
                                times = times_list)
      }
    )
  )

# Plot of the trajectories for the DGMs under the various scenarios. For every
# DGM, we plot the mean trajectories for the control and treatment groups, for
# each outcome.
plot_dgm_trajectories <- function(gls_fitted, params, times) {
  outcome_strata = stringr::str_split_i(
    string = rownames(gls_fitted$data$m_tilde),
    pattern = "_",
    i = 2
  )
  treatment_strata = stringr::str_split_i(
    string = rownames(gls_fitted$data$m_tilde),
    pattern = "_",
    i = 5
  )
  
  # Set parameter values to what is specific in `params`.
  gls_fitted$coefficients = params
  # We also set the null_model to FALSE, so that the fitted model is not treated
  # as a null model after having added possibly non-null parameter values.
  gls_fitted$null_model = FALSE
  
  plot_gls_fitted(
    gls_fitted = gls_fitted,
    times = times,
    outcome_strata = outcome_strata,
    treatment_strata = treatment_strata,
    first_stage = FALSE
  ) +
    facet_wrap( ~ outcome_strata)
}

scenarios_dgm_tbl %>%
  rowwise(everything()) %>%
  summarise(
    plot = plot_dgm_trajectories(
      gls_fitted = gls_fitted,
      params = params,
      times = rep(unlist(times_list), 2)
    ) %>% list(),
    plot_file_name = paste0(
      "results/simulations/figures/dgm-trajectories-",
      ref,
      "-",
      slowing_shared,
      "-",
      slowing_factor,
      "-",
      outcome,
      ".pdf"
    ),
    plot_file_name = {
      ggsave(
        filename = plot_file_name,
        plot = plot +
          ggtitle(
            paste0(
              "DGM Trajectories: ",
              ref,
              " - ",
              slowing_shared,
              " - ",
              slowing_factor,
              "-",
              outcome
            )
          ),
        width = 8,
        height = 6
      )
      plot_file_name
    }
  ) %>%
  ungroup()


# Save intermediate objects ---------

saveRDS(scenarios_dgm_tbl, file = "results/simulations/intermediate-objects/scenarios_dgm_tbl.rds")
