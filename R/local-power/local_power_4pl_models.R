#!/usr/bin/env Rscript

# Local power analysis with proportional slowing working models
# Two reference trajectory models: 4PL and Natural Cubic Spline
# Two scenarios: Shared vs Outcome-Specific Slowing
library(tidyverse)

output_dir <- "results/local-power/figures/"

# Set up all scenarios for the local alternatives.
scenario_grid <- expand_grid(
  J = c(3, 6, 12, 24, 48),
  K = c(3, 5, 7),
  rho_time = 0.8,
  rho_outcome = 0.75,
  local_alternative = c("shared", "outcome"),
  ref = "4PL"
)

# Add covariance matrices and omnibus contrast matrices to each local alternative.
scenario_grid <- scenario_grid %>%
  rowwise() %>%
  mutate(
    cov_list = list(build_kronecker_covariance(K = K, J = J, rho_time = rho_time, rho_outcome = rho_outcome)),
    Sigma = list(cov_list$Sigma),
    time_points = list(cov_list$time_points),
    A = list(build_omnibus_contrast_multi_outcome(K, J)),
    params_list = list(list(inflection = 3 + 3 * (1:J - 1) / (2 * J), slope = 1 + (1:J - 1) / J)),
    mean_vec = list(mean_vector(ref = ref, J = J, K = K, times = rep(time_points, 2), params_list = params_list)),
    Sigma = list(unlist(diag(sqrt(mean_vec * (1 - mean_vec))) %*% Sigma %*% diag(sqrt(mean_vec * (1 - mean_vec))))),
  ) %>%
  select(-cov_list)


# Compute local shift vectors for each scenario based on the reference model and local alternative type.
scenario_grid <- scenario_grid %>%
  rowwise() %>%
  mutate(
    h = ifelse(local_alternative == "shared", list(1), list((J / 2) * (1:J) / J)),
    local_shift = list(local_shift_vector_slowing_outcome(ref = ref, J = J, times = time_points, h = h, params_list = params_list))
  )

# Set up all scenarios for the tests whose local power is being evaluated.
working_model_contrast_f <- function(J, K, type, times = NULL, Sigma = NULL, ref = NULL, ...) {
  A = build_omnibus_contrast_multi_outcome(J, K)
  if (type == "omnibus") {
    A
  } else if (type == "linear") {
    return(build_linear_contrast_multi_outcome(J, K, times = times, Sigma = Sigma))
  } else if (type == "slowing_outcome") {
    jacobian <- jacobian_slowing_multiple_outcomes(J = J, K = K, times = times, ref = ref, slowing_only = TRUE, ...)
    return(build_contrast_matrix(jacobian, A = A, Sigma = Sigma))
  } else if (type == "slowing_shared") {
    jacobian <- jacobian_slowing_multiple_outcomes_shared(J = J, K = K, times = times, ref = ref, slowing_only = TRUE, ...)
    return(build_contrast_matrix(jacobian, A = A, Sigma = Sigma))
  } else if (type == "summing") {
    A_sum = build_summing_contrast_multi_outcome(J, K)
    return(A_sum)
  } else {
    stop("Unknown working model type: ", type)
  }
}

testing_grid <- expand_grid(
  ref_working_model = "4PL",
  working_model = c("slowing_shared", "slowing_outcome", "omnibus", "linear", "summing")
)

local_power_setup_grid <- cross_join(scenario_grid, testing_grid) %>%
  rowwise(everything()) %>%
  summarise(B_contrast = list(
    working_model_contrast_f(
      J = J,
      K = K,
      type = working_model,
      times = time_points,
      Sigma = Sigma,
      ref = ref_working_model,
      params_list = params_list
    )
  )) %>%
  ungroup()

h_grid <- c(seq(0, 0.1, length.out = 1e3), seq(0.1, 2, length.out = 1e3)) %>% unique()

local_power_grid <- local_power_setup_grid %>%
  rowwise(J, K, rho_time, rho_outcome, local_alternative, ref, ref_working_model, working_model) %>%
  reframe(
    compute_power_curve(
      B = B_contrast,
      h_grid = h_grid,
      alpha = 0.05,
      Sigma = Sigma,
      base_effect = local_shift
    )
  )

# For plotting, we need to rescale the x-axis with the local scaling. We put all local power plots on the same scale by treating 
# the local power curve for the omnibus test as reference. Specifically, we remove all h_grid values where the local power for the omnibus
# test is larger than 0.99.
local_power_grid_filtered <- local_power_grid %>%
  group_by(J, K, rho_time, rho_outcome, local_alternative, ref, ref_working_model) %>%
  mutate(
    max_h_omnibus = max(h[working_model == "omnibus" & power < 0.99], na.rm = TRUE)
  ) %>%
  filter(h <= max_h_omnibus) %>%
  ungroup()

# Local power curves.
local_power_grid_filtered %>%
  filter(local_alternative == "shared") %>%
  ggplot(aes(x = h, y = power, color = working_model)) +
  geom_line() +
  facet_grid(rows = vars(K), cols = vars(J), labeller = label_both, scales = "free") +
  labs(
    title = "Local Power Curves - Shared Slowing",
    x = "Local Shift (h)",
    y = "Power",
    color = "Working Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, "local_power_curves_shared.pdf"), width = 10, height = 6, dpi = 300)

local_power_grid_filtered %>%
  filter(local_alternative == "outcome") %>%
  ggplot(aes(x = h, y = power, color = working_model)) +
  geom_line() +
  facet_grid(rows = vars(K), cols = vars(J), labeller = label_both, scales = "free") +
  labs(
    title = "Local Power Curves - Outcome-Specific Slowing",
    x = "Local Shift (h)",
    y = "Power",
    color = "Working Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, "local_power_curves_outcome.pdf"), width = 10, height = 6, dpi = 300)


# Relative local power curves.
local_power_grid_filtered %>%
  group_by(J, K, rho_time, rho_outcome, local_alternative, ref, ref_working_model, h) %>%
  mutate(
    power_omnibus = power[working_model == "omnibus"],
    relative_power = power / power_omnibus
  ) %>%
  ungroup() %>%
  filter(local_alternative == "outcome") %>%
  ggplot(aes(x = h, y = relative_power, color = working_model)) +
  geom_line() +
  facet_grid(rows = vars(K), cols = vars(J), labeller = label_both, scales = "free_x") +
  labs(
    title = "Relative Local Power Curves - Outcome-Specific Slowing",
    x = "Local Shift (h)",
    y = "Relative Power (Power / Omnibus Power)",
    color = "Working Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, "relative_local_power_curves_outcome.pdf"), width = 10, height = 6, dpi = 300)

local_power_grid_filtered %>%
  group_by(J, K, rho_time, rho_outcome, local_alternative, ref, ref_working_model, h) %>%
  mutate(
    power_omnibus = power[working_model == "omnibus"],
    relative_power = power / power_omnibus
  ) %>%
  ungroup() %>%
  filter(local_alternative == "shared") %>%
  ggplot(aes(x = h, y = relative_power, color = working_model)) +
  geom_line() +
  facet_grid(rows = vars(K), cols = vars(J), labeller = label_both, scales = "free_x") +
  labs(
    title = "Relative Local Power Curves - Shared Slowing",
    x = "Local Shift (h)",
    y = "Relative Power (Power / Omnibus Power)",
    color = "Working Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, "relative_local_power_curves_shared.pdf"), width = 10, height = 6, dpi = 300)



local_power_grid %>%
  filter(J == 12, K == 5) %>%
  group_by(J, K, rho_time, rho_outcome, local_alternative, ref, ref_working_model, working_model) %>%
  summarise(
    df = mean(df),
    sd = sd(ncp/(h**2), na.rm = TRUE),
    ncp = mean(ncp)
  )


