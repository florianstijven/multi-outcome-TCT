# Packages and functions ---------

library(tidyverse)

# Load results
results_tbl <- readRDS("results/simulations/intermediate-objects/p_values_tbl.rds")

# Analysis ----------

## Tables

# Compute type I error rates and power for each analysis type. 
error_rates_tbl = results_tbl %>%
  group_by(outcome,
           ref,
           slowing_shared,
           slowing_factor,
           test_ref,
           test_slowing_shared) %>%
  summarise(type_I_error_rate = mean(p_value < 0.05, na.rm = TRUE),
            proportion_NA = mean(is.na(p_value)),
            .groups = "drop")
## Figures

# Distribution of p-values for each analysis type.
results_tbl %>%
  filter(outcome == "CDR-SB", slowing_shared == "common") %>%
  ggplot(aes(x = p_value))


results_tbl %>%
  filter(outcome == "CDR-SB", slowing_shared == "common") %>%
  ggplot(aes(x = p_value, fill = test_ref, linetype = test_slowing_shared)) +
  geom_histogram(alpha = 0.5, position = "identity") +
  facet_grid(ref ~ slowing_factor) +
  theme_bw() +
  labs(
    title = "Distribution of p-values by analysis type",
    x = "p-value",
    y = "Count"
  ) +
  theme(legend.position = "bottom")
