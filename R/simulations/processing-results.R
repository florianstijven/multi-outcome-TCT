# Packages and functions ---------

library(tidyverse)

# Load results
results_tbl <- readRDS("results/simulations/intermediate-objects/p_values_tbl.rds") %>%
  select(outcome, ref, slowing_shared, slowing_factor, alternative_type, p_values)

# Convert results_tbl into nicer format with the results of one analysis per
# row.
results_tbl <- results_tbl %>%
  rowwise(everything()) %>%
  reframe(p_values) %>%
  rowwise(everything()) %>%
  reframe(p_values) %>%
  select(-p_values)

# Analysis ----------

## Tables

## Figures

# Distribution of p-values for each analysis type.
results_tbl %>%
  filter(outcome == "CDR-SB", slowing_shared == "common") %>%
  ggplot(aes(x = p_value))


results_tbl %>%
  filter(outcome == "MMSE", slowing_shared == "common") %>%
  ggplot(aes(x = p_value, fill = test_ref, linetype = test_slowing_shared)) +
  geom_histogram(alpha = 0.5, position = "identity") +
  facet_grid(ref ~ slowing_factor) +
  theme_bw() +
  labs(
    title = "Distribution of p-values by analysis type",
    x = "p-value",
    y = "Count"
  ) +
  scale_fill_brewer(palette = "Set1", name = "Alternative Type") +
  theme(legend.position = "bottom")
