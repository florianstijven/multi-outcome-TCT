library(tidyverse)
library(A4LEARN)

figures_dir <- "results/local-power/figures/"

# ============================================================================
# Data Preparation
# ============================================================================ 

clinical_data <- A4LEARN::ADQS %>%
  filter(EPOCH %in% c("BLINDED TREATMENT", "SCREENING"))

clinical_data <- clinical_data %>%
  mutate(
    wk_end = str_locate(VISIT, "wk")[, 2],
    weeks_since_randomization = str_sub(VISIT, wk_end + 1, wk_end + 3),
    weeks_since_randomization = as.numeric(weeks_since_randomization),
    weeks_since_randomization = ifelse(str_detect(VISIT, "Baseline"), 0, weeks_since_randomization),
    weeks_since_randomization = ifelse(str_detect(VISIT, "Screening"), -1, weeks_since_randomization)
    ) %>%
  # Remove missing values in weeks_since_randomization. These are visits
  # corresponding to early termination visits.
  filter(!is.na(weeks_since_randomization)) %>%
  # Remove data from week 156 because it has fewer than 20 observations per
  # group for some reason.
  filter(weeks_since_randomization != 156) %>%
  # Remove some tests
  filter(!QSTESTCD %in% c("CDGLOBAL", "CDSOB")) %>%
  # Remove subjects with missing treatment information
  filter(!is.na(TX))

# If there are baseline measurements for a given test,; then we will use the
# baseline measurement as the week 0 measurement. Otherwise, we will use the
# screening measurement as the week 0 measurement.
clinical_data <- clinical_data %>%
  group_by(BID) %>%
  filter(!(
    weeks_since_randomization == 0 &
      any(weeks_since_randomization == 0)
  )) %>%
  mutate(weeks_since_randomization = pmax(0, weeks_since_randomization))
  
# Compute time- and treatment group specific means and standard deviations for
# each test.
clinical_data_summary <- clinical_data %>%
  group_by(QSTESTCD, TX, weeks_since_randomization) %>%
  summarise(
    mean_score = mean(QSSTRESN, na.rm = TRUE),
    sd_score = sd(QSSTRESN, na.rm = TRUE),
    n = n(),
    se_score = sd_score / sqrt(n)
  ) %>%
  ungroup()

# Data set with MMSE subitem scores.
MMSE_tbl = A4LEARN::mmse %>%
  pivot_longer(cols = c(8:20, 28:39),
               names_to = "item",
               values_to = "score") %>%
  mutate(score = score == "Correct", score = as.numeric(score)) %>%
  left_join(
    clinical_data %>% select(BID, VISITCD, weeks_since_randomization, TX, ADURW) %>%
      group_by(BID, VISITCD) %>%
      slice_head() %>%
      mutate(VISCODE = as.double(VISITCD)),
    by = c("BID", "VISCODE")
  ) %>%
  filter(!is.na(TX))

MMSE_summary_tbl <- MMSE_tbl %>%
  group_by(item, TX, weeks_since_randomization) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    n = n(),
    se_score = sd_score / sqrt(n)
  ) %>%
  ungroup()

problematic_items_MMSE <- c("MMAREA", "MMBALLDL", "MMFLAGDL", "MMHOSPIT", "MMBALL", "MMHAND")

# Exclude the problematic items from the MMSE data set.
MMSE_tbl <- MMSE_tbl %>%
  filter(!item %in% problematic_items_MMSE)

# Data set with CDRSB subitem scores.
CDRSB_tbl = A4LEARN::cdr %>%
  pivot_longer(cols = c(10:17),
               names_to = "item",
               values_to = "score") %>%
  left_join(
    clinical_data %>% select(BID, VISITCD, weeks_since_randomization, TX, ADURW) %>%
      group_by(BID, VISITCD) %>%
      slice_head() %>%
      mutate(VISCODE = as.double(VISITCD)),
    by = c("BID", "VISCODE")
  ) %>%
  filter(!is.na(TX))

CDRSB_summary_tbl <- CDRSB_tbl %>%
  group_by(item, TX, weeks_since_randomization) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    n = n(),
    se_score = sd_score / sqrt(n)
  ) %>%
  ungroup()

problematic_items_CDRSB <- c("CDGLOBAL", "CDSOB")

# Exclude the problematic items from the CDRSB data set.
CDRSB_tbl <- CDRSB_tbl %>%
  filter(!item %in% problematic_items_CDRSB)

# ============================================================================
# Data Exploration
# ============================================================================ 

clinical_data_summary %>%
  ggplot(aes(x=weeks_since_randomization, y=mean_score, color=TX)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean_score-se_score, ymax=mean_score+se_score), width=0.2) +
  facet_wrap(. ~ QSTESTCD, scales = "free") +
  xlab("Weeks since Randomization") +
  ylab("Mean Score") +
  theme(legend.position = "bottom")

ggsave(filename = file.path(figures_dir, "mean_scores_over_time.pdf"), width = double_width, height = double_height, dpi = res, unit = unit)

MMSE_summary_tbl %>%
  filter(weeks_since_randomization >= 0) %>%
  ggplot(aes(x = weeks_since_randomization, y = mean_score, color = TX)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean_score - 1.96 * se_score, ymax=mean_score + 1.96 * se_score), width=0.2) +
  facet_wrap(. ~ item) +
  xlab("Weeks since Randomization") +
  ylab("Mean Score") +
  theme(legend.position = "bottom")

ggsave(filename = file.path(figures_dir, "mean_scores_over_time_MMSE.pdf"), width = double_width, height = double_height, dpi = res, unit = unit)


CDRSB_summary_tbl %>%
  ggplot(aes(x = weeks_since_randomization, y = mean_score, color = TX)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(
    ymin = mean_score - 1.96 * se_score,
    ymax = mean_score + 1.96 * se_score
  ),
  width = 0.2) +
  facet_wrap(. ~ item, scales = "free_y") +
  xlab("Weeks since Randomization") +
  ylab("Mean Score") +
  theme(legend.position = "bottom")

ggsave(filename = file.path(figures_dir, "mean_scores_over_time_CDRSB.pdf"), width = double_width, height = double_height, dpi = res, unit = unit)

problematic_items_CDRSB <- c("CARE")

# ============================================================================
# Local Power
# ============================================================================ 

source("R/local-power/helper-functions-local-power.R")

# Function to extract the time points, mean vector, and covariance matrix for a
# given outcome.
extract_information_local_power <- function(data_set) {
  times <- data_set %>%
    pull(weeks_since_randomization) %>%
    unique()
  
  m_tilde <- data_set %>%
    pivot_wider(
      names_from = item,
      values_from = score,
      names_prefix = "SCORE_"
    ) %>%
    group_by(TX, weeks_since_randomization) %>%
    summarise(across(
      starts_with("SCORE_"),
      mean,
      na.rm = TRUE,
      .names = "{col}"
    )) %>%
    arrange(TX, weeks_since_randomization)
  
  m_tilde <- m_tilde %>% 
    pivot_wider(
      names_from = c("weeks_since_randomization", "TX"),
      values_from = starts_with("SCORE_"),
      names_glue = "{.value}_WEEK_{weeks_since_randomization}_{TX}"
    )
  
  m_tilde <- m_tilde %>% t()
  
  Sigma <- data_set %>%
    group_by(BID) %>%
    filter(all(!is.na(score))) %>%
    ungroup() %>%
    pivot_wider(names_from = c("item"),
                values_from = "score",
                names_prefix = "SCORE_") %>%
    select(
      c(starts_with("SCORE_"), "weeks_since_randomization", "TX", "BID")
    ) %>%
    pivot_wider(
      names_from = c("weeks_since_randomization", "TX"),
      values_from = starts_with("SCORE_"),
      names_glue = "{.value}_WEEK_{weeks_since_randomization}_{TX}"
    ) %>%
    select(contains("WEEK_")) %>%
    select(rownames(m_tilde)) %>%
    cov(use = "pair")
  
  # 50% of the covariance values should be NA, not more.
  if (mean(is.na(Sigma)) != 0.5) {
    stop("Unexpected proportion of missing values in covariance matrix.")
  }
  # The NAs are set to zero.
  Sigma[is.na(Sigma)] <- 0
  mean(is.na(Sigma))
  
  # Check whether rows/columns match between covariance matrices and mean vector.
  if (any(rownames(Sigma) != rownames(m_tilde))) {
    stop("Row names of covariance matrix do not match row names of mean vector.")
  }
  
  return(
    list(
      times = times,
      m_tilde = m_tilde,
      Sigma = Sigma
    )
  )
}

fit_4PL_model <- function(data_set, lower, upper) {
  fourPL_fit <- data_set %>%
    filter(TX == "Placebo") %>%
    group_by(item) %>%
    summarise(fourPL_fit = list(
      nls(
        score ~ lower + (upper - lower) /
          (1 + (exp(
            -slope * (weeks_since_randomization - inflection)
          ))),
        data = cur_data(),
        start = list(slope = 0.001, inflection = 1000),
        control = nls.control(maxiter = 500,, warnOnly = TRUE),
        algorithm = "port"
      )
    ), .groups = "drop")
  
  fourPL_fit <- fourPL_fit %>%
    mutate(
      inflection = map_dbl(fourPL_fit, ~ coef(.x)[2]),
      slope = map_dbl(fourPL_fit, ~ coef(.x)[1]),
      lower = .env$lower,
      upper = .env$upper
    ) %>%
    select(-fourPL_fit)
  
  return(fourPL_fit)
}

# Goodness of fit plot for the fitted 4PL model for each item in the data set.
plot_4PL_fit <- function(data_set, fourPL_fit) {
  # Extract range of time points.
  times <- data_set %>%
    pull(ADURW) %>%
    unique()
  grid_times <- seq(min(times), max(times), length.out = 100)
  
  # Obtain predictions from the fitted 4PL model for each item in the data set.
  fourPL_fit <- fourPL_fit %>%
    rowwise(everything()) %>%
    reframe(
      tibble(
        weeks_since_randomization = grid_times,
        predicted_score = lower + (upper - lower) / (1 + exp(-slope * (grid_times - inflection)))
      )
    )
  
  data_set %>%
    filter(TX == "Placebo") %>%
    ggplot(aes(x = ADURW, y = score)) +
    geom_line(stat = "smooth", alpha = 0.4, color = "blue") +
    geom_line(data = fourPL_fit %>% rename(ADURW = weeks_since_randomization), aes(x = ADURW, y = predicted_score)) +
    facet_wrap(. ~ item) +
    xlab("Weeks since Randomization")
  }

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



scenarios_tbl <- tibble(
  outcome = c("CDRSB", "MMSE"),
  data_set = list(CDRSB_tbl, MMSE_tbl)
)

scenarios_tbl <- scenarios_tbl %>%
  mutate(lower = 0,
         upper = ifelse(outcome == "CDRSB", 3, 1)) %>%
  mutate(
    information = map(data_set, extract_information_local_power),
    fourPL_fit = pmap(.l = list(data_set = data_set, lower = lower, upper = upper), fit_4PL_model)
  ) %>%
  mutate(
    times = map(information, ~ .x$times),
    m_tilde = map(information, ~ .x$m_tilde),
    Sigma = map(information, ~ .x$Sigma),
    params_list = map(fourPL_fit, ~ list(inflection = .x$inflection, slope = .x$slope)),
    J = map_dbl(fourPL_fit, ~ nrow(.x)),
    K = map_dbl(times, ~ length(.x) - 1),
    A = list(build_omnibus_contrast_multi_outcome(J = J, K = K))
  )

# Save goodness of fit plots for the fitted 4PL model for each item in the data
# set.
scenarios_tbl %>%
  rowwise() %>%
  summarise(
    plot_4PL_fit = list(plot_4PL_fit(data_set = data_set, fourPL_fit = fourPL_fit)),
    plot_filename = file.path(figures_dir, paste0("4PL_fit_", outcome, ".pdf")),
    save_plot = {
      ggsave(
        filename = plot_filename,
        plot = plot_4PL_fit,
        width = double_width,
        height = double_height,
        dpi = res,
        unit = unit
      )
      TRUE
    }
  ) %>%
  ungroup()

# Ensure the covariance matrices are positive semi-definite. If not, adjust them
# to be positive semi-definite.
scenarios_tbl <- scenarios_tbl %>%
  rowwise() %>%
  mutate(
    Sigma = list(as.matrix(Matrix::nearPD(Sigma)$mat))
  ) %>%
  ungroup()

# Compute local shift vectors for each scenario based on the reference model and local alternative type.
scenarios_tbl <- scenarios_tbl %>%
  rowwise() %>%
  mutate(
    h = 1,
    local_shift = list(local_shift_vector_slowing_outcome(ref = ref, J = J, times = times, h = h, params_list = params_list))
  )

testing_grid <- expand_grid(
  ref_working_model = "4PL",
  working_model = c("slowing_shared", "slowing_outcome", "summing")
)

local_power_setup_grid <- cross_join(scenarios_tbl, testing_grid) %>%
  rowwise(everything()) %>%
  summarise(B_contrast = list(
    working_model_contrast_f(
      J = J,
      K = K,
      type = working_model,
      times = times,
      Sigma = Sigma,
      ref = ref_working_model,
      params_list = params_list
    )
  )) %>%
  ungroup()

h_grid <- c(seq(0, 4, length.out = 1e3), seq(4, 20, length.out = 1e3)) %>% unique()

local_power_grid <- local_power_setup_grid %>%
  rowwise(outcome, working_model) %>%
  reframe(
    compute_power_curve(
      B = B_contrast,
      h_grid = h_grid,
      alpha = 0.05,
      Sigma = Sigma,
      base_effect = local_shift
    )
  )

local_power_grid %>%
  filter(ifelse(outcome == "CDRSB", TRUE, h <= 3)) %>%
  ggplot(aes(x = h, y = power, color = working_model)) +
  geom_line() +
  facet_grid(. ~ outcome, scales = "free") +
  labs(
    title = "Local Power Curves - Shared Slowing",
    x = "Local Shift (h)",
    y = "Power",
    color = "Working Model"
  ) +
  theme(legend.position = "bottom")

ggsave(filename = file.path(figures_dir, "local_power_curves_A4LEARN.pdf"), width = double_width, height = single_height, dpi = res, unit = unit)



