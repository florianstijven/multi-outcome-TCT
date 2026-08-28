# Packages ----

library(tidyverse)
library(A4LEARN)

# A4LEARN Data ----

# For now, we will based the DGM on the complete cases in the A4LEARN study.
# This should be changed later.

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

## MMSE --------

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

problematic_items_MMSE <- c("MMAREA",
                            "MMBALLDL",
                            "MMFLAGDL",
                            "MMHOSPIT",
                            "MMBALL",
                            "MMHAND")

# Exclude the problematic items from the MMSE data set.
MMSE_tbl <- MMSE_tbl %>%
  filter(!item %in% problematic_items_MMSE)

# Remove subjects that don't have complete data for all MMSE subitems for all
# time points.
time_points_MMSE <- sort(unique(MMSE_tbl$weeks_since_randomization))
MMSE_tbl_complete_cases <- MMSE_tbl %>%
  group_by(BID) %>%
  filter(all(time_points_MMSE %in% weeks_since_randomization) &
           all(!is.na(score))) %>%
  ungroup()

MMSE_tbl_complete_cases_control <- MMSE_tbl_complete_cases %>%
  filter(TX == "Placebo")




MMSE_summary_tbl <- MMSE_tbl_complete_cases %>%
  group_by(item, TX, weeks_since_randomization) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    n = n(),
    se_score = sd_score / sqrt(n)
  ) %>%
  ungroup()



MMSE_summary_tbl %>%
  filter(weeks_since_randomization >= 0) %>%
  ggplot(aes(x = weeks_since_randomization, y = mean_score, color = TX)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(
    ymin = mean_score - 1.96 * se_score,
    ymax = mean_score + 1.96 * se_score
  ),
  width = 0.2) +
  facet_wrap(. ~ item) +
  xlab("Weeks since Randomization") +
  ylab("Mean Score") +
  theme(legend.position = "bottom")

## CDR-SB -----------

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

problematic_items_CDRSB <- c("CDGLOBAL", "CDSOB")

# Exclude the problematic items from the CDRSB data set.
CDRSB_tbl <- CDRSB_tbl %>%
  filter(!item %in% problematic_items_CDRSB)

# Remove subjects that don't have complete data for all MMSE subitems for all
# time points.
time_points_CDRSB <- sort(unique(CDRSB_tbl$weeks_since_randomization))
CDRSB_tbl_complete_cases <- CDRSB_tbl %>%
  group_by(BID) %>%
  filter(all(time_points_CDRSB %in% weeks_since_randomization) &
           all(!is.na(score))) %>%
  ungroup()

CDRSB_tbl_complete_cases_control <- CDRSB_tbl_complete_cases %>%
  filter(TX == "Placebo")


CDRSB_summary_tbl <- CDRSB_tbl_complete_cases %>%
  group_by(item, TX, weeks_since_randomization) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    n = n(),
    se_score = sd_score / sqrt(n)
  ) %>%
  ungroup()



CDRSB_summary_tbl %>%
  filter(weeks_since_randomization >= 0) %>%
  ggplot(aes(x = weeks_since_randomization, y = mean_score, color = TX)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(
    ymin = mean_score - 1.96 * se_score,
    ymax = mean_score + 1.96 * se_score
  ),
  width = 0.2) +
  facet_wrap(. ~ item) +
  xlab("Weeks since Randomization") +
  ylab("Mean Score") +
  theme(legend.position = "bottom")

# Helper Functions ---------

# Function that samples iid observations from the A4LEARN empirical distribution
# in the placebo arm.
sample_A4LEARN_control <- function(endpoint) {
  if (!(endpoint %in% c("CDR-SB", "MMSE"))) {
    stop("Unknown endpoint: ", endpoint)
  }
  
  if (endpoint == "CDR-SB") {
    sampled_ids <- sample(
      unique(CDRSB_tbl_complete_cases_control$BID),
      size = length(unique(CDRSB_tbl_complete_cases_control$BID)),
      replace = TRUE
    )
    
    purrr::map2(.x = sampled_ids, .y = seq_along(sampled_ids), function(id, new_id) {
      CDRSB_tbl_complete_cases_control %>%
        filter(BID == id) %>%
        mutate(BID = new_id)
    }) %>%
      list_rbind()
  } else if (endpoint == "MMSE") {
    sampled_ids <- sample(
      unique(MMSE_tbl_complete_cases_control$BID),
      size = length(unique(MMSE_tbl_complete_cases_control$BID)),
      replace = TRUE
    )
    
    purrr::map2(.x = sampled_ids, .y = seq_along(sampled_ids), function(id, new_id) {
      MMSE_tbl_complete_cases_control %>%
        filter(BID == id) %>%
        mutate(BID = new_id)
    }) %>%
      list_rbind()
  }
}

sample_A4LEARN_null <- function(endpoint) {
  bind_rows(
    sample_A4LEARN_control(endpoint) %>%
      mutate(TX = "Placebo"),
    sample_A4LEARN_control(endpoint) %>% mutate(TX = "Experimental", BID = BID + nrow(.))
  )
}

data_set <- sample_A4LEARN_null("MMSE")
data_set <- sample_A4LEARN_null("CDR-SB")

analyze_A4LEARN <- function(data_set) {
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
  
  data_set_wide_cc <- data_set %>%
    group_by(BID) %>%
    filter(all(!is.na(score))) %>%
    ungroup() %>%
    pivot_wider(
      names_from = c("item"),
      values_from = "score",
      names_prefix = "SCORE_"
    ) %>%
    select(c(
      starts_with("SCORE_"),
      "weeks_since_randomization",
      "TX",
      "BID"
    )) %>%
    pivot_wider(
      names_from = c("weeks_since_randomization", "TX"),
      values_from = starts_with("SCORE_"),
      names_glue = "{.value}_WEEK_{weeks_since_randomization}_{TX}"
    ) %>%
    select(contains("WEEK_")) %>%
    select(rownames(m_tilde))
  
  Sigma <- data_set_wide_cc %>%
    cov(use = "pair")
  
  n_cc <- nrow(data_set_wide_cc)
  
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
  
  return(list(
    times = times,
    m_tilde = m_tilde,
    Sigma = Sigma,
    n = n_cc,
    Sigma_n = 2 * Sigma / n_cc
  ))
}

compute_treatment_shift <- function(model, params, times) {
  # In this function, we assume that the first length(times[[1]]) elements of
  # the outcome vector correspond to the control group; the next
  # length(times[[1]]) elements correspond to the experimental group. Next, the
  # same for the second outcome with length(times[[2]]) elements for the control
  # group and length(times[[2]]) elements for the experimental group, and so on.
  mu <- model$mean_fn(params)
  
  shift <- c()
  for (j in seq_along(times)) {
    n_time_points <- length(times[[j]])
    mu_j_indices <- ((j - 1) * 2 * n_time_points + 1):(j * 2 * n_time_points)
    
    control_mean <- mu[mu_j_indices][1:n_time_points]
    experimental_mean <- mu[mu_j_indices][(n_time_points + 1):(2 * n_time_points)]
    
    shift <- c(shift, control_mean - control_mean, experimental_mean - control_mean)
  }
  
  return(shift)
}

shared_prop_slowing_gamma1_f <- function(gamma0_list, shared_slowing_param) {
  params <- c()
  for (j in seq_along(gamma0_list)) {
    params <- c(params, gamma0_list[[j]], shared_slowing_param)
  }
  
  return(params)
}

outcome_dependent_prop_slowing_gamma1_f <- function(gamma0_list, shared_slowing_param_lower, shared_slowing_param_upper) {
  gamma1_vec <- seq(from = shared_slowing_param_lower, to = shared_slowing_param_upper, length.out = length(gamma0_list))
  params <- c()
  for (j in seq_along(gamma0_list)) {
    params <- c(params, gamma0_list[[j]], gamma1_vec[[j]])
  }
  
  return(params)
}
