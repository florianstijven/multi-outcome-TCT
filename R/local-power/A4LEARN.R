library(tidyverse)
library(A4LEARN)


A4LEARN::ADQS$QSTEST %>%
  unique()

ADQS_PACC <- A4LEARN::ADQS %>%
  dplyr::filter(MITTFL== 1) %>%
  dplyr::filter(EPOCH == "BLINDED TREATMENT" | AVISIT == "006") %>%
  dplyr::filter(QSTESTCD == "PACC") %>%
  rename(PACC = QSSTRESN) %>%
  select(BID, ASEQNCS, TX, ADURW, TX, AGEYR, 
         AAPOEGNPRSNFLG, EDCCNTU, SUVRCER, QSVERSION, PACC) %>%
  mutate(TX = factor(TX, levels = c("Placebo", "Solanezumab"))) %>%
  arrange(BID, ADURW) %>%
  na.omit()

ggplot(ADQS_PACC, aes(x=ADURW, y=PACC, color=TX)) +
  geom_line(aes(group = BID), alpha=0.2) +
  geom_smooth(se = FALSE) +
  theme(legend.position = "inside", legend.position.inside = c(0.2, 0.2)) +
  xlab("Weeks since Randomization") +
  scale_x_continuous(breaks = seq(0, max(ADQS_PACC$ADURW), by = 24))


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


clinical_data_summary <- clinical_data %>%
  group_by(QSTESTCD, TX, weeks_since_randomization) %>%
  summarise(
    mean_score = mean(QSSTRESN, na.rm = TRUE),
    sd_score = sd(QSSTRESN, na.rm = TRUE),
    n = n(),
    se_score = sd_score / sqrt(n)
  ) %>%
  ungroup()


clinical_data_summary %>%
  ggplot(aes(x=weeks_since_randomization, y=mean_score, color=TX)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean_score-se_score, ymax=mean_score+se_score), width=0.2) +
  facet_wrap(. ~ QSTESTCD, scales = "free") +
  xlab("Weeks since Randomization") +
  ylab("Mean Score")

clinical_data %>%
  ggplot(aes(x=ADURW, y=QSSTRESN, color=TX)) +
  # geom_line(aes(group = BID), alpha=0.1) +
  geom_smooth(se = TRUE) +
  facet_wrap(. ~ QSTESTCD, scales = "free")


view((A4LEARN::mmse))

A4LEARN::mmse %>%
  pivot_longer(cols = c(8:20, 28:39), names_to = "MMSE_Item", values_to = "Score") %>%
  mutate(Score = Score == "Correct") %>%
  group_by(MMSE_Item) %>%
  summarise(Mean_Score = mean(Score, na.rm = TRUE)) %>%
  print(n = 100)

MMSE_tbl = A4LEARN::mmse %>%
  pivot_longer(cols = c(8:20, 28:39),
               names_to = "MMSE_Item",
               values_to = "Score") %>%
  mutate(Score = Score == "Correct", Score = as.numeric(Score)) %>%
  left_join(
    clinical_data %>% select(BID, VISITCD, TX, ADURW, EPOCH, weeks_since_randomization) %>%
      mutate(VISCODE = as.double(VISITCD)),
    by = c("BID", "VISCODE")
  ) %>%
  filter(!is.na(TX))

MMSE_summary_tbl <- MMSE_tbl %>%
  group_by(MMSE_Item, TX, weeks_since_randomization) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    sd_score = sd(Score, na.rm = TRUE),
    n = n(),
    se_score = sd_score / sqrt(n)
  ) %>%
  ungroup()


MMSE_tbl %>%
  ggplot(aes(x = ADURW, y = Score, color = TX)) +
  geom_smooth(se = FALSE, method = "lm", formula = y ~ poly(x, 2)) +
  facet_wrap(. ~ MMSE_Item)

MMSE_summary_tbl %>%
  filter(weeks_since_randomization >= 0) %>%
  ggplot(aes(x = weeks_since_randomization, y = mean_score, color = TX)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=mean_score - 1.96 * se_score, ymax=mean_score + 1.95 * se_score), width=0.2) +
  facet_wrap(. ~ MMSE_Item, scales = "free_y") +
  xlab("Weeks since Randomization") +
  ylab("Mean Score") +
  theme(legend.position = "bottom")


CDRSB_tbl = A4LEARN::cdr %>%
  pivot_longer(cols = c(10:17),
               names_to = "CDR_Item",
               values_to = "Score") %>%
  left_join(
    clinical_data %>% select(BID, VISITCD, weeks_since_randomization, TX, ADURW) %>%
      group_by(BID, VISITCD) %>%
      slice_head() %>%
      mutate(VISCODE = as.double(VISITCD)),
    by = c("BID", "VISCODE")
  ) %>%
  filter(!is.na(TX))

CDRSB_summary_tbl <- CDRSB_tbl %>%
  group_by(CDR_Item, TX, weeks_since_randomization) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    sd_score = sd(Score, na.rm = TRUE),
    n = n(),
    se_score = sd_score / sqrt(n)
  ) %>%
  ungroup()

CDRSB_summary_tbl %>%
  ggplot(aes(x = weeks_since_randomization, y = mean_score, color = TX)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(
    ymin = mean_score - 1.96 * se_score,
    ymax = mean_score + 1.95 * se_score
  ),
  width = 0.2) +
  facet_wrap(. ~ CDR_Item, scales = "free_y") +
  xlab("Weeks since Randomization") +
  ylab("Mean Score") +
  theme(legend.position = "bottom")
  

CDRSB_tbl %>%
  ggplot(aes(x = ADURW, y = Score, color = TX)) +
  geom_smooth(se = FALSE) +
  facet_wrap(. ~ CDR_Item)

times <- CDRSB_tbl %>%
  pull(weeks_since_randomization) %>%
  unique()

m_tilde_CDRSB <- CDRSB_tbl %>%
  pivot_wider(
    names_from = CDR_Item,
    values_from = Score,
    names_prefix = "SCORE_"
  ) %>%
  group_by(TX, weeks_since_randomization) %>%
  summarise(across(
    starts_with("SCORE_"),
    mean,
    na.rm = TRUE,
    .names = "mean_{col}"
  )) %>%
  arrange(TX, weeks_since_randomization)

m_tilde_CDRSB <- m_tilde_CDRSB %>% pull(m_tilde_CDRSB)


Sigma_CDRSB_Placebo <- CDRSB_tbl %>%
  filter(TX == "Placebo") %>%
  pivot_wider(names_from = "CDR_Item",
              values_from = "Score",
              names_prefix = "SCORE_") %>%
  select(
    -c(
      EPOCH,
      SUBSTUDY,
      VISCODE,
      NDREASON,
      CDSPSRCE,
      CDSPVERS,
      AVISIT,
      BPID,
      CDADTC_DAYS_T0,
      CDADTC_DAYS_CONSENT,
      CDDY,
      CDOLEEVENT,
      CDRSB,
      DONE,
      CDPTSRCE,
      ADURW,
      VISITCD,
      CDEVENT
    )
  ) %>%
  pivot_wider(
    names_from = "weeks_since_randomization",
    values_from = starts_with("SCORE_"),
    names_glue = "{.value}_WEEK_{weeks_since_randomization}"
  ) %>%
  select(contains("WEEK_")) %>%
  cov(use = "pair")

Sigma_CDRSB_TRT <- CDRSB_tbl %>%
  filter(TX == "Solanezumab") %>%
  pivot_wider(names_from = "CDR_Item",
              values_from = "Score",
              names_prefix = "SCORE_") %>%
  select(
    -c(
      EPOCH,
      SUBSTUDY,
      VISCODE,
      NDREASON,
      CDSPSRCE,
      CDSPVERS,
      AVISIT,
      BPID,
      CDADTC_DAYS_T0,
      CDADTC_DAYS_CONSENT,
      CDDY,
      CDOLEEVENT,
      CDRSB,
      DONE,
      CDPTSRCE,
      ADURW,
      VISITCD,
      CDEVENT
    )
  ) %>%
  pivot_wider(
    names_from = "weeks_since_randomization",
    values_from = starts_with("SCORE_"),
    names_glue = "{.value}_WEEK_{weeks_since_randomization}"
  ) %>%
  select(contains("WEEK_")) %>%
  cor(use = "pair")
