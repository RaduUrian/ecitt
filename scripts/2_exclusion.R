########## SETUP ##########

rm(list = ls())

# Sets pseudo random number generator seed to fixed value for reproducibility
set.seed(1)

renv::restore()

options(pkgType = "binary")

# Check if package 'here' is installed; if not, install. 
if (!require("here")) 
{ 
  install.packages("here")
  library(here)
}

if (!require("pacman")) 
{ 
  install.packages("pacman")
  library(here)
}

# Increase timeout for downloading packages in case of slow download speed
options(timeout = 300)

# Load all required packages and update to ensure compatibility
pacman::p_load(
  here,
  dplyr,
  tidyr,
  lubridate,
  hms,
  phsmethods,
  readxl,
  stringr,
  tibble
)

renv::snapshot()

options(scipen = 999) # Remove scientific notations

options(es.use_symbols = TRUE)

ecitt_data <- read.csv(here("data", "ecitt_data_merged.csv"))

##########################################################
## 1. Participant-level summaries
##########################################################

part_summary <- ecitt_data %>%
  group_by(id) %>%
  summarise(
    n_trials = n(),
    n_prpt   = sum(trialVar == "prpt", na.rm = TRUE),
    n_inh    = sum(trialVar == "inhb", na.rm = TRUE),
    acc_prpt = mean(accu[trialVar == "prpt"], na.rm = TRUE),
    acc_inh  = mean(accu[trialVar == "inhb"], na.rm = TRUE),
    fail_all_inh = acc_inh == 0
  ) %>%
  ungroup()

##########################################################
## 2. Apply participant-level exclusion criteria
##########################################################

part_flags <- part_summary %>%
  mutate(
    excl_prpt_lt60   = acc_prpt < 0.60 | is.na(acc_prpt),
    excl_lt20_trials = n_trials < 20,
    excl_lt2_inh     = n_inh < 2,
    
    include_accuracy = !(excl_prpt_lt60 | excl_lt20_trials | excl_lt2_inh),
    
    # For RT analyses, also exclude those who fail all inhibition trials
    include_rt       = include_accuracy & !fail_all_inh
  )

###############################################
## 3. Mark invalid trials (global exclusions)
###############################################

ecitt_data <- ecitt_data %>%
  mutate(
    invalid_global = case_when(
      respTime < 300 
      | reason.invalid == 3 
      | reason.invalid == 6 
      | reason.invalid == 9
      | reason.invalid == 1
      | reason.invalid == 5 ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  relocate(invalid_global,
           .after = id)

# trial-wise exclusions
trial_exclusion_flags <- ecitt_data %>% 
  mutate(
    lt_300ms = respTime < 300,
    gt_5000ms = respTime > 5000,
    code_1 = reason.invalid == 1,
    code_3 = reason.invalid == 3,
    code_5 = reason.invalid == 5,
    code_6 = reason.invalid == 6,
    code_9 = reason.invalid == 9
  )

levels(ecitt_data$reason.invalid)

trial_exclusion_summary <- trial_exclusion_flags %>%
  summarise(
    lt_300ms_n = sum(lt_300ms, na.rm = TRUE),
    lt_300ms_pct = 100*mean(lt_300ms, na.rm = TRUE),
    gt_5000ms_n = sum(gt_5000ms, na.rm = TRUE),
    gt_5000ms_pct = 100*mean(gt_5000ms, na.rm = TRUE),
    code_1_n = sum(code_1, na.rm = TRUE),
    code_1_pct = 100*mean(code_1, na.rm = TRUE),
    code_3_n = sum(code_3, na.rm = TRUE),
    code_3_pct = 100*mean(code_3, na.rm = TRUE),
    code_5_n = sum(code_5, na.rm = TRUE),
    code_5_pct = 100*mean(code_5, na.rm = TRUE),
    code_6_n = sum(code_6, na.rm = TRUE),
    code_6_pct = 100*mean(code_6, na.rm = TRUE),
    code_9_n = sum(code_9, na.rm = TRUE),
    code_9_pct = 100*mean(code_9, na.rm = TRUE),
  )

# total excluded trials for rt analysis
excluded_trials_rt <- trial_exclusion_flags %>% 
  filter(lt_300ms == TRUE 
         | gt_5000ms == TRUE 
         | code_1 == TRUE 
         | code_3 == TRUE
         | code_5 == TRUE
         | code_6 == TRUE
         | code_9 == TRUE)

# total excluded trials for acc analysis
excluded_trials_acc <- trial_exclusion_flags %>% 
  filter(lt_300ms == TRUE
         | code_1 == TRUE 
         | code_3 == TRUE
         | code_5 == TRUE
         | code_6 == TRUE
         | code_9 == TRUE)

##########################################################
## 4. Accuracy dataset: only global invalid trials removed
##########################################################

ecitt_data_acc <- ecitt_data %>% 
  filter(!invalid_global)

##########################################################
## 5. RT dataset: global invalid + RT > 5000 ms removed
##########################################################

ecitt_data_rt <- ecitt_data %>%
  mutate(invalid_rt = invalid_global | respTime > 5000) %>%
  filter(!invalid_rt)

##########################################################
## 6. Final datasets
##########################################################

ecitt_data_acc_final <- ecitt_data_acc %>%
  left_join(part_flags %>% select(id, include_accuracy),
            by = "id") %>%
  filter(include_accuracy)

ecitt_data_rt_final <- ecitt_data_rt %>%
  left_join(part_flags %>% select(id, include_rt),
            by = "id") %>%
  filter(include_rt & accu == 1)

##########################################################
## 7. Counts of excluded participants and trials for manuscript
##########################################################

exclusion_summary <- part_flags %>% 
  summarize(
    excl_lt2_inh_n = sum(excl_lt2_inh, na.rm = TRUE),
    excl_lt2_inh_pct = 100*mean(excl_lt2_inh, na.rm = TRUE),
    excl_lt20_trials_n = sum(excl_lt20_trials, na.rm = TRUE),
    excl_lt20_trials_pct = 100*mean(excl_lt20_trials, na.rm = TRUE),
    excl_prpt_lt60_n = sum(excl_prpt_lt60, na.rm = TRUE),
    excl_prpt_lt60_pct = 100*mean(excl_prpt_lt60, na.rm = TRUE),
    fail_all_inh_n = sum(fail_all_inh, na.rm = TRUE),
    fail_all_inh_pct = 100*mean(fail_all_inh, na.rm = TRUE)
  )

excluded_participants <- part_flags %>% 
  filter(excl_lt2_inh == TRUE | excl_lt20_trials == TRUE | excl_prpt_lt60 == TRUE | fail_all_inh == TRUE)

##########################################################
# 9. Final dataset
##########################################################


# select relevant columns for analysis

ecitt_data_rt_final <- ecitt_data_rt_final %>%
  select(
    id,
    age.y,
    gender.y,
    trialVar,
    accu,
    respTime,
    income,
    domlanguagepercent,
    nondomlanguagepercent,
    category25,
    trialNo
  )

ecitt_data_acc_final <- ecitt_data_acc_final %>%
  select(
    id,
    age.y,
    gender.y,
    trialVar,
    accu,
    respTime,
    income,
    domlanguagepercent,
    nondomlanguagepercent,
    category25,
    trialNo
    )

# centering age

ecitt_data_rt_final$age_c <- as.numeric(scale(ecitt_data_rt_final$age.y, scale = FALSE))
ecitt_data_acc_final$age_c <- as.numeric(scale(ecitt_data_acc_final$age.y, scale = FALSE))

# centering nondonlanguagepercent

ecitt_data_rt_final$nondom_c <- scale(ecitt_data_rt_final$nondomlanguagepercent, scale = FALSE)
ecitt_data_acc_final$nondom_c <- scale(ecitt_data_acc_final$nondomlanguagepercent, scale = FALSE)

# export csv
write.csv(ecitt_data_acc_final,
          here("data", "ecitt_data_acc_final.csv"),
          row.names = FALSE)

write.csv(ecitt_data_rt_final,
          here("data", "ecitt_data_rt_final.csv"),
          row.names = FALSE)
