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
  tibble,
  naniar,
  devtools
)

devtools::install_github("dustinfife/flexplot")

require(flexplot)

renv::snapshot()

options(scipen = 999) # Remove scientific notations

options(es.use_symbols = TRUE)

########## DATA WRANGLING ##########

# Load raw data

ecitt_list <- list.files(path = here("data", "ecitt_data_individual_complete"),
                       recursive = TRUE,
                       pattern = "\\.csv",
                       full.names = TRUE)
# compile data
ecitt_data <- lapply(ecitt_list, read.csv) 

ecitt_data <- lapply(ecitt_data, as.data.frame)

# fix type of reason.invalid columns
ecitt_data <- lapply(ecitt_data, function(df) {
  if ("reason.invalid" %in% names(df)) {
    df$reason.invalid <- as.character(df$reason.invalid)
  }
  df
})

# fix type of reason.terminated columns
ecitt_data <- lapply(ecitt_data, function(df) {
  if ("reason.terminated" %in% names(df)) {
    df$reason.terminated <- as.character(df$reason.terminated)
  }
  df
})

# bind rows
ecitt_data <- bind_rows(ecitt_data)

# id column

ecitt_data <- ecitt_data %>%
  mutate(id = as.factor(as.character(partRef))) %>% 
  drop_na(id)

# save compiled raw data
write.csv(ecitt_data,
          here("data", "ecitt_compiled_data.csv"),
          row.names = FALSE,)

# load demographic data

demographic_data <- read.csv(here("data", "ECITT data SPSS oct. 28.csv"))

demographic_data <- demographic_data %>% 
  rename(c("id" = "ID")) %>% 
  mutate(id = as.factor(id))

n_distinct(demographic_data$id)

n_distinct(ecitt_data$id)

compiled_ids <- ecitt_data %>% distinct(id)
demo_ids     <- demographic_data %>% distinct(id)

compiled_missing_in_demo <- anti_join(compiled_ids, demo_ids, by = "id")
compiled_missing_in_demo

demo_missing_in_compiled <- anti_join(demo_ids, compiled_ids, by = "id")
demo_missing_in_compiled

# merge demographic and task data

ecitt_data_merged <- full_join(ecitt_data,
                    demographic_data,
                    by = "id")

# missing data analysis

flexplot(1 ~ respTime,
         data = ecitt_data_merged)

# fix final dataset
ecitt_data_merged_final <-  ecitt_data_merged %>% 
  mutate(
    category25 = factor(category25, levels = c(1,2),
                        labels = c("uni", "bi"))
  ) %>% 
  mutate(
    reason.invalid = as.factor(reason.invalid)
  ) %>% 
  relocate(id, 
           .before = tester) %>% 
  naniar::replace_with_na_all(condition = ~. == "999") %>% 
  mutate(
    trialVar = case_when(
      trialVar == "prpt" ~ "prpt",
      trialVar == "inhb" ~ "inhb"
    )
  )
  
# save merged data

write.csv(ecitt_data_merged_final,
          here("data", "ecitt_data_merged.csv"),
          row.names = FALSE)












































































