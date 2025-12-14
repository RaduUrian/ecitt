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
  lme4,
  sjPlot,
  glmmTMB,
  performance,
  emmeans,
  marginaleffects,
  clubSandwich
)

require(flexplot)

renv::snapshot()

options(scipen = 999) # Remove scientific notations

options(es.use_symbols = TRUE)

ecitt_data_acc <- read.csv(here("data", "ecitt_data_acc_final.csv"))

ecitt_data_rt <- read.csv(here("data", "ecitt_data_rt_final.csv"))

# gender as factor

ecitt_data_acc$gender.y <- as.factor(ecitt_data_acc$gender.y)
ecitt_data_rt$gender.y <- as.factor(ecitt_data_rt$gender.y)

##########################################################
# GLMM for accuracy (cont. nondomlang)
##########################################################
acc_main <- glmer(accu ~ age_c + gender.y + income + nondom_c + trialVar + trialNo + (1 | id),
                  data = ecitt_data_acc,
                  family = binomial(link = "logit"),
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 0
)

summary(acc_main)

tab_model(acc_main,
          digits = 3)

performance::check_model(acc_main)

# average marginal means follow-up

avg_slopes(acc_main,
           variables = "age_c",
           type = "response",
           re.form = NA)

avg_predictions(
  acc_main,
  by = "gender.y",
  type = "response",
  re.form = NA
)

avg_comparisons(
  acc_main,
  variables = "gender.y"
)

avg_slopes(acc_main,
           variables = "income",
           type = "response",
           re.form = NA)

avg_slopes(acc_main,
           variables = "nondom_c",
           type = "response",
           re.form = NA)

avg_predictions(
  acc_main,
  by = "trialVar",
  type = "response",
  re.form = NA
)

avg_comparisons(
  acc_main,
  variables = "trialVar"
)

avg_slopes(acc_main,
           variables = "trialNo",
           type = "response",
           re.form = NA)



##########################################################
# GLMM for accuracy (categorical nondomlang)
##########################################################
acc_main_cat <- glmer(accu ~ age.y + gender.y + income + category25 + trialVar + (1 | id),
                  data = ecitt_data_acc,
                  family = binomial(link = "logit"),
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 0
)

summary(acc_main_cat)

tab_model(acc_main_cat)

performance::check_model(acc_main_cat)

##########################################################
# GLMM for RT (continuous nondomlang)
##########################################################
flexplot(1 ~ respTime, ecitt_data_rt)

rt_main_gamma <- glmmTMB(respTime ~ age.y + gender.y + income + nondomlanguagepercent + trialVar + trialNo + (1 | id),
                  data = ecitt_data_rt,
                  family = Gamma(link = "log")
)

performance::check_model(rt_main_gamma)

rt_main_gamma_var <- glmmTMB(respTime ~ age.y + gender.y + income + nondomlanguagepercent + trialVar + trialNo + (1 | id),
                         data = ecitt_data_rt,
                         dispformula = ~ trialVar,
                         family = Gamma(link = "log")
)

performance::check_model(rt_main_gamma_var)

rt_main_log_trans <- lmer(log(respTime) ~ age.y + gender.y + income + nondomlanguagepercent + trialVar + trialNo + (1 | id),
                         data = ecitt_data_rt
)

performance::check_model(rt_main_log_trans)

summary(rt_main_log_trans)

tab_model(rt_main_log_trans,
          digits = 3)

# exponentiate coefficients for multiplicative effect

# Fixed-effect coefficients on log scale
fe <- fixef(rt_main_log_trans)

# 95% CIs for fixed effects on log scale
ci <- confint(rt_main_log_trans, parm = names(fe), method = "Wald")
# If profile CIs converge, you can use:
# ci <- confint(rt_main_log_trans, parm = names(fe))  # profile

# Combine and exponentiate to get multiplicative effects
mult_effects <- data.frame(
  term = names(fe),
  estimate_log = unname(fe),
  ci_lower_log = ci[, 1],
  ci_upper_log = ci[, 2]
) %>%
  mutate(
    estimate_mult = 100*(exp(estimate_log)-1),
    ci_lower_mult = 100*(exp(ci_lower_log)-1),
    ci_upper_mult = 100*(exp(ci_upper_log)-1)
  ) %>% 
  as.data.frame()

mult_effects

# average marginal effects

avg_predictions(
  rt_main_log_trans,
  by = "trialVar",
  transform = "exp",
  re.form = NA
)

##########################################################
# GLMM for RT (binary nondomlang)
##########################################################

rt_main_cat_log_trans <- lmer(log(respTime) ~ age.y + gender.y + income + category25 + trialVar + trialNo + (1 | id),
                          data = ecitt_data_rt
)

performance::check_model(rt_main_cat_log_trans)

tab_model(rt_main_cat_log_trans)

# citations

citation("marginaleffects")


























































