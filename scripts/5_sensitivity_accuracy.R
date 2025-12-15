########## SETUP ##########

rm(list = ls())

# Sets pseudo random number generator seed to fixed value for reproducibility
set.seed(1)

renv::restore()

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

if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)}

devtools::install_github("DejanDraschkow/mixedpower") 

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
  clubSandwich,
  simr,
  RLRsim,
  mixedpower
)

require(flexplot)

renv::snapshot()

options(scipen = 999) # Remove scientific notations

options(es.use_symbols = TRUE)

ecitt_data_acc <- read.csv(here("data", "ecitt_data_acc_final.csv"))

# power_FLP <- mixedpower(model = acc_demog, data = ecitt_data_rt_cc,
#                         fixed_effects = c("age.y", "gender.y", "nondomlanguagepercent", "income", "trialNo", "trialVar"),
#                         simvar = "id", steps = c(10, 20, 30, 40, 50, 60),
#                         critical_value = 2, n_sim = 10)
# 
# multiplotPower(power_FLP,
#                filename = here("output", "property_power_original_sesoi.png"))

########### wrangling ######################

# gender as factor

ecitt_data_acc$gender.y <- as.factor(ecitt_data_acc$gender.y)

ecitt_data_acc_cc <- ecitt_data_acc %>%
  filter(complete.cases(.))



##################### refit ###############

acc_demog <- glmer(accu ~ age_c + gender.y + income + nondom_c + trialVar + trialNo + (1 | id),
                   data = ecitt_data_acc_cc,
                   family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "bobyqa"),
                   nAGQ = 0
)

################ Accuracy Model Sensitivity Analysis ################

# Check your model coefficients first
fixef(acc_demog)
summary(acc_demog)

######################## AGE

# Get your observed effect
observed_beta <- fixef(acc_demog)["age_c"]
observed_OR <- exp(observed_beta)

cat("Observed age effect: OR =", round(observed_OR, 3), "\n")

# Create range around observed effect
# Test from 50% smaller to 150% larger than observed
multipliers <- seq(0.1, 1.5, by = 0.1)
targets_age <- observed_beta * multipliers

res_acc_age <- do.call(
  rbind,
  lapply(targets_age, function(b) {
    
    mod_b <- acc_demog
    mod_b@beta[match("age_c", names(fixef(acc_demog)))] <- b
    
    ps <- powerSim(
      mod_b,
      fixed("age_c", "z"),
      nsim = 200
    )
    
    ps_summary <- summary(ps)
    
    data.frame(
      beta  = b,
      OR    = exp(b),
      multiplier = b / observed_beta,  # How much larger/smaller than observed
      power = ps_summary$mean,
      lower = ps_summary$lower,
      upper = ps_summary$upper
    )
  })
)

# Plot relative to observed effect
png(here("output","power_curve_age_acc.png"), width = 1000, height = 800)
plot(
  res_acc_age$OR, res_acc_age$power,
  pch = 16, cex = 1.5,
  ylim = c(0, 1),
  xlim = range(res_acc_age$OR) * c(0.95, 1.05),  # Add xlim for better spacing
  xlab = "Odds Ratio per unit age",
  ylab = "Power",
  main = paste0("Sensitivity Analysis: Age Effect on Accuracy"),
  las = 1
)

# Add confidence intervals (fixed: removed negative sign)
segments(res_acc_age$OR, res_acc_age$lower, 
         res_acc_age$OR, res_acc_age$upper,
         lwd = 2)

# Add reference lines
abline(h = 0.80, lty = 2, col = "red", lwd = 2)
grid()
# Add power percentages (fixed: removed negative sign)
text(res_acc_age$OR, res_acc_age$power, 
     labels = paste0(round(res_acc_age$power * 100, 1), "%"),
     pos = 3, cex = 0.8, offset = 0.5)

dev.off()

print(res_acc_age)

######################## income

# Get your observed effect
observed_beta <- fixef(acc_demog)["income"]
observed_OR <- exp(observed_beta)

cat("Observed income effect: OR =", round(observed_OR, 3), "\n")

# Create range around observed effect
# Test from 50% smaller to 150% larger than observed
multipliers <- seq(0.1, 1.5, by = 0.1)
targets_income <- observed_beta * multipliers

res_acc_income <- do.call(
  rbind,
  lapply(targets_income, function(b) {
    
    mod_b <- acc_demog
    mod_b@beta[match("income", names(fixef(acc_demog)))] <- b
    
    ps <- powerSim(
      mod_b,
      fixed("income", "z"),
      nsim = 10
    )
    
    ps_summary <- summary(ps)
    
    data.frame(
      beta  = b,
      OR    = exp(b),
      multiplier = b / observed_beta,  # How much larger/smaller than observed
      power = ps_summary$mean,
      lower = ps_summary$lower,
      upper = ps_summary$upper
    )
  })
)

# Plot relative to observed effect
png(here("output","power_curve_income_acc.png"), width = 1000, height = 800)
plot(
  res_acc_income$OR, res_acc_income$power,
  pch = 16, cex = 1.5,
  ylim = c(0, 1),
  xlim = range(res_acc_income$OR) * c(0.95, 1.05),  # Add xlim for better spacing
  xlab = "Odds Ratio per unit income",
  ylab = "Power",
  main = paste0("Sensitivity Analysis: income Effect on Accuracy"),
  las = 1
)

# Add confidence intervals (fixed: removed negative sign)
segments(res_acc_income$OR, res_acc_income$lower, 
         res_acc_income$OR, res_acc_income$upper,
         lwd = 2)

# Add reference lines
abline(h = 0.80, lty = 2, col = "red", lwd = 2)
grid()
# Add power percentincomes (fixed: removed negative sign)
text(res_acc_income$OR, res_acc_income$power, 
     labels = paste0(round(res_acc_income$power * 100, 1), "%"),
     pos = 3, cex = 0.8, offset = 0.5)

dev.off()

print(res_acc_income)







































