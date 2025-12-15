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

ecitt_data_rt <- read.csv(here("data", "ecitt_data_rt_final.csv"))

# power_FLP <- mixedpower(model = rt_demog, data = ecitt_data_rt_cc,
#                         fixed_effects = c("age.y", "gender.y", "nondomlanguagepercent", "income", "trialNo", "trialVar"),
#                         simvar = "id", steps = c(10, 20, 30, 40, 50, 60),
#                         critical_value = 2, n_sim = 10)
# 
# multiplotPower(power_FLP,
#                filename = here("output", "property_power_original_sesoi.png"))

########### wrangling ######################

# gender as factor

ecitt_data_acc$gender.y <- as.factor(ecitt_data_acc$gender.y)
ecitt_data_rt$gender.y <- as.factor(ecitt_data_rt$gender.y)

ecitt_data_rt2 <- ecitt_data_rt %>%
  mutate(
    logRT = log(respTime)
  )

ecitt_data_rt_cc <- ecitt_data_rt2 %>%
  filter(complete.cases(.))

# Standardize continuous L2 predictors (recommended for sensitivity)
# ecitt_data_rt <- ecitt_data_rt %>%
#   mutate(
#     age_z    = as.numeric(scale(age.y)),
#     income_z = as.numeric(scale(income)),                   # you said income is numeric (2–7)
#     nondom_z = as.numeric(scale(nondomlanguagepercent))
#   )
# 
# ecitt_data_acc <- ecitt_data_acc %>%
#   mutate(
#     age_z    = as.numeric(scale(age_c)),                    # or scale(age.y) if you refit with age.y
#     income_z = as.numeric(scale(income)),
#     nondom_z = as.numeric(scale(nondom_c))                  # or scale(raw percent) if you refit
#   )

##################### refit models ###############
rt_demog <- lmer(logRT ~ age.y + gender.y + income + nondomlanguagepercent + trialVar + trialNo + (1 | id),
                 data = ecitt_data_rt_cc
)

acc_demog <- glmer(accu ~ age_c + gender.y + income + nondom_c + trialVar + trialNo + (1 | id),
                   data = ecitt_data_acc,
                   family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "bobyqa"),
                   nAGQ = 0
)

################ RT Model sensitivity ################

####################### age

# Define effect sizes: -1%, -2%, -3%, -4%, -5% change per unit age
targets <- log(seq(0.99, 0.90, by = -0.01))  # Test -1% to -10% in 1% increments

# Run power simulations across different effect sizes
res_rt_age <- do.call(
  rbind,
  lapply(targets, function(b) {
    
    # Copy model (avoid side effects)
    mod_b <- rt_demog
    mod_b@beta[match("age.y", names(fixef(rt_demog)))] <- b
    
    # Run power simulation
    ps <- powerSim(
      mod_b,
      fixed("age.y", "t"),
      nsim = 200
    )
    
    # Extract power summary
    ps_summary <- summary(ps)
    
    # Extract results - powerSim returns percentage values
    data.frame(
      beta  = b,
      mult  = exp(b),
      pct   = (exp(b) - 1) * 100,
      power = ps_summary$mean,     # Convert percentage to proportion
      lower = ps_summary$lower,
      upper = ps_summary$upper
    )
  })
)

# View results
print(res_rt_age)

png(here("output","power_curve_age_rt.png"), width = 1000, height = 800)

# Create power curve plot
plot(
  -res_rt_age$pct, res_rt_age$power,
  pch = 16,
  cex = 1.5,
  ylim = c(0, 1),
  xlim = c(0, max(-res_rt_age$pct) * 1.1),
  xlab = "% faster RT per month",
  ylab = "Power",
  main = "Sensitivity Analysis: Age Effect on RT"
)

# Add confidence intervals
segments(-res_rt_age$pct, res_rt_age$lower, 
         -res_rt_age$pct, res_rt_age$upper,
         lwd = 2)

# Add reference line for 80% power
abline(h = 0.80, lty = 2, col = "red", lwd = 2)
grid()

# Optional: Add power values as text labels
text(-res_rt_age$pct, res_rt_age$power, 
     labels = paste0(round(res_rt_age$power * 100, 1), "%"),
     pos = 3, cex = 0.8, offset = 0.5)

dev.off()

####################### income

# Define effect sizes: -1%, -2%, -3%, -4%, -5% change per unit age
targets <- log(seq(.99, .90, by = -.01))  # Test -1% to -10% in 1% increments

# Run power simulations across different effect sizes
res_rt_income <- do.call(
  rbind,
  lapply(targets, function(b) {
    
    # Copy model (avoid side effects)
    mod_b <- rt_demog
    mod_b@beta[match("income", names(fixef(rt_demog)))] <- b
    
    # Run power simulation
    ps <- powerSim(
      mod_b,
      fixed("income", "t"),
      nsim = 200
    )
    
    # Extract power summary
    ps_summary <- summary(ps)
    
    # Extract results - powerSim returns percentage values
    data.frame(
      beta  = b,
      mult  = exp(b),
      pct   = (exp(b) - 1) * 100,
      power = ps_summary$mean,     # Convert percentage to proportion
      lower = ps_summary$lower,
      upper = ps_summary$upper
    )
  })
)

# View results
print(res_rt_income)

png(here("output","power_curve_income_rt.png"), width = 1000, height = 800)

# Create power curve plot
plot(
  -res_rt_income$pct, res_rt_income$power,
  pch = 16,
  cex = 1.5,
  ylim = c(0, 1),
  xlim = c(0, max(-res_rt_income$pct) * 1.1),
  xlab = "% faster RT per income bracket",
  ylab = "Power",
  main = "Sensitivity Analysis: Income Effect on RT"
)

# Add confidence intervals
segments(-res_rt_income$pct, res_rt_income$lower, 
         -res_rt_income$pct, res_rt_income$upper,
         lwd = 2)

# Add reference line for 80% power
abline(h = 0.80, lty = 2, col = "red", lwd = 2)
grid()

# Optional: Add power values as text labels
text(-res_rt_income$pct, res_rt_income$power, 
     labels = paste0(round(res_rt_income$power * 100, 1), "%"),
     pos = 3, cex = 0.8, offset = 0.5)

dev.off()

####################### nondomlanguage

# Define effect sizes: -1%, -2%, -3%, -4%, -5% change per unit age
targets_nondonlang <- log(seq(.999, .989, by = -.001))  # Test -1% to -10% in 1% increments

# Run power simulations across different effect sizes
res_rt_nondomlanguagepercent <- do.call(
  rbind,
  lapply(targets, function(b) {
    
    # Copy model (avoid side effects)
    mod_b <- rt_demog
    mod_b@beta[match("nondomlanguagepercent", names(fixef(rt_demog)))] <- b
    
    # Run power simulation
    ps <- powerSim(
      mod_b,
      fixed("nondomlanguagepercent", "t"),
      nsim = 200
    )
    
    # Extract power summary
    ps_summary <- summary(ps)
    
    # Extract results - powerSim returns percentage values
    data.frame(
      beta  = b,
      mult  = exp(b),
      pct   = (exp(b) - 1) * 100,
      power = ps_summary$mean,     # Convert percentage to proportion
      lower = ps_summary$lower,
      upper = ps_summary$upper
    )
  })
)

# View results
print(res_rt_nondomlanguagepercent)

png(here("output","power_curve_nondomlanguagepercent_rt.png"), width = 1000, height = 800)

# Create power curve plot
plot(
  -res_rt_nondomlanguagepercent$pct, res_rt_nondomlanguagepercent$power,
  pch = 16,
  cex = 1.5,
  ylim = c(0, 1),
  xlim = c(0, max(-res_rt_nondomlanguagepercent$pct) * 1.1),
  xlab = "% faster RT per 1% non-dominant language exposure",
  ylab = "Power",
  main = "Sensitivity Analysis: non-dominant language exposure effect on RT"
)

# Add confidence intervals
segments(-res_rt_nondomlanguagepercent$pct, res_rt_nondomlanguagepercent$lower, 
         -res_rt_nondomlanguagepercent$pct, res_rt_nondomlanguagepercent$upper,
         lwd = 2)

# Add reference line for 80% power
abline(h = 0.80, lty = 2, col = "red", lwd = 2)
grid()

# Optional: Add power values as text labels
text(-res_rt_nondomlanguagepercent$pct, res_rt_nondomlanguagepercent$power, 
     labels = paste0(round(res_rt_nondomlanguagepercent$power * 100, 1), "%"),
     pos = 3, cex = 0.8, offset = 0.5)

dev.off()

################ Gender sensitivity analysis ################

# First, check what levels you have and which is the reference
levels(ecitt_data_rt_cc$gender.y)
# Check current coefficients
fixef(rt_demog)

# Define effect sizes for gender comparisons
# Example: testing different effect sizes for one gender contrast
# Let's say gender.y has levels: "Male" (reference), "Female", "Other"
# You'll have coefficients like gender.yFemale and gender.yOther

# Option 1: Test power for ONE specific contrast (e.g., Female vs Male)
targets_gender2 <- log(seq(.99, .80, by = -.02))

res_gender2 <- do.call(
  rbind,
  lapply(targets_gender2, function(b) {
    
    mod_b <- rt_demog
    mod_b@beta[match("gender.y2", names(fixef(rt_demog)))] <- b
    
    ps <- powerSim(
      mod_b,
      fixed("gender.y2", "t"),
      nsim = 200
    )
    
    ps_summary <- summary(ps)
    
    data.frame(
      contrast = "Gender 2 vs Gender 1",
      beta  = b,
      mult  = exp(b),
      pct   = (exp(b) - 1) * 100,
      power = ps_summary$mean,
      lower = ps_summary$lower,
      upper = ps_summary$upper
    )
  })
)

png(here("output","power_curve_gender_rt.png"), width = 1000, height = 800)

# Plot 1: Power curve for gender.y2
plot(
  -res_gender2$pct, res_gender2$power,
  pch = 16,
  cex = 1.5,
  ylim = c(0, 1),
  xlim = range(-res_gender2$pct) * c(0.9, 1.1),
  xlab = "% faster RT for girls vs boys",
  ylab = "Power",
  main = "Power Analysis: girls vs. boys"
)
segments(-res_gender2$pct, res_gender2$lower, 
         -res_gender2$pct, res_gender2$upper, lwd = 2)
abline(h = 0.80, lty = 2, col = "red", lwd = 2)
grid()

dev.off()

####################### trial number

# Define effect sizes: -1%, -2%, -3%, -4%, -5% change per unit age
targets_trialNo <- log(seq(.999, .989, by = -.001))  # Test -1% to -10% in 1% increments

# Run power simulations across different effect sizes
res_rt_trialNo <- do.call(
  rbind,
  lapply(targets_trialNo, function(b) {
    
    # Copy model (avoid side effects)
    mod_b <- rt_demog
    mod_b@beta[match("trialNo", names(fixef(rt_demog)))] <- b
    
    # Run power simulation
    ps <- powerSim(
      mod_b,
      fixed("trialNo", "t"),
      nsim = 200
    )
    
    # Extract power summary
    ps_summary <- summary(ps)
    
    # Extract results - powerSim returns percentage values
    data.frame(
      beta  = b,
      mult  = exp(b),
      pct   = (exp(b) - 1) * 100,
      power = ps_summary$mean,     # Convert percentage to proportion
      lower = ps_summary$lower,
      upper = ps_summary$upper
    )
  })
)

# View results
print(res_rt_trialNo)

png(here("output","power_curve_trialNo_rt.png"), width = 1000, height = 800)

# Create power curve plot
plot(
  -res_rt_trialNo$pct, res_rt_trialNo$power,
  pch = 16,
  cex = 1.5,
  ylim = c(0, 1),
  xlim = c(0, max(-res_rt_trialNo$pct) * 1.1),
  xlab = "% faster RT per additional trial",
  ylab = "Power",
  main = "Sensitivity Analysis: trial number effect on RT"
)

# Add confidence intervals
segments(-res_rt_trialNo$pct, res_rt_trialNo$lower, 
         -res_rt_trialNo$pct, res_rt_trialNo$upper,
         lwd = 2)

# Add reference line for 80% power
abline(h = 0.80, lty = 2, col = "red", lwd = 2)
grid()

# Optional: Add power values as text labels
text(-res_rt_trialNo$pct, res_rt_trialNo$power, 
     labels = paste0(round(res_rt_trialNo$power * 100, 1), "%"),
     pos = 3, cex = 0.8, offset = 0.5)

dev.off()

####################### trial type

# First, check the actual coefficient name in your model
fixef(rt_demog)
# You should see something like "trialVarprpt" (not just "trialVar")

# Check the levels
levels(as.factor(ecitt_data_rt_cc$trialVar))
# Based on your fixef output earlier, it looks like 'inhb' is reference and 'prpt' is the contrast

# Define effect sizes: -1% to -10% change
targets_trialVar <- log(seq(.99, .89, by = -.01))

# Run power simulations across different effect sizes
res_rt_trialVar <- do.call(
  rbind,
  lapply(targets_trialVar, function(b) {
    
    # Copy model (avoid side effects)
    mod_b <- rt_demog
    # Use the EXACT coefficient name from fixef(rt_demog)
    mod_b@beta[match("trialVarprpt", names(fixef(rt_demog)))] <- b
    
    # Run power simulation
    ps <- powerSim(
      mod_b,
      fixed("trialVarprpt", "t"),  # Use exact coefficient name
      nsim = 200  # Increased from 10 for more stable estimates
    )
    
    # Extract power summary
    ps_summary <- summary(ps)
    
    # Extract results
    data.frame(
      beta  = b,
      mult  = exp(b),
      pct   = (exp(b) - 1) * 100,
      power = ps_summary$mean,
      lower = ps_summary$lower,
      upper = ps_summary$upper
    )
  })
)

# View results
print(res_rt_trialVar)

# Create power curve plot
png(here("output","power_curve_trialVar_rt.png"), width = 1000, height = 800)

plot(
  -res_rt_trialVar$pct, res_rt_trialVar$power,
  pch = 16,
  cex = 1.5,
  ylim = c(0, 1),
  xlim = c(0, max(-res_rt_trialVar$pct) * 1.1),
  xlab = "% faster RT for prepotent vs. inhibitory trials",
  ylab = "Power",
  main = "Sensitivity Analysis: Trial Type Effect on RT"
)

# Add confidence intervals
segments(-res_rt_trialVar$pct, res_rt_trialVar$lower, 
         -res_rt_trialVar$pct, res_rt_trialVar$upper,
         lwd = 2)

# Add reference line for 80% power
abline(h = 0.80, lty = 2, col = "red", lwd = 2)
grid()

# Add power values as text labels
text(-res_rt_trialVar$pct, res_rt_trialVar$power, 
     labels = paste0(round(res_rt_trialVar$power * 100, 1), "%"),
     pos = 3, cex = 0.8, offset = 0.5)

dev.off()


























