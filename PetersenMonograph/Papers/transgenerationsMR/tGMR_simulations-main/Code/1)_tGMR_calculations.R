# Purpose: calculate tGMR estimates 

# Use pedigree results from https://github.com/krshedd/Chilkat-Chinook-tGMR

# Contact: sam.rosenbaum@umontana.edu
# ==============================================================================
# Load libraries
library(stringr)
library(tidyverse)
library(ggplot2)
library(dplyr)
# ==============================================================================

# Create binomial tGMR (Bailey modification) function: =========================
calc_binom <- function(n1, n2, m2) {
  N <- (n1*(n2+1)) / (m2+1)
  var <- (n1^2 * (n2+1) * (n2-m2)) / ((m2+1)^2 * (m2+2))
  cv <- sqrt(var)/N
  upper_binom <- N + (1.96 * (sqrt(var)))
  lower_binom <- N - (1.96 * (sqrt(var)))
  # Output list
  out <- NULL
  out$N <- N
  out$var <- var
  out$cv <- cv
  out$upper_binom <- upper_binom
  out$lower_binom <- lower_binom
  return(out)
}

# Create hypergeometric tGMR (Chapman modification) function: ==================
calc_hyper <- function(n1, n2, m2) {
  N <- ((n1+1)*(n2+1))/(m2+1) - 1
  var <- ((n1+1)*(n2+1)*(n1-m2)*(n2-m2)) / ((m2+1)^2 * (m2+2))
  cv <- sqrt(var)/N
  upper_hyper <- N + (1.96 * (sqrt(var)))
  lower_hyper <- N - (1.96 * (sqrt(var)))
  # Output list
  out <- NULL
  out$N <- N
  out$var <- var
  out$cv <- cv
  out$upper_hyper <- upper_hyper
  out$lower_hyper <- lower_hyper
  return(out)
}

# ==============================================================================
# ==============================================================================
# ==============================================================================

#  Calculate tGMR estimates using all adults sampled (event 1 + event 2):
# ==============================================================================
# First, calculate binomial tGMR estimate for all adults (event 1 + event 2).

# Use parameters from most recent COLONY output https://github.com/krshedd/Chilkat-Chinook-tGMR.

# n1 = 583
# n2 = 682 x 2 = 1,364
# m2 = 237

bin_combined <- calc_binom(n1 =583, n2 = 1364, m2 = 237)
str(bin_combined)

bin_combined_N <- round(bin_combined$N)
bin_combined_N # 3,344

bin_combined_lower95 <- round(bin_combined$lower_binom)
bin_combined_lower95 # 2,958

bin_combined_upper95 <- round(bin_combined$upper_binom)
bin_combined_upper95 # 3,729

bin_combined_CV <- bin_combined$cv
bin_combined_CV # 0.06

# Second, calculate hypergeometric for all adults (combined event 1 + event 2): 
# Use parameters from most recent COLONY output https://github.com/krshedd/Chilkat-Chinook-tGMR:

# n1 = 583
# n2 = 745
# m2 = 146

hyp_combined <- calc_hyper(n1 = 583, n2 = 745, m2 = 148)
str(hyp_combined)

hyp_combined_N <- round(hyp_combined$N)
hyp_combined_N # 2,923

hyp_combined_lower95 <- round(hyp_combined$lower_hyper)
hyp_combined_lower95 # 2,562

hyp_combined_upper95 <- round(hyp_combined$upper_hyper)
hyp_combined_upper95 # 3,284

hyp_combined_CV <- hyp_combined$cv
hyp_combined_CV # 0.06

# ==============================================================================
# ==============================================================================
# ==============================================================================

# Mainstem (event 1):
# ==============================================================================
# Calculate a binomial tGMR estimate for mainstem adults (event 1): ============
# Use parameters from most recent COLONY output https://github.com/krshedd/Chilkat-Chinook-tGMR

# n1 = 295
# n2 = 682 x 2 = 1,364
# m2 = 83

bin_mainstem <- calc_binom(n1 = 295, n2 = 1364, m2 = 83)
str(bin_mainstem)

bin_mainstem_N <- round(bin_mainstem$N)
bin_mainstem_N # 4,794

bin_mainstem_lower95 <- round(bin_mainstem$lower_binom)
bin_mainstem_lower95 # 3,806

bin_mainstem_upper95 <- round(bin_mainstem$upper_binom)
bin_mainstem_upper95 # 5,781

bin_mainstem_CV <- bin_mainstem$cv
bin_mainstem_CV # 0.11

# Calculate a hypergeometric tGMR estimate for mainstem adults (event 1): ======
# Use parameters from most recent COLONY output https://github.com/krshedd/Chilkat-Chinook-tGMR

# n1 = 295
# n2 = 705
# m2 = 56

hyp_mainstem <- calc_hyper(n1 = 295, n2 = 705, m2 = 56)
str(hyp_mainstem)

hyp_mainstem_N <- round(hyp_mainstem$N)
hyp_mainstem_N # 3,665

hyp_mainstem_lower95 <- round(hyp_mainstem$lower_hyper)
hyp_mainstem_lower95 # 2,852

hyp_mainstem_upper95 <- round(hyp_mainstem$upper_hyper)
hyp_mainstem_upper95 # 4,478

hyp_mainstem_CV <- hyp_mainstem$cv
hyp_mainstem_CV # 0.11
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Tributary (event 2):
# ==============================================================================
# Calculate a binomial tGMR estimate for tributary adults (event 2): ===========
# Use parameters from most recent COLONY output https://github.com/krshedd/Chilkat-Chinook-tGMR

# n1 = 306
# n2 = 682 x 2 = 1364
# m2 = 165

bin_tributary <- calc_binom(n1 =306, n2 = 1364, m2 = 165)
str(bin_tributary)

bin_tributary_N <- round(bin_tributary$N)
bin_tributary_N # 2,516

bin_tributary_lower95 <- round(bin_tributary$lower_binom)
bin_tributary_lower95 # 2,159

bin_tributary_upper95 <- round(bin_tributary$upper_binom)
bin_tributary_upper95 # 2,874

bin_tributary_CV <- bin_tributary$cv
bin_tributary_CV # 0.07

# Calculate a hypergeometric tGMR estimate for tributary adults (event 2): =====
# Use parameters from most recent COLONY output https://github.com/krshedd/Chilkat-Chinook-tGMR:

# n1 = 306
# n2 = 721
# m2 = 96

hyp_tributary <- calc_hyper(n1 = 306, n2 = 721, m2 = 96)
str(hyp_tributary)

hyp_tributary_N <- round(hyp_tributary$N)
hyp_tributary_N # 2,284

hyp_tributary_lower95 <- round(hyp_tributary$lower_hyper)
hyp_tributary_lower95 # 1,936

hyp_tributary_upper95 <- round(hyp_tributary$upper_hyper)
hyp_tributary_upper95 # 2,632

hyp_tributary_CV <- hyp_tributary$cv
hyp_tributary_CV # 0.08
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Calculate a traditional MR estimate for >= age 1.2 Chinook:
# ==============================================================================
# Use parameters from ADFG's mark-recapture data (contact brian.elliott1@alaska.gov for data):

# Nc = 3,769
# SE = 532
# CV = 532 / 3,769 = 0.1411515
# Upper 95% = 3,769 + (1.96 x 532) = 4,811.72
# Lower 95% = 3,769 - (1.96 x 532) = 2,726.28

traditional_N <- 3769
traditional_lower95 <- 2726
traditional_upper95 <- 4812
traditional_CV <- 0.14
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Make dataframes to hold estimates: ======================================
Estimator_Type <- c("Combined_Binomial", "Combined_Hypergeometric", "Mainstem_Binomial", "Mainstem_Hypergeometric", "Tributary_Binomial", "Tributary_Hypergeometric", "Traditional")

Point_Estimate <- c(bin_combined_N, hyp_combined_N, bin_mainstem_N, hyp_mainstem_N, bin_tributary_N, hyp_tributary_N, traditional_N)

Lower_95 <- c(bin_combined_lower95, hyp_combined_lower95, bin_mainstem_lower95, hyp_mainstem_lower95, bin_tributary_lower95, hyp_tributary_lower95, traditional_lower95)

Upper_95 <- c(bin_combined_upper95, hyp_combined_upper95, bin_mainstem_upper95, hyp_mainstem_upper95, bin_tributary_upper95, hyp_tributary_upper95, traditional_upper95)

CV <-c (bin_combined_CV, hyp_combined_CV, bin_mainstem_CV, hyp_mainstem_CV, bin_tributary_CV, hyp_tributary_CV, traditional_CV)

# Make a dataframe:
estimator_2020_dat <- as.data.frame(cbind(Estimator_Type, Point_Estimate, Lower_95, Upper_95, CV))

# Convert to integers:
estimator_2020_dat$Point_Estimate <- as.integer(estimator_2020_dat$Point_Estimate)
estimator_2020_dat$Lower_95 <- as.integer(estimator_2020_dat$Lower_95)
estimator_2020_dat$Upper_95 <- as.integer(estimator_2020_dat$Upper_95)
estimator_2020_dat$CV <- as.numeric(estimator_2020_dat$CV)

str(estimator_2020_dat)
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Subset dataframe to only include hypergeometric and traditional estimators:
hyp_estimator_dat <- filter(estimator_2020_dat, Estimator_Type %in% c("Combined_Hypergeometric", "Mainstem_Hypergeometric", "Tributary_Hypergeometric", "Traditional"))

# Subset dataframe to only include combined hypergeometric and traditional estimators:
hyp_combined_dat <- filter(estimator_2020_dat, Estimator_Type %in% c("Combined_Hypergeometric", "Traditional"))

# Subset dataframe to only include mainstem & tributary hypergeometric and traditional estimators:
hyp_main_trib_dat <- filter(estimator_2020_dat, Estimator_Type %in% c("Mainstem_Hypergeometric", "Tributary_Hypergeometric", "Traditional"))

# Subset dataframe to only include mainstem & tributary hypergeometric estimators:
hyp_dropout_dat <- filter(estimator_2020_dat, Estimator_Type %in% c("Mainstem_Hypergeometric", "Tributary_Hypergeometric"))

# Subset dataframe to only include combined hypergeometric estimators AND traditional estimator:
combo_dat <- filter(estimator_2020_dat, Estimator_Type %in% c("Traditional", "Combined_Hypergeometric", "Combined_Binomial"))
# ==============================================================================

# Visualize the data:
hyp_estimator_dat %>%
  mutate(Estimator_Type = fct_relevel(Estimator_Type, "Traditional", "Mainstem_Hypergeometric", "Tributary_Hypergeometric", "Combined_Hypergeometric")) %>%
  ggplot() +
  geom_point(aes(Estimator_Type, Point_Estimate), size = 1.5) +
  geom_errorbar(aes(Estimator_Type, Point_Estimate, ymin = Lower_95, ymax = Upper_95, width = 0.1)) +
  scale_x_discrete(labels=c("Traditional Mark-Recapture", "Mainstem tGMR", "Tributary tGMR", "Combined tGMR")) +
  ylab("Escapement") +
  xlab("Estimator") +
  scale_y_continuous(limits = c(0, 5000)) +
  theme_classic() 
# ==============================================================================