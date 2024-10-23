# ==============================================================================
# Purpose: Run dropout simulations with parameters calculated in step 2 of this pipeline.

# Contact: sam.rosenbaum@umontana.edu
# ==============================================================================
# Load libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
# ==============================================================================
# Store simulation parameters as objects:
nc_adults <- 3769
samp_adults_MS <- 295
samp_adults_TRIB <- 306
offspring <- 480000
samp_offspring <- 682
adult_ages <- c("1.1", "1.2", "1.3", "1.4")
age_structure <- c(0.18, 0.12, 0.60, 0.10)
MS_selex <- c(0.18, 0.12, 0.60, 0.10)
TRIB_selex <- c(0.01, 0.34, 0.29, 0.36)
RS_prob <- c(0.07, 0.17, 0.28, 0.38)
# ==============================================================================
# Load function "tGMRsim" from step 4 and begin simulations:
# ==============================================================================
# Mainstem 0% dropout
mainstem_0_drop <-tGMRsim(Nc_adult = nc_adults, 
                   n_adult = samp_adults_MS, 
                   Nc_offspring = offspring, 
                   n_offspring = samp_offspring, 
                   ages = adult_ages, 
                   age_prob = age_structure, 
                   select_prob = MS_selex, 
                   rs_prob = RS_prob, 
                   drop_prob = c(0.05, 0.05, 0.05, 0.05),
                   drop_percent = 0.00,
                   iterations = 1000,
                   scenario = 1,
                   delta_selex = 0,
                   delta_rs = 0)
# ==============================================================================
# Mainstem 0% dropout (but with **tributary selectivity**)
mainstem_0_drop_tribselex <- tGMRsim(Nc_adult = nc_adults, 
                           n_adult = samp_adults_MS, 
                           Nc_offspring = offspring, 
                           n_offspring = samp_offspring, 
                           ages = adult_ages, 
                           age_prob = age_structure, 
                           select_prob = TRIB_selex, 
                           rs_prob = RS_prob, 
                           drop_prob = c(0.05, 0.05, 0.05, 0.05),
                           drop_percent = 0.00,
                           iterations = 1000,
                           scenario = 1,
                           delta_selex = 0,
                           delta_rs = 0)
# ==============================================================================
# Tributary 0% dropout
tributary_0_drop <- tGMRsim(Nc_adult = nc_adults, 
                   n_adult = samp_adults_TRIB, 
                   Nc_offspring = offspring, 
                   n_offspring = samp_offspring, 
                   ages = adult_ages, 
                   age_prob = age_structure, 
                   select_prob = TRIB_selex, 
                   rs_prob = RS_prob, 
                   drop_prob = c(0.05, 0.05, 0.05, 0.05),
                   drop_percent = 0.00,
                   iterations = 1000,
                   scenario = 2,
                   delta_selex = 0,
                   delta_rs = 0)
# ==============================================================================
# ==============================================================================
# Tributary 5% dropout
tributary_5_drop <- tGMRsim(Nc_adult = nc_adults, 
                            n_adult = samp_adults_TRIB, 
                            Nc_offspring = offspring, 
                            n_offspring = samp_offspring, 
                            ages = adult_ages, 
                            age_prob = age_structure, 
                            select_prob = TRIB_selex, 
                            rs_prob = RS_prob, 
                            drop_prob = c(0.05, 0.05, 0.05, 0.05),
                            drop_percent = 0.05,
                            iterations = 1000,
                            scenario = 3,
                            delta_selex = 0,
                            delta_rs = 0)
# ==============================================================================
# ==============================================================================
# Tributary 10% dropout
tributary_10_drop <- tGMRsim(Nc_adult = nc_adults, 
                             n_adult = samp_adults_TRIB, 
                             Nc_offspring = offspring, 
                             n_offspring = samp_offspring, 
                             ages = adult_ages, 
                             age_prob = age_structure, 
                             select_prob = TRIB_selex, 
                             rs_prob = RS_prob, 
                             drop_prob = c(0.05, 0.05, 0.05, 0.05),
                             drop_percent = 0.10,
                             iterations = 1000,
                             scenario = 4,
                             delta_selex = 0,
                             delta_rs = 0)
# ==============================================================================
# ==============================================================================
# Tributary 15% dropout
tributary_15_drop <- tGMRsim(Nc_adult = nc_adults, 
                             n_adult = samp_adults_TRIB, 
                             Nc_offspring = offspring, 
                             n_offspring = samp_offspring, 
                             ages = adult_ages, 
                             age_prob = age_structure, 
                             select_prob = TRIB_selex, 
                             rs_prob = RS_prob, 
                             drop_prob = c(0.05, 0.05, 0.05, 0.05),
                             drop_percent = 0.15,
                             iterations = 1000,
                             scenario = 5,
                             delta_selex = 0,
                             delta_rs = 0)
# ==============================================================================
# Tributary 20% dropout
tributary_20_drop <- tGMRsim(Nc_adult = nc_adults, 
                             n_adult = samp_adults_TRIB, 
                             Nc_offspring = offspring, 
                             n_offspring = samp_offspring, 
                             ages = adult_ages, 
                             age_prob = age_structure, 
                             select_prob = TRIB_selex, 
                             rs_prob = RS_prob, 
                             drop_prob = c(0.05, 0.05, 0.05, 0.05),
                             drop_percent = 0.20,
                             iterations = 1000,
                             scenario = 6,
                             delta_selex = 0,
                             delta_rs = 0)
# ==============================================================================
# Tributary 25% dropout
tributary_25_drop <- tGMRsim(Nc_adult = nc_adults, 
                             n_adult = samp_adults_TRIB, 
                             Nc_offspring = offspring, 
                             n_offspring = samp_offspring, 
                             ages = adult_ages, 
                             age_prob = age_structure, 
                             select_prob = TRIB_selex, 
                             rs_prob = RS_prob, 
                             drop_prob = c(0.05, 0.05, 0.05, 0.05),
                             drop_percent = 0.25,
                             iterations = 1000,
                             scenario = 7,
                             delta_selex = 0,
                             delta_rs = 0)
# ==============================================================================
# Tributary 30% dropout
tributary_30_drop <- tGMRsim(Nc_adult = nc_adults, 
                             n_adult = samp_adults_TRIB, 
                             Nc_offspring = offspring, 
                             n_offspring = samp_offspring, 
                             ages = adult_ages, 
                             age_prob = age_structure, 
                             select_prob = TRIB_selex, 
                             rs_prob = RS_prob, 
                             drop_prob = c(0.05, 0.05, 0.05, 0.05),
                             drop_percent = 0.30,
                             iterations = 1000,
                             scenario = 8,
                             delta_selex = 0,
                             delta_rs = 0)
# ==============================================================================
# Tributary 35% dropout
tributary_35_drop <- tGMRsim(Nc_adult = nc_adults, 
                             n_adult = samp_adults_TRIB, 
                             Nc_offspring = offspring, 
                             n_offspring = samp_offspring, 
                             ages = adult_ages, 
                             age_prob = age_structure, 
                             select_prob = TRIB_selex, 
                             rs_prob = RS_prob, 
                             drop_prob = c(0.05, 0.05, 0.05, 0.05),
                             drop_percent = 0.35,
                             iterations = 1000,
                             scenario = 9,
                             delta_selex = 0,
                             delta_rs = 0)
# ==============================================================================
# Tributary 40% dropout
tributary_40_drop <- tGMRsim(Nc_adult = nc_adults, 
                             n_adult = samp_adults_TRIB, 
                             Nc_offspring = offspring, 
                             n_offspring = samp_offspring, 
                             ages = adult_ages, 
                             age_prob = age_structure, 
                             select_prob = TRIB_selex, 
                             rs_prob = RS_prob, 
                             drop_prob = c(0.05, 0.05, 0.05, 0.05),
                             drop_percent = 0.40,
                             iterations = 1000,
                             scenario = 10,
                             delta_selex = 0,
                             delta_rs = 0)
# ==============================================================================
# Visualize results
# ==============================================================================
# ==============================================================================
# Create a dataframe to hold all simulation results
ms_tr_dat <- as.data.frame(rbind(mainstem_0_drop, tributary_0_drop, tributary_5_drop, tributary_10_drop, tributary_15_drop, tributary_20_drop, tributary_25_drop, tributary_30_drop, tributary_35_drop, tributary_40_drop))

drop_tr_dat <- as.data.frame(rbind(tributary_0_drop, tributary_5_drop, tributary_10_drop, tributary_15_drop, tributary_20_drop, tributary_25_drop, tributary_30_drop, tributary_35_drop, tributary_40_drop))

# tributary dropout simulations (compared to empirical hypergeometric estimate)
ggplot(data = drop_tr_dat, aes(x = scenario, y = mean)) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_confint, ymax = upper_confint, width = 0.1)) +
  ylab("Simulated tGMR Estimate") +
  xlab("Dropout Scenarios") +
  ylim(0, 5000) + 
  geom_hline(yintercept=2284, linetype = "dashed", color = "red", linewidth = 0.3) +
scale_x_continuous(labels = c('0%', '5%', '10%', '15%', '20%', '25%', '30%', "35%", '40%'), breaks = c(2:10)) +
  theme_classic() +
  theme(text=element_text(size=16))

# ==============================================================================
# tributary dropout simulations (compared to empirical binomial estimate)
ggplot(data = drop_tr_dat, aes(x = scenario, y = mean)) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_confint, ymax = upper_confint, width = 0.1)) +
  ylab("Simulated tGMR Estimate") +
  xlab("Dropout Scenarios") +
  ylim(0, 5000) + 
  geom_hline(yintercept=2516, linetype = "dashed", color = "red", linewidth = 0.3) +
  scale_x_continuous(labels = c('0%', '5%', '10%', '15%', '20%', '25%', '30%', "35%", '40%'), breaks = c(2:10)) +
  theme_classic() 
# ==============================================================================

str(ms_tr_dat)
# Plot only mainstem vs tributary sims w/o dropout
ms_tr_nodrop <- ms_tr_dat %>%
  filter(scenario %in% c(1,2)) %>%
  ggplot(aes(x = factor(scenario), y = mean)) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_confint, ymax = upper_confint, width = 0.1)) +
  ylab("tGMR Estimate") +
  xlab("Simulations with 0% Dropout") +
  ylim(0, 5500) + 
  scale_x_discrete(labels = c("Mainstem", "Tributary")) +
  theme_classic() +
  theme(text=element_text(size=14))

# Plot only mainstem vs tributary sims WITH 30% dropout
ms_tr_30drop <- ms_tr_dat %>%
  filter(scenario %in% c(1, 8)) %>%
  ggplot(aes(x = factor(scenario), y = mean)) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_confint, ymax = upper_confint, width = 0.1)) +
  ylab("") +
  xlab("Simulations with 30% Dropout") +
  scale_x_discrete(labels = c("Mainstem", "Tributary")) +
  ylim(0, 5500) + 
  theme_classic() +
  theme(text=element_text(size=14))

# Plot only mainstem vs tributary sims WITH 25% dropout
ms_tr_25drop <- ms_tr_dat %>%
  filter(scenario %in% c(1, 7)) %>%
  ggplot(aes(x = factor(scenario), y = mean)) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_confint, ymax = upper_confint, width = 0.1)) +
  ylab("") +
  xlab("Simulations with 25% Dropout") +
  scale_x_discrete(labels = c("Mainstem", "Tributary")) +
  ylim(0, 5500) + 
  theme_classic() +
  theme(text=element_text(size=14))
# ==============================================================================
# Plot empirical hypergeometric tGMR mainstem vs tributary 
est_hyp <- c("E1_hyp", "E2_hyp")
estim_hyp <- c(3665, 2284)
low_hyp <- c(2852, 1936)
upp_hyp <- c(4478, 2632)

tGMR_hyp <- as.data.frame(cbind(est_hyp, estim_hyp, low_hyp, upp_hyp))

tGMR_hyp$estim_hyp <- as.integer(tGMR_hyp$estim_hyp)
tGMR_hyp$low_hyp <- as.integer(tGMR_hyp$low_hyp)
tGMR_hyp$upp_hyp <- as.integer(tGMR_hyp$upp_hyp)

ms_tr_empirical <- ggplot(data = tGMR_hyp, aes(est_hyp, estim_hyp)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = low_hyp, ymax = upp_hyp, width = 0.1)) +
  ylab("") +
  xlab("Empirical tGMR") +
  scale_x_discrete(labels=c("Mainstem", "Tributary")) +
  ylim(0, 5500) +
  theme_classic() +
  theme(text=element_text(size=14))
# ==============================================================================
# Plot empirical binomial tGMR mainstem vs tributary 
est_bin <- c("E1_bin", "E2_bin")
estim_bin <- c(4794, 2516)
low_bin <- c(3806, 2159)
upp_bin <- c(5781, 2874)

tGMR_bin <- as.data.frame(cbind(est_bin, estim_bin, low_bin, upp_bin))

tGMR_bin$estim_bin <- as.integer(tGMR_bin$estim_bin)
tGMR_bin$low_bin <- as.integer(tGMR_bin$low_bin)
tGMR_bin$upp_bin <- as.integer(tGMR_bin$upp_bin)

bin_empirical <- ggplot(data = tGMR_bin, aes(est_bin, estim_bin)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = low_bin, ymax = upp_bin, width = 0.1)) +
  ylab("") +
  xlab("Empirical binomial tGMR") +
  scale_x_discrete(labels=c("Mainstem", "Tributary")) +
  ylim(0, 6000) +
  theme_classic() +
  theme(text=element_text(size=16))
# ==============================================================================
# Compare all 3
ggarrange(ms_tr_nodrop, ms_tr_30drop,  ms_tr_empirical, nrow = 1, common.legend = TRUE)

# Compare 0% dropout
drop_0_fig <- ggarrange(ms_tr_empirical, ms_tr_nodrop, nrow = 1, common.legend = TRUE)

# Compare 30% dropout
drop_30_fig <- ggarrange(ms_tr_empirical, ms_tr_30drop, nrow = 1, common.legend = TRUE)
