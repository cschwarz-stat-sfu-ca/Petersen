# ==============================================================================
# Purpose: Create a function to simulate the influence of age-specific differences in reproductive success and selectivity on tGMR abundance estimation

# Creators: Sam Rosenbaum, Sam May, & Curry Cunningham

# Contact: sam.rosenbaum@umontana.edu
# ==============================================================================
# Load Library ===============================================================
library(tidyverse)
# ==============================================================================
# Below is code for a simple and flexible individual based model (IBM) that can simulate the accuracy and precision of trans-generational genetic mark-recapture (tGMR) under various demographic and sampling scenarios. Potential tGMR users are encouraged to harvest this IBM and adapt it for their specific applications. 

# The function *inputs* are: 
# 1) assumed adult census size present at the time of sampling *Nc_adult*, 
# 2) number of adults sampled *n_adults*, 
# 3) assumed number of offspring present at the time of sampling *Nc_offspring*, 
# 4) number of offspring sampled *n_offspring*, 
# 5) ages of the adults *adults* (as a vector), 
# 6) probability of an adult being a given age *age_prob* (as a vector), 
# 7) probability of an adult being sampled given its age *select_prob* (as a vector), 
# 8) probability of an adult being reproductively successful given its age *rs_prob* (as a vector), 
# 9) the probability of an adult dropping out of the system before sampling occurs given its age *drop_prob* (as a vector),
# 10) the percent of adults 'dropped' from the system before sampling *drop_percent* (as a vector), 
# 11) the number of iterations the IBM will run *iterations*, 
# 12) and a way to associate results with a given set of parameter values (scenario).

# Inputs *delta_selex* and *delta_rs* were added for step 5 of this pipeline, and can be ignored for most users.

# The key *outputs* are:
# 1) averaged tGMR estimate *mean*, 
# 2) averaged number of parent-offspring pairs found across all iterations *POPs*, 
# 3) averaged variance *var*, 
# 4) averaged standard error *se*, 
# 6) averaged coefficient of variation *cv*, 
# 7) averaged upper 95% confidence interval *upper_confint*, 
# 8) averaged lower 95% confidence interval *lower_confint*, 
# 9) and averaged bias *bias*.
# ==============================================================================
# tGMRsim function =============================================================
tGMRsim <- function(Nc_adult = 3702, 
                    n_adult = 305, 
                    Nc_offspring = 480000, 
                    n_offspring = 682, 
                    ages = c("1.1", "1.2", "1.3", "1.4"), 
                    age_prob = c(0.1, 0.14, 0.64, 0.12), 
                    select_prob = c(0.03, 0.30, 0.38, 0.29), 
                    rs_prob = c(0.07, 0.17, 0.28, 0.38), 
                    drop_prob = c(0.1, 0.1, 0.1, 0.1),
                    drop_percent = 0.00,
                    iterations = 1000,
                    scenario = 1,
                    delta_selex = 0,
                    delta_rs = 0) {
  
  # Initialize output_dat ======================================================
  output_dat <- data.frame(matrix(ncol = 3, nrow = iterations))
  colnames(output_dat) <- c('Iteration', 'Estimate', 'Scenario')
  # ============================================================================
  
  
for (i in 1:iterations) {
  

  # Initialize adult dataframe =================================================
  adult_dat <- data.frame(matrix(ncol = 3, nrow = Nc_adult))
  colnames(adult_dat) <- c('Adult_ID', 'Adult_Age', 'Adult_Sampled')
  
  # Assign Adult_ID's ==========================================================
  adult_dat$Adult_ID <- 1:Nc_adult
  
  # Assign Adult_Age based on estimated age-structure of adults
  adult_dat$Adult_Age <- sample(ages, size = Nc_adult, replace = TRUE, prob = age_prob)
  # ============================================================================
  # Provide an opportunity for adult dropout to occur (optional and can be age-specific or not)
  
  adult_dat <- adult_dat %>% mutate(Adult_Drop_Prob = ifelse(Adult_Age=="1.1",drop_prob[1],
                                                             ifelse(Adult_Age=="1.2",drop_prob[2],
                                                                    ifelse(Adult_Age=="1.3",drop_prob[3],
                                                                           ifelse(Adult_Age=="1.4",drop_prob[4],NA)))))
  
  # Identify adults that are going to be dropped (vector of IDs)
  dropped_adults <- sample(adult_dat$Adult_ID, size = (Nc_adult * drop_percent), replace = FALSE, prob = adult_dat$Adult_Drop_Prob)
  
  adult_dat <- adult_dat %>% dplyr::filter(!Adult_ID %in% dropped_adults)
  # ============================================================================
  # Sample adults based on selectivity-by-age parameters =======================
  # ============================================================================
  # Using selectivity-by-age, assign adults either a 1 (sampled) or 0 (unsampled)
  adult_dat$Adult_Sampled <- 0
  adult_dat<-adult_dat %>% mutate(Adult_Age_Prob = ifelse(Adult_Age=="1.1",select_prob[1],
                                                          ifelse(Adult_Age=="1.2",select_prob[2],
                                                                 ifelse(Adult_Age=="1.3",select_prob[3],
                                                                        ifelse(Adult_Age=="1.4",select_prob[4],NA)))))

  sampled_adults<-sample(adult_dat$Adult_ID, size = n_adult, replace = FALSE, prob = adult_dat$Adult_Age_Prob)
  adult_dat$Adult_Sampled[which(adult_dat$Adult_ID%in%sampled_adults)] <- 1
  # End with 581 sampled adults and 2,599 unsampled adults.
  # ============================================================================
  # Initialize offspring dataframe =============================================
  offspring_dat <- data.frame(matrix(ncol = 5, nrow = Nc_offspring))
  colnames(offspring_dat) <- c('Offspring_ID', 'Offspring_Sampled', 'Adult_ID', 'Adult_Age', 'Adult_Sampled')
  
  # Assign Offspring_ID's 1:480,000
  offspring_dat$Offspring_ID <- 1:Nc_offspring
  # ============================================================================
  # Sample offspring randomly
  offspring_dat$Offspring_Sampled <- 0
  offspring_dat$Offspring_Sampled[sample(offspring_dat$Offspring_ID, size = n_offspring, replace = FALSE)] <- 1
  
  # Now assign parents to offspring using RS-by-age parameters
  adult_dat <- adult_dat %>% mutate(Adult_RS_Prob = ifelse(Adult_Age=="1.1",rs_prob[1],
                                                          ifelse(Adult_Age=="1.2",rs_prob[2],
                                                                 ifelse(Adult_Age=="1.3",rs_prob[3],
                                                                        ifelse(Adult_Age=="1.4",rs_prob[4],NA)))))
  
  # ============================================================================
  offspring_dat$Adult_ID <- sample(adult_dat$Adult_ID, size = Nc_offspring, replace = TRUE, prob = adult_dat$Adult_RS_Prob)
  
  
  # Now match adult age to adult ID in offspring_dat
  offspring_dat$Adult_Age <- adult_dat$Adult_Age[match(offspring_dat$Adult_ID, adult_dat$Adult_ID)]
  
  # And finally, match the adult sampled info to offspring_dat
  offspring_dat$Adult_Sampled <- adult_dat$Adult_Sampled[match(offspring_dat$Adult_ID, adult_dat$Adult_ID)]
  # ============================================================================
  # Create a new dataframe that contains only *sampled* offspring and adults ===
  POP_dat <- offspring_dat %>% filter(Offspring_Sampled == 1 & Adult_Sampled == 1)
  # ============================================================================
  # Calculate tGMR estimate ====================================================
  tGMR_estimate <- (n_adult * (n_offspring + 1)) / (nrow(POP_dat) + 1)
  # ============================================================================
  # Store estimates in output_dat ==============================================
  output_dat$Iteration[i] <- i
  output_dat$Estimate[i] <- tGMR_estimate 
  output_dat$POPs[i] <- nrow(POP_dat)
  output_dat$Scenario[i] <- scenario
  output_dat$n_adult[i] <- n_adult
  output_dat$n_offspring[i] <- n_offspring
  output_dat$N_adult[i] <- Nc_adult
    
}
  output_dat <- output_dat %>%  summarize(mean = mean(Estimate), 
                          POPs = mean(POPs), 
                          n_adult = mean(n_adult), 
                          n_offspring = mean(n_offspring),
                          var = (n_adult^2 * (n_offspring+1) * (n_offspring-POPs)) / ((POPs+1)^2 * (POPs+2)),
                          cv = sqrt(var)/mean, 
                          se = sqrt(var),
                          upper_confint = mean + (1.96 * se), 
                          lower_confint = mean - (1.96 * se),
                          bias = (mean - Nc_adult) / Nc_adult)
  
  output_dat$Nc_adult <- Nc_adult  
  output_dat$n_adult <- n_adult
  output_dat$Nc_offspring <- Nc_offspring 
  output_dat$n_offspring <- n_offspring  
  output_dat$drop_percent<- drop_percent 
  output_dat$iterations <- iterations 
  output_dat$scenario <- scenario 
  output_dat$delta_selex <- delta_selex
  output_dat$delta_rs <- delta_rs
  
  return(output_dat)
  }
  
# ==============================================================================