# Purpose: Simulate how asymptotic patterns of age-specific patterns of reproductive success and sampling selectivity influence tGMR abundance estimation.

# Contact: sam.rosenbaum@umontana.edu
# ==============================================================================
# Load libraries
library(ggplot2)
library(dplyr)
# ==============================================================================
# Visualize asymptotic patterns that we will draw our age-specific reproductive success and sampling selectivity values from:
age <- rep(seq(3,6,length.out=1000),3)
a50 <- 4.5
delta <- rep(c(-3,0,3),each=1000)
age_pred<-data.frame(age=age,a50=4.5,D=delta)

age_pred$selex<-(1+exp(-1*age_pred$D*(age_pred$age-age_pred$a50)))^-1

eq.labs<-data.frame(D = c(-3,0,3),eq = as.character(as.expression(paste("delta","=",c(-3,0,3)))))

delta_fig <- age_pred %>% ggplot(aes(x = age, y = selex)) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        panel.grid = element_blank()
  ) +
  geom_line() +
  facet_grid(.~D) +
  geom_text(data = eq.labs, x = 2.4,y=0.25,inherit.aes=F,
            aes(label=paste("\u03B4","=",D))) +
  ylab("Prob.") +
  xlab("Age")
# ==============================================================================
# Load asymptotic function designed by Curry Cunningham:
selex_asymp <- function(age, a50, delta) {
  selex <- (1+exp(-1*delta*(age-a50)))^-1
  return(selex)
}
# ==============================================================================
# Simulate asymptotic selectivity:
# ==============================================================================
# Set up object to hold selectivity predictions across a range of deltas ===
selex_pred_dat <- data.frame(matrix(ncol = 5, nrow = 7))
colnames(selex_pred_dat) <- c("Delta", "Age_3", "Age_4", "Age_5", "Age_6")
selex_pred_dat$Delta<- -3:3
# ==============================================================================
# Create a loop to populate "selex_pred_dat" ===================================
# Loop through delta values -3:3 and extract the 4 selectivity at age values.

for (i in 1:7) {
  selex_pred_dat[i,2:5] <- selex_asymp(age=c(3:6), a50=4.5, delta=selex_pred_dat$Delta[i])
}
# ==============================================================================
nc_adults <- 3769
samp_adults_combo <- 583
offspring <- 480000
samp_offspring <- 682
adult_ages <- c("1.1", "1.2", "1.3", "1.4")
age_structure <- c(0.1, 0.14, 0.64, 0.12)
MS_selex <- c(0.41, 0.19, 0.21, 0.19)
TRIB_selex <- c(0.03, 0.30, 0.38, 0.29)
RS_prob <- c(0.07, 0.17, 0.28, 0.38)
# Create a loop to populate select_prob ===================================
# NOTE: This assumes increasing RS prob *rs_prob* = c(0.07, 0.16, 0.28, 0.39)
for(i in 1:7) {
  out_selex <- tGMRsim(select_prob = selex_pred_dat[i,2:5],
                       rs_prob = RS_prob,
                       Nc_adult = nc_adults, 
                       n_adult = samp_adults_combo, 
                       Nc_offspring = offspring, 
                       n_offspring = samp_offspring, 
                       ages = adult_ages, 
                       age_prob = age_structure,
                       drop_prob = c(0.2, 0.2, 0.2, 0.2),
                       drop_percent = 0,
                       iterations = 1000,
                       delta_selex = selex_pred_dat$Delta[i])
  if(i == 1) {
    output_selex <- out_selex
  }
  if(i > 1) {
    output_selex <- rbind(output_selex, out_selex)
  }
}
# ==============================================================================
delta_selex_fig <- output_selex %>% 
  ggplot(aes(x = delta_selex, y = bias)) +
  geom_ribbon(aes(ymin = ((lower_confint - Nc_adult) / Nc_adult) - bias, ymax = ((upper_confint - Nc_adult) / Nc_adult) + bias), alpha = 0.5) +
  geom_line() + 
  xlab(expression(delta[Selectivity])) +
  ylab("") +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  theme_classic() +
  scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
  ylim(-0.5, 1.2)
# ==============================================================================
# Plot simulated estimate and 95% for each delta and overlay a horizontal line with the true population size
selex_estimate_fig <- output_selex %>%
  ggplot(aes(x = delta_selex, y = mean)) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_confint, ymax = upper_confint, width = 0.1)) +
  ylab("") +
  xlab(expression(delta[Selectivity])) +
  ylim(0, 7000) + 
  geom_hline(yintercept = 3702, color = "red", lty = 2) +
  theme_classic() +
  scale_x_continuous(breaks = seq(-3, 3, by = 1))
# ==============================================================================
# ==============================================================================
# Simulate asymptotic reproductive success:
# Set up object to hold reproductive success predictions across a range of deltas 
rs_pred_dat <- data.frame(matrix(ncol = 5, nrow = 7))
colnames(rs_pred_dat) <- c("Delta", "Age_3", "Age_4", "Age_5", "Age_6")
rs_pred_dat$Delta<- -3:3
# ==============================================================================
# Create a loop to populate "rs_pred_dat" ===================================
# Loop through delta values -3:3 and extract the 4 selectivity at age values.

for (i in 1:7) {
  rs_pred_dat[i,2:5] <- selex_asymp(age=c(3:6), a50=4.5, delta=rs_pred_dat$Delta[i])
}
# ==============================================================================
nc_adults <- 3769
samp_adults_combo <- 583
offspring <- 480000
samp_offspring <- 682
adult_ages <- c("1.1", "1.2", "1.3", "1.4")
age_structure <- c(0.1, 0.14, 0.64, 0.12)
selex <- c(0.1, 0.4, 0.7, 0.10)
# Create a loop to populate rs_prob ============================================
# NOTE: This assumes increasing *select_prob* = c(0.1, 0.4, 0.7, 0.10)
for(i in 1:7) {
  out_rs <- tGMRsim(select_prob = selex,
                    rs_prob = rs_pred_dat[i,2:5],
                    Nc_adult = nc_adults, 
                    n_adult = samp_adults_combo, 
                    Nc_offspring = offspring, 
                    n_offspring = samp_offspring, 
                    ages = adult_ages, 
                    age_prob = age_structure,
                    drop_prob = c(0.2, 0.2, 0.2, 0.2),
                    drop_percent = 0,
                    iterations = 1000,
                    delta_rs = rs_pred_dat$Delta[i])
  if(i == 1) {
    output_rs <- out_rs
  }
  if(i > 1) {
    output_rs <- rbind(output_rs, out_rs)
  }
}
# ==============================================================================
delta_rs_fig <- output_rs %>% 
  ggplot(aes(x = delta_rs, y = bias)) +
  geom_ribbon(aes(ymin = ((lower_confint - Nc_adult) / Nc_adult) - bias, ymax = ((upper_confint - Nc_adult) / Nc_adult) + bias), alpha = 0.5) +
  geom_line() + 
  xlab(expression(delta[ReproductiveSuccess])) +
  ylab("Bias") +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  theme_classic() +
  scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
  ylim(-0.5, 1.2)
# ==============================================================================
# Plot simulated estimate and 95% for each delta and overlay a horizontal line with the true population size
rs_estimate_fig <- output_rs %>%
  ggplot(aes(x = delta_rs, y = mean)) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_confint, ymax = upper_confint, width = 0.1)) +
  ylab("Escapement") +
  xlab(expression(delta[ReproductiveSuccess])) +
  ylim(0, 7000) + 
  geom_hline(yintercept = 3702, color = "red", lty = 2) +
  theme_classic() +
  scale_x_continuous(breaks = seq(-3, 3, by = 1))
# ==============================================================================
# ==============================================================================
# Visualize everything together:
asymp_fig <- ggarrange(delta_rs_fig, delta_selex_fig, nrow = 1, common.legend = TRUE)

delta_estimates <- ggarrange(rs_estimate_fig, selex_estimate_fig, nrow = 1, common.legend = TRUE)

ggarrange(delta_fig, asymp_fig, delta_estimates, nrow = 3, common.legend = TRUE)
# ==============================================================================
