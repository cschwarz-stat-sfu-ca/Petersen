# ==============================================================================
# Purpose: Calculate parameters for mainstem and spawning ground (tributary) tGMR simulations

# Contact: sam.rosenbaum@umontana.edu
# ==============================================================================
# Load Libraries
library(stringr)
library(tidyverse)
library(ggplot2)
library(visreg)
library(lubridate)
library(coefplot)
library(PNWColors)
# ==============================================================================
# Read in pedigree outputs from https://github.com/krshedd/Chilkat-Chinook-tGMR:

# Load in ParentPair data from COLONY pedigree
parentpair_data <- read.csv(file = "/Users/samrosenbaum/Desktop/tGMR_reviews/Code/ParentPair.csv", header = TRUE, sep = ",")

# Load in field metadata
field_dat <- read.csv(file = "/Users/samrosenbaum/Desktop/tGMR_reviews/Code/adult_parr_field_age.csv")
# ==============================================================================
# Create a binary reproductive success column (1 = successful, 0 = unsuccessful)

# Load in all adult IDs 
Chilkat_adults <- read.csv(file ="/Users/samrosenbaum/Desktop/tGMR_reviews/Code/Chilkat_adults_data.csv")
adults <- as.data.frame(Chilkat_adults)
adult_ids <- as.data.frame(adults$SillySource) # Just want the ids
names(adult_ids)[1] <- 'id' # rename the column to 'id'

# Find unique parents assigned to offspring
rs_dad <- as.data.frame(unique(parentpair_data$InferredDad))
rs_mum <- as.data.frame(unique(parentpair_data$InferredMum))

# Remove the "*" and "#" symbols so that only actual adults remain
rs_dad <- rs_dad[-1,]
rs_dad <- as.data.frame(rs_dad)
names(rs_dad)[1] <- 'id' # rename the column so we can join later

rs_mum <- rs_mum[-1,]
rs_mum <- as.data.frame(rs_mum)
names(rs_mum)[1] <- 'id' # rename the column so we can join later

# Sex doesn't matter, so put all rs adults into one column
rs_adults <- rbind(rs_dad, rs_mum)
sum(duplicated(rs_adults$id)) # No duplicates
nrow(rs_adults) # 148

# Make new object that holds only unique ids
rs_adults_unique <- as.data.frame(unique(rs_adults$id))
names(rs_adults_unique)[1] <- 'id' # rename the column so we can join later
nrow(rs_adults_unique) # 148, good

# Create a column signifying reproductive success as '1'
rs <- rep(1, 148)
rs_adults_unique$rs <- rs # Add the column

# Join the rs column to the dataframe holding all ids
all_adults <- adult_ids %>% left_join(rs_adults_unique, by = "id")

# Change na's to 0's to signify lack of reproductive success
all_adults[is.na(all_adults)] <- 0
str(all_adults)
# ==============================================================================
# Create a column that we can join with
field_dat$id <- paste(field_dat$silly, field_dat$fish_id_LOKI, sep = "_")

# Join the field data to the rs data
rs_field <- all_adults %>% left_join(field_dat, by = "id") %>% 
  filter(!is.na(fish_id_LOKI))

# Remove columns that won't be used 
param_dat <- rs_field[,-c(3, 4, 6, 8, 11:14, 16, 17, 19:24)]

# Convert the variables to their proper type
param_dat$MEF <- as.numeric(param_dat$MEF) # length as numeric

param_dat$fact_rs <- as.factor(param_dat$rs) # fact_rs as a factor

param_dat$event <- as.factor(param_dat$event) # event as a factor

param_dat <- param_dat %>% 
  mutate(date = mdy(date),
         doy = yday(date)) # doy as date format

param_dat$AEC <- as.factor(param_dat$AEC) # AEC as a factor

# Now we have a data frame with only the variables we want
str(param_dat)
unique(param_dat$AEC)
# ==============================================================================
# Dataframe holding only age 1.1 adults
chil_11 <- param_dat[param_dat$AEC == 11,]
nrow(chil_11) # 42
# 42/434 = 0.096 = 10% of the sampled adults
mean(chil_11$MEF) # 372.5 mm

# Dataframe holding only age 1.2 adults
chil_12 <- param_dat[param_dat$AEC == 12,]
nrow(chil_12) # 60
# 60/434 =  0.138 = 14% of sampled adults
mean(chil_12$MEF) # 604.0833 mm

# Dataframe holding only age 1.3 adults
chil_13 <- na.omit(param_dat[param_dat$AEC == 13,])
nrow(chil_13) # 280
# 280/434 = 0.645 = 65% of the sampled adults
mean(chil_13$MEF, ) # 785.5357 mm

# Dataframe holding only age 1.4 adults
chil_14 <- param_dat[param_dat$AEC == 14,]
nrow(chil_14) # 52
# 52/434 = 0.119 = 12% of the sampled adults
mean(chil_14$MEF, ) # 847.1154 mm

# Dataframe holding all individuals with an assigned age
chil.dat <- rbind(chil_11, chil_12, chil_13, chil_14)
nrow(chil.dat) # 434

str(chil.dat)
unique(chil.dat$AEC)
# ==============================================================================
# Selectivity Figures ==========================================================

# Age counts by location
ggplot(chil.dat, aes(x=AEC, fill = location)) +
  geom_bar(alpha = 0.8) +
  facet_wrap(~location) +
  xlab("Age") + 
  ylab("Count") +
  theme_classic()

# Age counts by event
ggplot(chil.dat, aes(x=AEC, fill = event)) +
  geom_bar(alpha = 0.8) +
  facet_wrap(~event) +
  xlab("Age") + 
  ylab("Count") +
  theme_classic()
# ==============================================================================
# Reproductive Success Figures ================================================+

# Age counts by RS 
chil.dat %>% filter(fact_rs==1) %>%
  ggplot(aes(x=AEC, fill = AEC)) +
  geom_bar(alpha = 0.8) +
  xlab("Age") + 
  ylab("Count (Empirical)") +
  theme_classic()

# Age counts by RS
ggplot(chil.dat, aes(x=AEC, fill = fact_rs)) +
  geom_bar(alpha = 0.8) +
  facet_wrap(~fact_rs) +
  xlab("Age") + 
  ylab("Count") +
  theme_classic()

# ==============================================================================
# ==============================================================================
# "True" (Mainstem) Age-Structure for Chilkat Chinook:

# 42/434 = 0.096 = 10% of the sampled adults were 1.1
# 60/434 =  0.138 = 14% of sampled adults were 1.2
# 280/434 = 0.645 = 64% of the sampled adults were 1.3
# 52/434 = 0.119 = 12% of the sampled adults were 1.4

# ==============================================================================
# ==============================================================================
# ==============================================================================
# Create a mainstem dataframe ==================================================

mainstem <- chil.dat %>% 
  filter(location %in% c("fishwheel", "gillnet"))

# Subset by gear used in the mainstem

fishwheel <- chil.dat %>% 
  filter(location %in% c("fishwheel"))
nrow(fishwheel)

gillnet <- chil.dat %>% 
  filter(location %in% c("gillnet"))
nrow(gillnet)

# Calculate number of fish in mainstem sample
nrow(mainstem) # 222 fish

# Calculate mainstem age counts
nrow(mainstem[mainstem$AEC == 11,]) # 40
nrow(mainstem[mainstem$AEC == 12,]) # 27
nrow(mainstem[mainstem$AEC == 13,]) # 132
nrow(mainstem[mainstem$AEC == 14,]) # 23

# Calculate mainstem age-structure
# n_adults of age x / n_adults sampled
ms_age_1.1 <- 40 / 222 # 0.18
ms_age_1.2 <- 27 / 222 # 0.12
ms_age_1.3 <- 132 / 222 # 0.60
ms_age_1.4 <- 23 / 222 # 0.10

# MAINSTEM AGE-STRUCTURE = c(0.18, 0.12, 0.60, 0.10)
mainstem_age_prob <- c(0.18, 0.12, 0.60, 0.10)

# Calculate Mainstem SELECTIVITY (age-specific relative vulnerability of capture)

# Divide mainstem age-structure by 'true' age structure to get 'relative selectivity'
# In this case, we're assuming the mainstem age-structure is the truth
relative_selex_ms <- c(0.18, 0.12, 0.60, 0.10)

# Divide relative selectivity by its sum to standardize
standardized_selex_ms <- relative_selex_ms / sum(relative_selex_ms)
standardized_selex_ms # 0.18 0.12 0.60 0.10
ms_selex <- c(0.18, 0.12, 0.60, 0.10)

# ==============================================================================
# GLM with age as the predictor for mainstem
mainstem_mod <- glm(fact_rs ~ AEC, data = mainstem, family = binomial (link="logit"))
summary(mainstem_mod) # AIC = 213.55

visreg(mainstem_mod, scale="response")
coefplot(mainstem_mod, intercept = FALSE)

mainstem_pred.data <- data.frame(AEC=c("11","12", "13", "14"), fact_rs=NA)
mainstem_pred.data$fact_rs <- predict(mainstem_mod, newdata = mainstem_pred.data, type="response")
mainstem_pred.data$se.fit <- predict(mainstem_mod, newdata = mainstem_pred.data, type="response",se.fit = T)$se.fit

mainstem %>% group_by(AEC) %>% summarize(mean_AEC=mean(as.numeric(fact_rs),na.rm=T))

mainstem_pred.data$AEC <- as.factor(mainstem_pred.data$AEC)

str(mainstem_pred.data)

# Mainstem reproductive success
mainstem_rs_prob <- c(0.08, 0.11, 0.22, 0.26)
# ==============================================================================
# ==============================================================================
# ==============================================================================

# Spawning Grounds =============================================================
# ==============================================================================
# ==============================================================================
# Create a spawning grounds dataframe ==================================================

spawning_grounds <- chil.dat %>% 
  filter(location %in% c("Tahni", "Kelsall", "Klehini"))

# Calculate number of fish in spawning grounds sample
nrow(spawning_grounds) # 212 fish

# Calculate mainstem fish age-structure
nrow(spawning_grounds[spawning_grounds$AEC == 11,]) # 2
nrow(spawning_grounds[spawning_grounds$AEC == 12,]) # 33
nrow(spawning_grounds[spawning_grounds$AEC == 13,]) # 148
nrow(spawning_grounds[spawning_grounds$AEC == 14,]) # 29

# ==============================================================================
# Calculate tributary age counts
# n_adults of age x / n_adults sampled
tr_age_1.1 <- 2 / 212 # 0.01
tr_age_1.2 <- 33 / 212 # 0.16
tr_age_1.3 <- 148 / 212 # 0.69
tr_age_1.4 <- 29 / 212 # 0.14

# Tributary age-structure = c(0.01, 0.12, 0.69, 0.10)
tributary_age_prob <- c(0.01, 0.12, 0.69, 0.10)

# Calculate tributary SELECTIVITY (age-specific relative vulnerability of capture) 
# Plot out age comps (do this for Curry!)

# Divide tributary age-structure by 'true' (mainstem here) age structure to get 'relative selectivity'
relative_selex_trib <- tributary_age_prob / mainstem_age_prob

# Divide relative selectivity by its sum to standardize
standardized_selex_trib <- relative_selex_trib / sum(relative_selex_trib)
standardized_selex_trib # 0.01733102 0.31195841 0.35875217 0.31195841
trib_selex <- c(0.02, 0.31, 0.36, 0.31)

# ==============================================================================
par(mfcol=c(3,2),mar=c(4,4,4,4))
barplot(mainstem_age_prob, ylim = c(0.0, 0.8))
barplot(ms_selex, ylim = c(0.0, 0.8))
barplot(mainstem_age_prob, ylim = c(0.0, 0.8))
barplot(mainstem_age_prob, ylim = c(0.0, 0.8))
barplot(trib_selex, ylim = c(0.0, 0.8))
barplot(tributary_age_prob, ylim = c(0.0, 0.8))


mainstem_age_prob
ms_selex

tributary_age_prob
trib_selex

par(mfcol=c(3,2),mar=c(4,4,4,4))
barplot(mainstem_age_prob, ylim = c(0.0, 0.8))
barplot(mainstem_age_prob, ylim = c(0.0, 0.8))
barplot(ms_selex, ylim = c(0.0, 0.8))
barplot(mainstem_age_prob, ylim = c(0.0, 0.8))
barplot(tributary_age_prob, ylim = c(0.0, 0.8))
barplot(trib_selex, ylim = c(0.0, 0.8))

# ==============================================================================
# GLM with age as the predictor for spawning grounds
spawning_grounds_mod <- glm(fact_rs ~ AEC, data = spawning_grounds, family = binomial (link="logit"))
summary(spawning_grounds_mod) # AIC = 270.22

visreg(spawning_grounds_mod, scale="response")
coefplot(spawning_grounds_mod, intercept = FALSE)

spawning_grounds_pred.data <- data.frame(AEC=c("11","12", "13", "14"), fact_rs=NA)
spawning_grounds_pred.data$fact_rs <- predict(spawning_grounds_mod, newdata = spawning_grounds_pred.data, type="response")
spawning_grounds_pred.data$se.fit <- predict(spawning_grounds_mod, newdata = spawning_grounds_pred.data, type="response",se.fit = T)$se.fit

spawning_grounds %>% group_by(AEC) %>% summarize(mean_AEC=mean(as.numeric(fact_rs),na.rm=T))

spawning_grounds_pred.data

# Spawning grounds reproductive success
tributary_rs_prob <- c(0.00, 0.21, 0.33, 0.48)

# ==============================================================================
# ==============================================================================
# ==============================================================================

# Calculate reproductive success for ALL adults

# GLM with age as the predictor for all adults
all_adults_mod <- glm(fact_rs ~ AEC, data = chil.dat, family = binomial (link="logit"))
summary(spawning_grounds_mod) # AIC = 270.22

visreg(all_adults_mod, scale="response")
coefplot(all_adults_mod, intercept = FALSE)

all_adults_pred.data <- data.frame(AEC=c("11","12", "13", "14"), fact_rs=NA)
all_adults_pred.data$fact_rs <- predict(all_adults_mod, newdata = all_adults_pred.data, type="response")
all_adults_pred.data$se.fit <- predict(all_adults_mod, newdata = all_adults_pred.data, type="response",se.fit = T)$se.fit

chil.dat %>% group_by(AEC) %>% summarize(mean_AEC=mean(as.numeric(fact_rs),na.rm=T))

all_adults_pred.data

# All adults reproductive success
all_adults_rs_prob <- c(0.07, 0.17, 0.28, 0.38)

nrow(chil.dat)

