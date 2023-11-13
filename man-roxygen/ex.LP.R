# fit a simple Petersen model and get the estimated abundance
data(data_rodli)
fit <- Petersen::LP_fit(data=data_rodli, p_model=~..time)
fit$summary
# Now to get the estimated abundance
est <- Petersen::LP_est(fit, N_hat=~1)
est$summary

# repeat the fit with the Chapman correction
# we add an additional animal with history 11
rodli.chapman <- plyr::rbind.fill(data_rodli,
                                  data.frame(cap_hist="11",
                                             freq=1,
                                             comment="Added for Chapman"))

rodli.chapman
fit.chapman <- Petersen::LP_fit(data=rodli.chapman, p_model=~..time)
fit.chapman$summary
# Now to get the estimated abundance
est.chapman <- Petersen::LP_est(fit.chapman, N_hat=~1)
est.chapman$summary


# Example of simple stratification (by sex)
data(data_NorthernPike)
nop.red <- plyr::ddply(data_NorthernPike, c("cap_hist","Sex"), plyr::summarize,
                       freq=sum(freq))
nop.red # reduced capture history to speed execution time of example

# Fit the various models
nop.fit.sex.time   <- Petersen::LP_fit(nop.red, p_model=~-1+Sex:..time)
nop.fit.sex.time$summary

# estimate of overall abundance
nop.est.ALL   <- Petersen::LP_est(nop.fit.sex.time, N=~1)
nop.est.ALL$summary

# estimate of abundance for each sex
nop.est.by.sex <- Petersen::LP_est(nop.fit.sex.time, N=~-1+Sex)
nop.est.by.sex$summary


# Refer to vignettes for example using continuous variable (e.g. length) to model catchability
