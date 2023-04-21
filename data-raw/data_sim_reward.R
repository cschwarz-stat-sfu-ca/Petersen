## code to generate LPTL data for reward tags.
# This generates two types of animals, with either a regular tag which is subject to reporting bias, and/or
# a reward tag which is assumed to be 100% reporting, like a permanent tag.

library(Petersen)

dt_type <- Petersen:::valid_dt_type()[3]

data_sim_reward <-LP_TL_simulate(
       dt_type=dt_type,  #  permanent tag
       N=10000,
       cov1=function(N)         {rep(1,N)},
       cov2=function(cov1)      {rep(1,  length(cov1))},
       p1  =function(cov1, cov2){rep(.1, length(cov1))},
       pST =function(cov1, cov2){rep(.75,length(cov1))},
       rho1=function(cov1, cov2){rep(.70,length(cov1))},
       rho2=function(cov1, cov2){rep(1,  length(cov1))},  # permanent second tag
       p2  =function(cov1, cov2){rep(.1, length(cov1))},
       seed=45985, trace=FALSE)

# we don't have fish with both tags
data_sim_reward$cap_hist <- gsub("1P", "0P", data_sim_reward$cap_hist)

# condense the data
data_sim_reward <- plyr::ddply(data_sim_reward, "cap_hist", plyr::summarize, freq=sum(freq))


data_sim_reward
res1 <- LP_TL_fit(data_sim_reward, dt_type=dt_type, p_model=~1, rho_model=~1)
res1$summary

est1 <- LP_TL_est(res1)
est1$summary

est1$detail$data
est1$detail$data.expand

LP_AICc  (res1, res1)
LP_modavg(res1)
usethis::use_data(data_sim_reward, overwrite = TRUE)
