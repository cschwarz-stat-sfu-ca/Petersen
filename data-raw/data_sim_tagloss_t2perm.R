## code to generate LPTL data for "t2perm" example

library(Petersen)

dt_type <- Petersen:::valid_dt_type()[3]

data_sim_tagloss_t2perm <-LP_TL_simulate(
       dt_type=dt_type,  #  permanent tag
       N=10000,
       cov1=function(N)         {rep(1,N)},
       cov2=function(cov1)      {rep(1,  length(cov1))},
       p1  =function(cov1, cov2){rep(.1, length(cov1))},
       pST =function(cov1, cov2){rep(.25,length(cov1))},
       rho1=function(cov1, cov2){rep(.70,length(cov1))},
       rho2=function(cov1, cov2){rep(1,  length(cov1))},  # permanent second tag
       p2  =function(cov1, cov2){rep(.1, length(cov1))},
       seed=34535, trace=FALSE)

data_sim_tagloss_t2perm
res1 <- LP_TL_fit(data_sim_tagloss_t2perm, dt_type=dt_type, p_model=~1, rho_model=~1)
res1$summary

est1 <- LP_TL_est(res1)
est1$summary

est1$detail$data
est1$detail$data.expand

LP_AICc  (res1, res1)
LP_modavg(res1)
usethis::use_data(data_sim_tagloss_t2perm, overwrite = TRUE)
