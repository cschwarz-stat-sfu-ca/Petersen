## code to generate LPTL data for "twoD" example

data_sim_tagloss_twoD <-LPTL_simulate(
       dt_type=valid_dt_type()[2],  # two distinguishable tags
       N=10000,
       cov1=function(N)         {rep(1,N)},
       cov2=function(cov1)      {rep(1,  length(cov1))},
       p1  =function(cov1, cov2){rep(.1, length(cov1))},
       pST =function(cov1, cov2){rep(.25,length(cov1))},
       rho1=function(cov1, cov2){rep(.70,length(cov1))},
       rho2=function(cov1, cov2){rep(.80,length(cov1))},
       p2  =function(cov1, cov2){rep(.1, length(cov1))},
       seed=234523, trace=FALSE)

data_sim_tagloss_twoD
library(Petersen)
res1 <- LPTL_fit(data_sim_tagloss_twoD, dt_type=valid_dt_type()[2], p_model=~1, rho_model=~1)
res1$summary

est1 <- LPTL_est(res1)
est1$summary

est1$detail$data
est1$detail$data.expand


res2 <- LPTL_fit(data_sim_tagloss_twoD, dt_type=valid_dt_type()[2], p_model=~1, rho_model=~-1+..tag)
res2$summary

est2 <- LPTL_est(res2)
est2$summary

est2$detail$data
est2$detail$data.expand

LP_AICc( res1, res2)
LP_modavg(res1, res2)
usethis::use_data(data_sim_tagloss_twoD, overwrite = TRUE)
