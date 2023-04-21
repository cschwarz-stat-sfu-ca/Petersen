## code to prepare `seals_tagloss` dataset goes here
## This is the data in Eberhardt et al (1979), Wildlife Monographs, 63, p.8


#Double tagging experiment initiated in 1958 on the Pribilof Islands (Abegglen et al. 1958).
#  In that year, 34,923 fur seal pups were single tagged and 5,000 double tagged on St. Paul Island.
#  Returns were as follows in the 1961 kill of 48,458 3-year-old males:
#  Single tagged (originally single tagged) 2,098,
#  Double tagged (retaining both tags) 285 (md)
#  Double tagged (retaining only 1 tag) 140 (ms)



data.csv <- textConnection("
 cap_hist, freq, comment
 1000, 32825, singe tag;  not seen again
 1010,  2098, A single tag; recovered
 1111,   285, BB double tagged; recovered with double tags
 111X,   140, Double tagged; recovered with a single tag but cannot tell which tag was lost
 1100,  4575, double tagged; not seen again
 0010, 45935, U - apparently captured for first time at t2 but this includes fish that lost all tags"
 )

data_seals_tagloss <- read.csv(data.csv, header=TRUE, as.is=TRUE, strip.white=TRUE, )
data_seals_tagloss$cap_hist <- formatC(data_seals_tagloss$cap_hist, width=4, digits=0, format="f", flag="0")

#usethis::use_data(data_seals_tagloss, overwrite = TRUE)


seals.res1 <- LPTL_fit(data_seals_tagloss, dt_type="notD", p_model=~1, rho_model=~1, trace=FALSE)
seals.res1$summary
seals.res1$fit

seals.est1 <- LPTL_est(seals.res1, N_hat=~1, conf_level=0.95, trace=FALSE)
seals.est1$summary

seals.est1$detail$data.expand[1:2,c("..tag","rho","rho.se")]

# compute the number of tagged fish recaptured at time 2 with all missing tags
n1 <- sum(data_seals_tagloss$freq[substr(data_seals_tagloss$cap_hist,1,2) %in% c("11","10","01") ])
n1
pST <-sum(data_seals_tagloss$freq[substr(data_seals_tagloss$cap_hist,1,2) %in% c("10","01") ]) / n1 # % single tagged
rho <- seals.est1$detail$data.expand$rho[1]
rho
n2 <- sum(data_seals_tagloss$freq[substr(data_seals_tagloss$cap_hist,3,4) %in% c("11","10","01","1X") ])
n2
p2  <- n2 / seals.est1$summary$N_hat

m2.obs <-sum(data_seals_tagloss$freq[
  !substr(data_seals_tagloss$cap_hist,1,2) %in% c("00") &
  substr(data_seals_tagloss$cap_hist,3,4) %in% c("11","10","01","1X") ])
m2.obs

single.tag.all.lost.recap <- n1 * pST * (1-rho) * p2
double.tag.all.lost.recap <- n1 * (1-pST) * (1-rho)^2 * p2
m2.est <- m2.obs + single.tag.all.lost.recap + double.tag.all.lost.recap
m2.est

seals.petersen <- n1 * n2 / m2.est
seals.petersen
