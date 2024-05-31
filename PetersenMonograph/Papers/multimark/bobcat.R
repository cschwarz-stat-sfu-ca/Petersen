#################################################################################################
#################################################################################################
#################################################################################################

library(multimark)
library(Petersen)
library(plyr)

logit <- function (p){ log(p/(1-p))}
expit <- function(theta){ 1/(1+exp(-theta))}

set.seed(234324)

max.iter <- 100000
n.thin   <- 10



# Can we use bobcat data
data("bobcat")
bobcat

# classify into 2 events; first 4 sampling events; last 4 sampling events; if seen in any of the multiple occasions in an event

bobcat <- cbind(bobcat, V1=apply(bobcat[,c("occ1","occ2","occ3","occ4")],1,max))
bobcat <- cbind(bobcat, V2=apply(bobcat[,c("occ5","occ6","occ7","occ8")],1,max))
head(bobcat)

bobcat.df <- as.data.frame(bobcat)
bobcat.ch.red <- plyr::ddply(bobcat.df,c("V1","V2"), plyr::summarize, 
                             cap_hist=paste0(V1[1],V2[1]),
                             freq=length(V1))
bobcat.ch.red

bobcat.left.fit <- Petersen::LP_fit(bobcat.ch.red[grepl('1', bobcat.ch.red$cap_hist),], p_model=~..time)
bobcat.left.est <- Petersen::LP_est(bobcat.left.fit)
bobcat.left.est$summary

temp <- bobcat.ch.red[grepl('2', bobcat.ch.red$cap_hist),]
temp$cap_hist <- gsub('2','1',temp$cap_hist)
bobcat.right.fit <- Petersen::LP_fit(temp, p_model=~..time)
bobcat.right.est <- Petersen::LP_est(bobcat.right.fit)
bobcat.right.est$summary

# naive simple average
bobcat.naive.N.avg <- (bobcat.left.est$summary$N_hat + bobcat.right.est$summary$N_hat)/2
bobcat.naive.N.avg
bobcat.naive.N.avg.se <- sqrt((bobcat.left.est$summary$N_hat_SE^2 + bobcat.right.est$summary$N_hat_SE^2)/4)
bobcat.naive.N.avg.se

#bobcat.fit <- multimarkClosed(bobcat, mod.p=~time)

bobcat.fit <- multimarkClosed(
  as.matrix(bobcat[,c("V1","V2")]),
  data.type = "never",
  covs = data.frame(),
  mms = NULL,
  mod.p = ~time,
  mod.delta = ~type,
  parms = c("pbeta", "delta", "N"),
  nchains = 3,
  iter = max.iter,
  adapt = 10000,
  bin = 50,
  thin = n.thin,
  burnin = 2000,
  taccept = 0.44,
  tuneadjust = 0.95,
  proppbeta = 0.1,
  propzp = 1,
  propsigmap = 1,
  npoints = 500,
  maxnumbasis = 1,
  a0delta = 1,
  a0alpha = 1,
  b0alpha = 1,
  a = 25,
  mu0 = 0,
  sigma2_mu0 = 1.75,
  a0psi = 1,
  b0psi = 1,
  initial.values = NULL,
  known = integer(),
  printlog = FALSE
)

summary(bobcat.fit$mcmc)
plot(bobcat.fit$mcmc)

bobcat.fit$DM

# get the estimates of p and c
pc <- getprobsClosed(bobcat.fit, link = "logit")
summary(pc)


