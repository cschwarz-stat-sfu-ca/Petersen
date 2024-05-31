# Testing out multimark package with more sampling occasions


library(multimark)
library(Petersen)
library(plyr)

logit <- function (p){ log(p/(1-p))}
expit <- function(theta){ 1/(1+exp(-theta))}

set.seed(234324)

max.iter <- 100000
n.thin   <- 10


#################################################################################################
#################################################################################################
#################################################################################################
# Simple model with unequal p, no c
# Try with a wider spread of the p's
# Fit a model p~time using the covariate vector
# Try setting initial values

myCovs <- data.frame(mytime2=c(0,1,0),
                     mytime3=c(0,0,1))
myCovs

init.pbeta <- c(logit(.3), logit(.5)-logit(.3), logit(.7)-logit(.3))



test <- simdataClosed(
  N = 5000,
  noccas = 3,
  pbeta = logit(c(.3, .5, .7)),
  tau = 0, # behavioural effects (i.e. c)
  sigma2_zp = 0,
  delta_1 = 0.4,
  delta_2 = 0.4,
  alpha = 0,
  data.type = "never",
  link = "logit"
)


fit <- multimarkClosed(
  test$Enc.Mat,
  data.type = "never",
  covs = myCovs,
  mms = NULL,
  mod.p = ~mytime2+mytime3,
  mod.delta = ~1,
  parms = c("pbeta", "delta", "N","psi"),
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

summary(fit$mcmc)

temp<- summary(fit$mcmc)

fit$DM

expit(sum(temp$statistics[c(1),   "Mean"]))  # p1
expit(sum(temp$statistics[c(1,2), "Mean"]))  # p2
expit(sum(temp$statistics[c(1,3), "Mean"]))  # p3

# get the estimates of p and c
pc <- getprobsClosed(fit, link = "logit")
summary(pc)


plot(fit$mcmc)