# Testing out multimark package
# Testing if LRB dataset resolves confounding

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
# Fit a model p~time

test <- simdataClosed(
  N = 5000,
  noccas = 2,
  pbeta = logit(c(.2, .3)),
  tau = 0, # behavioural effects (i.e. c)
  sigma2_zp = 0,
  delta_1 = 0.4,
  delta_2 = 0.4,
  alpha = 0,
  data.type = "never",
  link = "logit"
)
test$Enc.Mat[1:10,]

# get the standard Petersen on Left/right along
test.ch <- as.data.frame(test$Enc.Mat)

test.ch.red <- plyr::ddply(test.ch,c("V1","V2"), plyr::summarize,
                           cap_hist=paste0(V1[1],V2[1]),
                           freq=length(V1))
test.ch.red

# Notice that we are expecting N*p1 and N*p2 captures at t=1,2
# however we find that total "recording" at each sample time is

sum(test.ch.red$freq[ test.ch.red$V1 > 0])
sum(test.ch.red$freq[ test.ch.red$V2 > 0])

# This indicates that some animals were detected on both sides, but we could not make
# the connection between the two sides



# However, if we look at the true encounter matrix, everything is ok
unique(test$trueEnc.Mat)

test.true.ch <- as.data.frame(test$trueEnc.Mat)

test.true.ch.red <- plyr::ddply(test.true.ch,c("V1","V2"), plyr::summarize,
                                cap_hist=paste0(V1[1],V2[1]),
                                freq=length(V1))
test.true.ch.red

sum(test.true.ch.red$freq[ test.true.ch.red$V1 > 0])
sum(test.true.ch.red$freq[ test.true.ch.red$V2 > 0])



# do fits based on left or right only

left.fit <- Petersen::LP_fit(test.ch.red[grepl('1', test.ch.red$cap_hist),], p_model=~..time)
left.est <- Petersen::LP_est(left.fit)
left.est$summary

temp <- test.ch.red[grepl('2', test.ch.red$cap_hist),]
temp$cap_hist <- gsub('2','1',temp$cap_hist)
right.fit <- Petersen::LP_fit(temp, p_model=~..time)
right.est <- Petersen::LP_est(right.fit)
right.est$summary

# naive simple average
naive.N.avg <- (left.est$summary$N_hat + right.est$summary$N_hat)/2
naive.N.avg
naive.N.avg.se <- sqrt((left.est$summary$N_hat_SE^2 + right.est$summary$N_hat_SE^2)/4)
naive.N.avg.se


fit <- multimarkClosed(
  test$Enc.Mat,
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

summary(fit$mcmc)

temp<- summary(fit$mcmc)


fit$DM

expit(sum(temp$statistics[c(1),   "Mean"]))  # intercept
expit(sum(temp$statistics[c(1,2), "Mean"]))  # p1
expit(sum(temp$statistics[c(1,3), "Mean"]))  # p2

# get the estimates of p and c
pc <- getprobsClosed(fit, link = "logit")
summary(pc)


# If you look at the R code on GitHub for Closed.R, you find
# p <- expittol(matrix(rep(DM$p%*%pbeta,each=n)*firstcap+rep(DM$c%*%pbeta,each=n)*(1-firstcap)+zp[which(H>1)],nrow=n,ncol=noccas))
# This seem to suggests that something weird is going on in the design matrices?


plot(fit$mcmc)

#################################################################################################
#################################################################################################
#################################################################################################
# Simple model with unequal p, no c
# Fit a model p~time

test <- simdataClosed(
  N = 5000,
  noccas = 2,
  pbeta = logit(c(.2, .3)),
  tau = 0, # behavioural effects (i.e. c)
  sigma2_zp = 0,
  delta_1 = 0.4,
  delta_2 = 0.4,
  alpha = 0,
  data.type = "always",
  link = "logit"
)
test$Enc.Mat[1:10,]

# get the standard Petersen on Left/right along
test.ch <- as.data.frame(test$Enc.Mat)

test.ch.red <- plyr::ddply(test.ch,c("V1","V2"), plyr::summarize,
                           cap_hist=paste0(V1[1],V2[1]),
                           freq=length(V1))
test.ch.red

# Notice that we are expecting N*p1 and N*p2 captures at t=1,2
# however we find that total "recording" at each sample time is

sum(test.ch.red$freq[ test.ch.red$V1 %in% c(1,4)])
sum(test.ch.red$freq[ test.ch.red$V2 %in% c(2,4)])

# This indicates that some animals were detected on both sides, but we could not make
# the connection between the two sides



# However, if we look at the true encounter matrix, everything is ok
unique(test$trueEnc.Mat)

test.true.ch <- as.data.frame(test$trueEnc.Mat)

test.true.ch.red <- plyr::ddply(test.true.ch,c("V1","V2"), plyr::summarize,
                                cap_hist=paste0(V1[1],V2[1]),
                                freq=length(V1))
test.true.ch.red

sum(test.true.ch.red$freq[ test.true.ch.red$V1 %in% c(1,4)])
sum(test.true.ch.red$freq[ test.true.ch.red$V2 %in% c(2,4)])



# do fits based on left or right only

temp <- test.ch.red
temp$cap_hist <- gsub('4','1',temp$cap_hist)
temp$cap_hist <- gsub('2','0',temp$cap_hist)
temp <- temp[grepl('1', temp$cap_hist),]
temp

left.fit <- Petersen::LP_fit(temp, p_model=~..time)
left.est <- Petersen::LP_est(left.fit)
left.est$summary

temp <- test.ch.red
temp$cap_hist <- gsub('4','2',temp$cap_hist)
temp$cap_hist <- gsub('1','0',temp$cap_hist)
temp$cap_hist <- gsub('2','1',temp$cap_hist)
temp <- temp[grepl('1', temp$cap_hist),]
temp
tempright.fit <- Petersen::LP_fit(temp, p_model=~..time)
right.est <- Petersen::LP_est(right.fit)
right.est$summary

# naive simple average
naive.N.avg <- (left.est$summary$N_hat + right.est$summary$N_hat)/2
naive.N.avg
naive.N.avg.se <- sqrt((left.est$summary$N_hat_SE^2 + right.est$summary$N_hat_SE^2)/4)
naive.N.avg.se


fit <- multimarkClosed(
  test$Enc.Mat,
  data.type = "always",
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

summary(fit$mcmc)

temp<- summary(fit$mcmc)


fit$DM

expit(sum(temp$statistics[c(1),   "Mean"]))  # intercept
expit(sum(temp$statistics[c(1,2), "Mean"]))  # p1
expit(sum(temp$statistics[c(1,3), "Mean"]))  # p2

# get the estimates of p and c
pc <- getprobsClosed(fit, link = "logit")
summary(pc)


# If you look at the R code on GitHub for Closed.R, you find
# p <- expittol(matrix(rep(DM$p%*%pbeta,each=n)*firstcap+rep(DM$c%*%pbeta,each=n)*(1-firstcap)+zp[which(H>1)],nrow=n,ncol=noccas))
# This seem to suggests that something weird is going on in the design matrices?


plot(fit$mcmc)










