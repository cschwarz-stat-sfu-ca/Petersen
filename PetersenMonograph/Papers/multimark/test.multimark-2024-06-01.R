# Testing out multimark package


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
# Simple model with equal p and equal delta and no c.
# fit a model with equal p as well.

test <- simdataClosed(
  N = 5000,
  noccas = 2,
  pbeta = logit(c(.25, .25)),
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

# Notice that we are expecting N*p1 and N*pn2 captures at t=1,2
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
  mod.p = ~1,
  mod.delta = ~type,
  parms = c("pbeta", "delta", "N"),
  nchains = 3,
  iter = 50000,
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

# Everything seems to be ok.
# Estimates of N, p, are bang on.
# Estimates of delta appear to be slight too large

plot(fit$mcmc)




#################################################################################################
#################################################################################################
#################################################################################################
# Simulated data equal p and equal delta and no c.
# Fit a model with a p~time effect

test <- simdataClosed(
  N = 5000,
  noccas = 2,
  pbeta = logit(c(.25, .25)),
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

# The estimates of N and p seem to be ok.
# Estimates of delta again appear to be biased upwards

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
# Try with a wider spread of the p's
# Fit a model p~time

test <- simdataClosed(
  N = 5000,
  noccas = 2,
  pbeta = logit(c(.5, .7)),
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

# Notice that we are expecting N*p1 an N*p2 captures at t=1,2
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
temp <- test.ch.red[grepl('1', test.ch.red$cap_hist),]

left.fit <- Petersen::LP_fit(temp, p_model=~..time)
left.est <- Petersen::LP_est(left.fit)
left.est$summary
left.est$detail$data.expand[1:2,]
c(.5,.7)*.6  # product of capture probs x delta_1+(1-delta1-delta2) (recaptured left or both)

temp <- test.ch.red[grepl('2', test.ch.red$cap_hist),]
temp$cap_hist <- gsub('2','1',temp$cap_hist)
Petersen::LP_summary_stats(temp)

right.fit <- Petersen::LP_fit(temp, p_model=~..time)
right.est <- Petersen::LP_est(right.fit)
right.est$summary
right.est$detail$data.expand[1:2,]

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
# Try with a wider spread of the p's
# Fit a model p~time using the covariate vector

myCovs <- data.frame(mytime=c(0,1))
myCovs

fit <- multimarkClosed(
  test$Enc.Mat,
  data.type = "never",
  covs = myCovs,
  mms = NULL,
  mod.p = ~mytime,
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

expit(sum(temp$statistics[c(1),   "Mean"]))  # p1
expit(sum(temp$statistics[c(1,2), "Mean"]))  # p2

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
# Try with a wider spread of the p's
# Fit a model p~time using the covariate vector
# Try setting initial values

myCovs <- data.frame(mytime=c(0,1))
myCovs

init.pbeta <- c(logit(.5), logit(.7)-logit(.5))


my.inits <- list( list(pbeta=init.pbeta, delta_1=.4, delta_2=.4, N=7000),   # N must exceed observed count
                  list(pbeta=init.pbeta, delta_1=.4, delta_2=.4, N=7000), 
                  list(pbeta=init.pbeta, delta_1=.4, delta_2=.4, N=7000))
my.inits <- list( list(), list(), list())

test <- simdataClosed(
  N = 5000,
  noccas = 2,
  pbeta = logit(c(.4, .8)),
  tau = 0, # behavioural effects (i.e. c)
  sigma2_zp = 0,
  delta_1 = 0.3,
  delta_2 = 0.6,
  alpha = 0,
  data.type = "never",
  link = "logit"
)


fit <- multimarkClosed(
  test$Enc.Mat,
  data.type = "never",
  covs = myCovs,
  mms = NULL,
  mod.p = ~mytime,
  mod.delta = ~1, #~type,
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
  initial.values = my.inits,
  known = integer(),
  printlog = FALSE
)

summary(fit$mcmc)

#head(fit$mcmc[[1]])
colMeans(fit$mcmc[[1]])
GGally::ggpairs( fit$mcmc[[1]])

# In loglikelhoodClosed we have
#  pstar <- 1-min(1.-tol,max(tol,prod(1-expit(DM$p%*%pbeta))))
# Let us compute pstar with the given data
temp <- fit$mcmc[[1]]
head(temp)
ps <- expit(temp[,1:2] %*% t(fit$DM$p))
pstar <- 1- apply(1-ps,1,prod)
mean(pstar)
# theoretical pstar
pstar.true <- 1- prod(1-c(.5, .7))
cat("Theoretical pstar ", pstar.true, "\n")



temp<- summary(fit$mcmc)

fit$DM

expit(sum(temp$statistics[c(1),   "Mean"]))  # p1
expit(sum(temp$statistics[c(1,2), "Mean"]))  # p2

# get the estimates of p and c
pc <- getprobsClosed(fit, link = "logit")
summary(pc)

.6/.5#  Hmmm... the estimates of p appear to be "wrong", but the estimates of
# N appear to be ok. This suggests that p* (prob of ever seeing an animal)
# has been computed properly. This further suggests that the summary function
# may be off?

# lets have our own look at the MCMC output

head(pc[[1]])
colMeans(pc[[1]])



plot(fit$mcmc)

sessionInfo()
