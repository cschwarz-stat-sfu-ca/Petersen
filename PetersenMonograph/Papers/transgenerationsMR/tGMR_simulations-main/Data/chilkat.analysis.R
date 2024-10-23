# Look at the Chilkat genetic data for transgenerational capture-recapture


library(msm)
library(Petersen)
library(plyr)
library(tidyr)

# get the parents
parents <- read.csv("Chilkat_adults_data.csv", header=TRUE, as.is=TRUE, strip.white=TRUE)
head(parents[,c("FK_FISH_ID","SillySource")])
parents <- parents[,"SillySource",drop=FALSE]
parents <- plyr::rename(parents, c("SillySource"="ParentID"))

juveniles <- read.csv("ParentPair.csv", header=TRUE, as.is=TRUE, strip.white=TRUE)
xtabs(~Probability, data=juveniles, exclude=NULL, na.action=na.pass)
juveniles <- juveniles[ juveniles$Probability > .35,]
juveniles$Probability <- NULL
head(juveniles)

juveniles.long <- tidyr::pivot_longer(juveniles,
                                     cols=c("InferredDad","InferredMum"),
                                     names_to="ParentSex",
                                     values_to="ParentID")
juveniles.long$ParentID <- gsub("*","", juveniles.long$ParentID, fixed=TRUE)
juveniles.long$ParentID <- gsub("#","", juveniles.long$ParentID, fixed=TRUE)
head(juveniles.long)


cap_hist <- parents
cap_hist$release <- 1
cap_hist$recaps  <- apply( outer(parents$ParentID, juveniles.long$ParentID,FUN="=="),1,sum)
cap_hist$freq <- 1
sum(cap_hist$recaps)

# This not fully "satisfied" because some of the juveniles that have no parents in the marked fish
# are duplicated as well, but this is not possible to identify given the data because
# the inferredMum and inferredDad columns are not given if not part of the marked sample.

no.recaps <- juveniles.long[ juveniles.long$ParentID=="",]
no.recaps$release <- 0
no.recaps$recaps  <- 1
no.recaps$freq    <- 1

cap_hist <- rbind(cap_hist,
                  no.recaps[,c("ParentID","release","recaps","freq")])
cap_hist$cap_hist <- paste0(cap_hist$release, as.numeric(cap_hist$recap>0))

cap_hist.red <- plyr::ddply(cap_hist, c("cap_hist","recaps"), plyr::summarize,
                            freq=sum(freq))
cap_hist.red

# Estimate p using a binomial model

cap_hist.red$n        <- cap_hist.red$freq * cap_hist.red$recaps
cap_hist.red$Recap    <- as.numeric(substr(cap_hist.red$cap_hist,1,1)=="1")*cap_hist.red$n
cap_hist.red$NotRecap <- cap_hist.red$n - cap_hist.red$Recap
select <- cap_hist.red$cap_hist %in% c("01","11")
fit <- glm( cbind(Recap,NotRecap) ~ 1, data=cap_hist.red[select,], family=binomial(link="logit"), start=-2)

sum(cap_hist.red[select,]$Recap) / sum(cap_hist.red[select,c("Recap",'NotRecap')])
summary(fit)

logit.p <- coef(fit)[1]
p <- 1/(1+exp(-logit.p))
p

inv.p <- 1/p
inv.p

inv.p# now call msm to get se of 1/p
inv.p.se <- msm::deltamethod( ~(1+exp(-x1)), mean=coef(fit)[1], cov=vcov(fit)[1,1], ses=TRUE)
inv.p.se

# compute N_hat as M /p


# Now for Binomial, all adult values
M <- sum(cap_hist[ cap_hist$cap_hist %in% c("10","11"),]$freq )
M

N_hat = M *inv.p
N_hat.se = M*inv.p.se
c(N_hat, N_hat.se)

C <- sum(cap_hist[ cap_hist$cap_hist %in% c("01","11"),]$freq * cap_hist[ cap_hist$cap_hist %in% c("01","11"),]$recaps)
C

R <- sum(cap_hist[ cap_hist$cap_hist %in% c("11"),]$freq * cap_hist[ cap_hist$cap_hist %in% c("11"),]$recaps)
R

C/R
R/C

N_hat <- M * C / R
N_hat

R/C # recapture probability


## Now for the hypergeometric using LP_fit

## Not clear how they get n2 (or R) for the number sampled at t2 under hypergeometric sampling.
## We need to remove duplicate juvenile data where the same POP was sampled.
## Rats we need to remove duplicated juvenile data that is NOT in the marked sample, but all of the
## InferredDads and InferredMums are * and # so I can't identify duplicated.

# Use the methods of Arnason, Schwarz and Gerrard

mprime <- sum(cap_hist.red$freq[ cap_hist.red$cap_hist=="11" & cap_hist.red$recaps >0])
mprime
m =     sum(cap_hist.red$freq  [ cap_hist.red$cap_hist=="11" & cap_hist.red$recaps >0] *
            cap_hist.red$recaps[ cap_hist.red$cap_hist=="11" & cap_hist.red$recaps >0])
m
n = C
n

# Solve equation (3) for Mhat
Mhat = M
diff <- function(Mhat){m/Mhat - sum(1/(Mhat-(0:(mprime-1))))}

diff(200)
Mhat <- uniroot(diff, c(200,500))
Mhat

Nhat <- Mhat$root * n / m
Nhat




#M marked adults seen a total of M.recaps times
M.recaps <- sum( cap_hist.red$n[cap_hist.red$cap_hist =="11"])
M.recaps

M/M.recaps
M +M/M.recaps * cap_hist.red$n[ cap_hist.red$cap_hist == "01"]




# I will assume that the same "fraction" of the juvenile without mom/dad identified
# are also duplicated

select <- juveniles$InferredDad != "*" | juveniles$InferredMum != "#" # juvenile that have a dad/mum identified
dup.recaps <- duplicated(juveniles[select,c("InferredDad","InferredMum")])
prop.dup.recaps <- mean(dup.recaps)



cap_hist.red
# Need to adjust the 01 frequence to account for duplicated recaps that must be removed

cap_hist.red$freq[ cap_hist.red$cap_hist== "01"] <- round(cap_hist.red$freq[ cap_hist.red$cap_hist== "01"] * (1-prop.dup.recaps))
cap_hist.red
1124*(1-prop.dup.recaps)

Petersen::LP_summary_stats(cap_hist)
1124+184

