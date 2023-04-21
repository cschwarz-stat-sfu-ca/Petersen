############################################################################### 
#####                Parametric bootstrap goodness of fit                ######
#####                 using deviance and Tukey Statistic                 ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
############################################################################### 

## Functions in this file
##   PSkSCR.parametric.bootstrap
##   PSkSCR.saturated.model.log.likelihood 
##   PSkSCR.parametric.bootstrap.deviance
##   PSkSCR.fitted.model.conditional.log.likelihood 
##   PSkSCR.fit.model.bootstrap 
##   PSkSCR.neg.log.likelihood.bootstrap
##   PSkSCR.log.prob.history.bootstrap


###############################################################################
###############################################################################

# Function: PSkSCR.parametric.bootstrap 
# inputs : model information after fiting the model and number of bootstrap 
#          samples needed(n.boots),
# output : Histogram of the bootstrap deviance and gof p-value and histogram
#          for Tukey statistic and p-value


PSkSCR.parametric.bootstrap = function(model.info,n.boots=100){
  
  # model.info <- MLE_PSkSCR_model_1 # to be deleted
  model.id <- model.info$model.id
  
  captureDM <- model.info$captureDM
  thetaDM <- model.info$thetaDM
  lambdaDM <- model.info$lambdaDM
  p_lossDM <- model.info$p_lossDM
  
  captureOFFSET <- model.info$captureOFFSET
  thetaOFFSET <- model.info$thetaOFFSET
  lambdaOFFSET <- model.info$lambdaOFFSET
  p_lossOFFSET <- model.info$p_lossOFFSET
  
  indicator <- model.info$indicator
  parm <- model.info$est
  
  ############ Deviance for observed data ########################
  
  # unique capture histories from all the histories in the data set is needed. 
  # Since dat is prepared within the "PSkSCR.fit.model" function each capture
  # history can have maximum of two dupticates. One to represent number of 
  # alive animals(positive counts)  and otherone represents the number of 
  # dead animals(negative counts). 
  
  unique.data <- model.info$rawdata
  
  ## Likelihood calculation for observed saturated model (for observed data), 
  #   this is conditioning on "n",  the total individuals captured
  obs.satu.model.logL <- PSkSCR.saturated.model.log.likelihood(unique.data$counts) 
  
  
  ######### Likelihood for fitted model(for observed data) ######
  
  obs.fitted.model.log.likelehood <- (-1)*model.info$NLL
  
  # probability for capture history of unseen animals
  # (i.e. for two sampling occasions "00", for 3 sampling occasion "000" , ..)
  p_00 <- sum(apply((1-model.info$est$p),1,prod) * model.info$est$lambda)
  
  N <- model.info$est$N # estimated population size
  tot.captured <- sum(abs(unique.data$counts))  # total individuals captured
  
  # remove following from the likelihood of the fitted model to get conditional 
  # log likelihood to calculate deviance 
  remove.fitted.log.likelihood <- lgamma(N+1)-lgamma(N-tot.captured+1)+
                                              (N-tot.captured)*log(p_00) + 
                                              tot.captured*log(1-p_00)
  
  ## Likelihood part for the deviance calculation for fitted model(for observed data)
  ## need to add log(n!) to match the deviance part of the saturated model
  obs.fitted.model.conditional.logL <- obs.fitted.model.log.likelehood -
                                       remove.fitted.log.likelihood + 
                                       lgamma(tot.captured+1) 
  
  
  # Deviance for observed data
  Deviance.obs <- -2*( obs.fitted.model.conditional.logL) - (-2*( obs.satu.model.logL))
  
  ######### Tukey Statistic value for observed data #####################
  
  tukey.observed.counts <- model.info$obs.exp.counts$Observed.counts
  tukey.expected.counts <- model.info$obs.exp.counts$Expected.counts
  
  tukey.statistic.observed.data <- sum ((sqrt(abs(tukey.observed.counts)) - 
                                           sqrt(abs(tukey.expected.counts)))^2)
  
  
  ############ Deviance for bootstrap data #########################
  
  cats <- unique.data$category
  
  boots.data <- NULL
  # all possible observable capture histories("00" is unobservable)
  boots.data$history <- c(unique.data$history,"OTHER") 
  boots.data$category <- cats
  boots.data$sign <- sign(model.info$obs.exp.counts$Expected.counts)
  
  parm <- model.info$est
  indicator <- model.info$indicator
 
  expected.log.prob  <- PSkSCR.log.prob.history(parm,indicator)
  # probabilities of the all observed capture histories and the "Unseen Capture History"
  expected.prob  <- exp(expected.log.prob) - (expected.log.prob==0)
  expected.prob.OTHER <- 1- sum(expected.prob)
  
  # Probability of all possible observable capture histories. 
  # Here, there is a combined capture history "OTHER"
  obs.hist.exp.prob <- c( expected.prob[1:length(expected.prob)-1], expected.prob.OTHER)
  
  # conditional probabilities for  observable capture histories
  cond.obs.hist.exp.prob <- obs.hist.exp.prob/(1-p_00)  
  
  sample.size <- sum(abs(unique.data$counts)) # bootstrap sample size
  
  # devince for  n.boots number of bootstrap samples
  Deviance.boots <- rdply(n.boots,
                          PSkSCR.parametric.bootstrap.deviance(cond.obs.hist.exp.prob,
                                                               sample.size ,
                                                               boots.data)) 
  names(Deviance.boots) = c("sample","DEV","Tukey")
  
  ####### calculate the gof chi-square p-value #####################
  
  p.value = mean(Deviance.boots$DEV >Deviance.obs)
  
  ####### calculate the gof Tukey statistic p-value ################
  
  tukey.p.value = mean(Deviance.boots$Tukey >tukey.statistic.observed.data)
  
  ####### create Histogram for deviance #########################################
  
  # mark the deviance for observed data
  cuts1 = data.frame(Values=c("Observed Deviance"), vals=c(Deviance.obs)) 
  # x position to show the p-value in the graph
  x.pos = min(Deviance.boots$DEV,Deviance.obs) + 0.3*(max(Deviance.boots$DEV,Deviance.obs) - 
                                          min(Deviance.boots$DEV,Deviance.obs)) 
  
  para.boots.histogram = ggplot(Deviance.boots, aes(x=DEV,..density.. ))+   
    geom_histogram(aes(y=..density.. ))+
    geom_vline(data=cuts1,aes(xintercept=vals,
                              linetype=Values,colour=Values),
               size=1,show_guide = FALSE)+
    xlab("Deviance") + ylab("Density")+
    ggtitle(paste("Parametric Bootstrap Sampling Distribution for Deviance \n model: ",
                  model.id, sep=" " ))+
    theme(axis.text=element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"), 
          plot.title = element_text(size = 14,face="bold"))+
    annotate("text", label = paste("Observed Deviance = ",
                                   round(Deviance.obs,2),"\np-value = ",
                                   p.value, sep=""),
             x = x.pos, hjust = 0, y = Inf, vjust = 1.5, color = "darkred")
  
  
  ####### create Histogram for Tukey statistic ################################
  
  # mark the deviance for observed data
  cuts2 = data.frame(Values=c("Observed Tukey Statistic"), vals=c(tukey.statistic.observed.data)) 
  # x position to show the p-value in the graph
  x.pos = min(Deviance.boots$Tukey,tukey.statistic.observed.data) + 
               0.3*(max(Deviance.boots$Tukey,tukey.statistic.observed.data) - 
                              min(Deviance.boots$Tukey,tukey.statistic.observed.data))
  
  para.boots.histogram.tukey = ggplot(Deviance.boots, aes(x=Tukey,..density.. ))+   
    geom_histogram(aes(y=..density.. ))+
    geom_vline(data=cuts2,aes(xintercept=vals,linetype=Values,
                              colour=Values),size=1,
               show_guide = FALSE)+
    xlab("Tukey Statistic") + ylab("Density")+
    ggtitle(paste("Parametric Bootstrap Sampling Distribution for Tukey Statistic \n model: ",
                  model.id, sep=" " ))+
    theme(axis.text=element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"), 
          plot.title = element_text(size = 14,face="bold"))+
    annotate("text", label = paste("Observed Tukey Statistic = ",
                                   round(tukey.statistic.observed.data,2),
                                   "\np-value = ", tukey.p.value, sep=""), 
             x = x.pos, hjust = 0, y = Inf, vjust = 1.5, color = "darkred")
  
  
  
  ## return gof p-value and the the histogram of the bootstrap sampling distribution ##
  result = NULL
  result$Chi.square.DEvince.p.value = p.value
  result$tukey.p.value <- tukey.p.value
  result$boots.histogram.deviance  = para.boots.histogram
  #result$boots.histogram.2 = para.boots.histogram.2 
  result$boots.histogram.tukey  = para.boots.histogram.tukey
  
  return(result)
  
  
} # end of "PSkSCR.parametric.bootstrap"



###############################################################################
########### find the conditional log likelihood for saturated model ###########
###############################################################################

#Function : PSkSCR.saturated.model.log.likelihood
#   input  : counts for unique histories ( might have zero entries for counts)
#   output : log likelihood value for saturated model

PSkSCR.saturated.model.log.likelihood <- function(counts){
  counts <- counts[counts!=0] # remove the entries with zero counts
  total <- sum(abs(counts))
  satu.model.logL <- lgamma(total+1) - sum(lgamma(abs(counts)+1)) + 
    sum( abs(counts)* log(abs(counts)/total)) 
  
  return(satu.model.logL)
  
} # end of PSkSCR.saturated.model.log.likelihood 


###############################################################################
##########  bootstrap deviance function for one replicate #####################
###############################################################################

# Function : parametric.bootstrap.deviance
# input :  observed data for unique capture histaries
# output : deviance for a bootstrap sample, and Tukey Statistic for bootstrap sample

PSkSCR.parametric.bootstrap.deviance = function(cond.obs.hist.exp.prob,sample.size, boots.data){
  
  bootstrap.data <- NULL
  
  # generate bootstrap sample for observable capture histories ( without "00")
  bootstrap.sample.counts <- boots.data$sign * rmultinom(1, size = sample.size, 
                                                         prob = cond.obs.hist.exp.prob) 
  
  # counts of group of observable histories,  but not observed histories
  other.count <- bootstrap.sample.counts[length(bootstrap.sample.counts)]
  
  length.data <- length(boots.data$history)
  bootstrap.data$history <- boots.data$history[1:(length.data-1)]
  bootstrap.data$counts  <- bootstrap.sample.counts[1:(length.data-1)]
  bootstrap.data$category <- boots.data$category
  
  
  ## Likelihood for saturated model for bootstrap data
  boots.satu.model.log.likelihood <- PSkSCR.saturated.model.log.likelihood(bootstrap.sample.counts)
  
  ## conditional likelihood for fitted model for bootstrap data to calculate deviance
  result = PSkSCR.fitted.model.conditional.log.likelihood(bootstrap.data,
                                                          model.id=model.id,
                                                          captureDM=captureDM,
                                                          thetaDM=thetaDM,
                                                          lambdaDM=lambdaDM,
                                                          p_lossDM=p_lossDM,
                                                          captureOFFSET=captureOFFSET,
                                                          thetaOFFSET=thetaOFFSET,
                                                          lambdaOFFSET=lambdaOFFSET,
                                                          p_lossOFFSET=p_lossOFFSET,
                                                          other.count=other.count)
  ## conditional likelihood for fitted model for bootstrap data to calculate deviance
  boots.fitted.model.conditional.LogL = result[1] 
  
  tukey.statistic.boots.data = result[2] # Tukey Statistic for bootstrap data 
  
  # devince for  parametric bootstrap sample
  Dev.boots = -2*boots.fitted.model.conditional.LogL - 
    (-2* boots.satu.model.log.likelihood) 
  
  return(c(Dev.boots,tukey.statistic.boots.data))
}


###############################################################################
## find the conditional log likelihood for fitted model of a bootstrap sample #
###############################################################################

# Function :  PSkSCR.fitted.model.conditional.log.likelihood 
#  Inputs: possible observable  histories and categories of a bootstrap sample 
#          (bootstrap.Data) along with model.id,captureDM,thetaDM,lambdaDM,
#           captureOFFSET,thetaOFFSET,lambdaOFFSET
#  output : conditional log likelihood  for the fitted model of a bootstrap 
#            sample and theTukeys statistic for bootstrap sample

PSkSCR.fitted.model.conditional.log.likelihood <- function(bootstrap.data,
                                                           model.id,
                                                           captureDM,
                                                           thetaDM,
                                                           lambdaDM,
                                                           p_lossDM,
                                                           captureOFFSET,
                                                           thetaOFFSET,
                                                           lambdaOFFSET,
                                                           p_lossOFFSET,
                                                           other.count){
  
  ## Likelihood for fitted model for parametric bootstrap data
  boots.fitted.model<- PSkSCR.fit.model.bootstrap(model.id,bootstrap.data,
                                                  captureDM,thetaDM,
                                                  lambdaDM,p_lossDM,
                                                  captureOFFSET,thetaOFFSET,
                                                  lambdaOFFSET,p_lossOFFSET,
                                                  other.count)
 
  # probability of the capture history "00"
  p_00 <- sum(apply((1-boots.fitted.model$est$p),1,prod) * boots.fitted.model$est$lambda)
  
  N <- boots.fitted.model$est$N # estimated population size
  tot.captured <- sum( abs(bootstrap.data$counts)) + other.count  # total individuals captured
  
  # remove following from the likelihood of the fitted model to calculate deviance
  remove.fitted.log.likelihood <- lgamma(N+1) - lgamma(N-tot.captured+1)+
                                  (N-tot.captured)*log(p_00) + tot.captured*log(1-p_00)
  
  ## conditional likelihood part for the deviance calculation for
  ## fitted model(for bootstrap data)
  boots.fitted.model.conditional.logL <- -(boots.fitted.model$NLL) - 
                                          remove.fitted.log.likelihood +lgamma(tot.captured+1)
  
  ######### Tukey Statistic for bootstrap data #####################
  tukey.boots.observed.counts <- abs(boots.fitted.model$obs.exp.counts$Observed.counts)
  tukey.boots.expected.counts <- abs(boots.fitted.model$obs.exp.counts$Expected.counts)
  
  tukey.statistic.boots.data <- sum ((sqrt(tukey.boots.observed.counts) - 
                                       sqrt(tukey.boots.expected.counts))^2)
  
  return(c(boots.fitted.model.conditional.logL,tukey.statistic.boots.data) )
  
} # end of fitted.model.log.likelehood


###############################################################################
####              Fit Model  for parametric bootstrap                     #####
#### There is extra History in the bootstrap data set compared to rawdata #####
####  possible, but uncaptured historied combined in to "other" history   #####
#### therefore CANNOT use the function "PSkSCR.fit.model"                 #####

PSkSCR.fit.model.bootstrap <- function(model.id ,bootstrap.data,
                                       captureDM, thetaDM,lambdaDM,p_lossDM,
                                       captureOFFSET,thetaOFFSET,
                                       lambdaOFFSET,p_lossOFFSET,
                                       other.count,
                                       dataprepare="NO"){
  
  ## dataprepare is similar in the function"PSkSCR.fit.model"
  ## dataprepare is needed to remove  "0" counts in bootstrap data
  
  if(dataprepare=="NO"){
    data <- bootstrap.data
    # histories of the given data
    hist <- data$history 
    # corresponding counts for the  histories of the given data
    counts <- data$counts 
    # replicate each history to tally the corresponding counts
    individual.hist <- rep(hist,abs(counts))
    # create -1 an +1 to tally the corresponding counts
    sign.counts <- rep( sign(counts),abs(counts))
    
    # create a two-way table with unique histories and capture indicators
    # (eg. Capture histories might duplicate two times if there are both positive
    #      and negative  counts for each history )
    summary <-  table(individual.hist,sign.counts)
    
    data.summary <- as.data.frame(summary) # create a data frame
    names(data.summary) <- c("history", "cap.ind", "frequency")
    
    new.cap.ind <- rep(0, nrow(data.summary))
    new.cap.ind[ data.summary$cap.ind==-1] <- -1
    new.cap.ind[ data.summary$cap.ind!=-1] <- 1
    
    data.summary$counts <- new.cap.ind* data.summary$frequency
    
    obs.summary <- data.summary[data.summary$counts != 0, ]
    
    data$history <- as.character(obs.summary$history)
    data$counts <- obs.summary$counts
  }
  
  
  ###### End Preparing data ##############
  ########################################
  
  
  Result <- NULL
  ncats <- length(data$category) # number of categories
  st <-  nchar(data$history[1]) # number of sampling occasions
  Result$st <- st
  
  Result$rawdata <- data # data set used for bootstrap. No "0" counts here
  
  # create required indicator variables in matrix form. These are required for 
  # calculating the likelihood
  indicator <- PSkSCR.create.indicator(data)
  Result$indicator <- indicator
  
  # Initial Estimates (in regular form)
  initial.est <- PSkSCR.initial.estimates(data,indicator)
  
  # initial.pack.parm.est are in logit/log form
  initial.pack.parm.est <- PSkSCR.pack.parm(initial.est$full,data,
                                            captureDM,thetaDM,lambdaDM,p_lossDM)  
  Result$initial.unpack.parm <- PSkSCR.unpack.parm(initial.pack.parm.est,data, 
                                                   captureDM, thetaDM,lambdaDM,p_lossDM,
                                                   captureOFFSET,thetaOFFSET,
                                                   lambdaOFFSET,p_lossOFFSET)

  
  # initial value for the log(N)
  log.N.init <- initial.pack.parm.est[length(initial.pack.parm.est)]  
  
  # lower bound for optimization routine
  l.b <- c(rep(logit(0.000001),length(initial.pack.parm.est)-1),max(2,log.N.init-2))
  # upper bound for optimization routine
  u.b <- c(rep(logit(.99999),length(initial.pack.parm.est)-1),log.N.init+2) 
  
  # optimize the negative log likelihood function(MLEs in logit/log form)
  res <- optim(par=initial.pack.parm.est, fn=PSkSCR.neg.log.likelihood.bootstrap, 
               method="L-BFGS-B",
               lower=l.b, upper=u.b , control=list(trace=3,maxit=1000),
               data=data, indicator=indicator, 
               captureDM=captureDM, thetaDM=thetaDM,lambdaDM=lambdaDM,p_lossDM=p_lossDM,
               captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
               lambdaOFFSET=lambdaOFFSET, p_lossOFFSET=p_lossOFFSET,
               other.count=other.count)
  
  Result$res <- res
  
  #extract the MLEs for the parameter (estimates in regular form)
  Result$est <- PSkSCR.unpack.parm(Result$res$par,data,
                                   captureDM,thetaDM,lambdaDM,p_lossDM,
                                   captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)  
  
  #Negative log likelihood value
  Result$NLL = PSkSCR.neg.log.likelihood.bootstrap(Result$res$par,data,indicator,
                                                   captureDM,thetaDM,lambdaDM,p_lossDM,
                                                   captureOFFSET,thetaOFFSET,
                                                   lambdaOFFSET,p_lossOFFSET,
                                                   other.count )
  
  # expected log probabilities
  expected.log.prob  <- PSkSCR.log.prob.history.bootstrap(Result$est,indicator)
  
  # expected probabilities
  expected.prob = exp(expected.log.prob) - (expected.log.prob==0) 
  
  # observable Histories
  data.histories <- c(Result$rawdata$history, "Other" )
  # counts bor all observable Histories
  observed.data.counts <- c( Result$rawdata$counts, other.count)
  
  N <- Result$est$N # MLE for Population size
  # expected values for observable Histories
  expected.data.counts <- N * expected.prob[1:(length(expected.prob)-1)]*sign(observed.data.counts)
  
  output.df <- data.frame(History = data.histories,
                         Observed.counts = observed.data.counts,
                         Expected.counts = expected.data.counts)
  
  # capture Histories and obeserved and expected counts
  Result$obs.exp.counts <- output.df
  
  
  
  return(Result)
  
  
} # end of the function "PSkSCR.fit.model.bootstrap"


###############################################################################

# The function "PSkSCR.neg.log.likelihood.bootstrap"  produce the negative 
# log likelihood value considering all parameters(here it consider the p_loss too)
# and considerinf allthe capture histories incluting the history "other"


PSkSCR.neg.log.likelihood.bootstrap <- function(logit.est, data,indicator,
                                      captureDM,thetaDM,lambdaDM,p_lossDM,
                                      captureOFFSET,thetaOFFSET,
                                      lambdaOFFSET, p_lossOFFSET, 
                                      other.count){
  
  parm <- PSkSCR.unpack.parm(logit.est,data,
                             captureDM,thetaDM,lambdaDM,p_lossDM,
                             captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
  
  N <- parm$N
  
  hist <- data$history # obsevable capture histories from data except the histories "other" 
  abs.ct   <- abs(data$counts)   # absolute values of counts from raw data
  other.ct <- abs(other.count)   #absolute value of the count og the capture history "other"
  ct.00 <- (N-sum(abs.ct)-other.ct ) # count for capture history "00"
  # counts of observable capture histories and capture history "00"
  ct <- c(abs.ct,other.ct,ct.00)
  
  
  # call the function 'PSkSCR.log.prob.history'. These probabilities and the counts
  # are in the exact order related to each observable capture history.
  # probabilities of observable capture histories and capture history of unseen animals
  # (i.e. for two sampling occasions "00", for 3 sampling occasion "000" , ..)
  log.prob.hist <- PSkSCR.log.prob.history.bootstrap(parm,indicator) 
  
  #######################################
  # log-ikelihood has 2 parts
  part1 <- lgamma(N+1) - sum(lgamma(ct+1)) # factorial part
  
  part2 <- sum(ct * log.prob.hist)
  
  # log-likelihood
  LLH <-  part1 +  part2  
  
  # negative log-likelihood
  NLLH <- - LLH  # -(log-likelihood),convert to negative value for optimization purpose
  
  return(NLLH)
  
} # end of "PSkSCR.neg.log.likelihood.bootstrap "



###############################################################################
####### This function is used only for parametric bootstrap step ##############
############################################################################### 

# histories are passed through indicator variables (matrices)

PSkSCR.log.prob.history.bootstrap <- function(parm,indicator){
  
        h <- indicator$h
        z <- indicator$z
  cat.ind <- indicator$cat.ind
        s <- indicator$s
       cs <- indicator$cs
    s.ind <- indicator$s.ind
     dead <- indicator$dead
  
       p <- parm$p
  lambda <- parm$lambda
   theta <- parm$theta
  p_loss <- parm$p_loss
  
  
  
  prob_hist<- rep(0, nrow(h))
  for( j in 1:nrow(h) ){
    for(i in 1:length(lambda)){
      temp <- lambda[i]*prod( (p[i,]^(h[j,]*z[j,])) * 
                                ((1-p_loss)^(h[j,]*dead[j,])) * 
                                (theta^(s[j,]))* 
                                ((1-p[i,])^((1-h[j,])*z[j,])) * 
                                (p_loss^(h[j,]*(1-dead[j,]))) *
                                ((1-theta)^((1-s[j,])*cs[j,])) ) * s.ind[j,i]
      prob_hist[j] <- temp + prob_hist[j] 
    }
  }
  
  
  # probability of observed capture histories
  prob.hist.observed <- prob_hist
  
  # probability for capture history of unseen animals
  # (i.e. for two sampling occasions "00", for 3 sampling occasion "000" , ..)
  prob.hist.00 <- sum(apply((1-p),1,prod) *lambda)
  
  # probability of the history "other"
  prob.other <- 1-( sum(prob.hist.observed) + prob.hist.00 )
  
  if(prob.other < 0){
    print(" Error in probability of histories calculations in bootstrap")
    stop
  }
  
  # probabilities of all observable capture histories including 
  # the history "other" and capture history "00"
  prob.hist <- c(prob.hist.observed, prob.other, prob.hist.00)
  
  # log- probabilities of observable capture histories and capture history "00"
  log_prob.hist <- log(prob.hist + (prob.hist==0))
  
  return(log_prob.hist)
  
  
} # end of function "PSkSCR.log.prob.history.bootstrap"

###############################################################################
###############################################################################

