###############################################################################
### Parametric bootstrap goodness of fit using deviance and Tukey Statistic ###
###############################################################################
###############################################################################

# Function: parametric.bootstrap 
# inputs : model information after fiting the model and number of bootstrap 
#          samples needed(n.boots),
# output : Histogram of the bootstrap deviance and gof p-value and histogram
#          for Tukey statistic and p-value

parametric.bootstrap = function(model.info,n.boots=1000){
  
  model.id = model.info$model.id
  
  captureDM = model.info$captureDM
  thetaDM = model.info$thetaDM
  lambdaDM = model.info$lambdaDM
  
  captureOFFSET = model.info$captureOFFSET
  thetaOFFSET = model.info$thetaOFFSET
  lambdaOFFSET = model.info$lambdaOFFSET
  
  ############ Deviance for observed data ########################
  
  # get unique capture histories from all the histories in the data set
  uniq.history = unique(model.info$rawdata$History)  
  #sum of the counts for each of history in Data$History
  obs.counts = outer(uniq.history, model.info$rawdata$History,
                     '==') %*% model.info$rawdata$count 
  
  unique.Data = NULL
  unique.Data$History = uniq.history  #  unique histories
  unique.Data$counts = obs.counts    # total counts for related unique histories
  unique.Data$category = model.info$rawdata$category
  
  ## Likelihood calculation for observed saturated model (for observed data), 
  #   this is conditioning on "n",  the total individuals captured
  obs.satu.model.logL = saturated.model.log.likelihood(unique.Data$counts) 
  
  
  ## Likelihood for fitted model(for observed data)
  obs.fitted.model.log.likelehood = -model.info$NLL
  
  # probability of the capture history "00"
  p_00 = sum(model.info$est$lambda*(1-model.info$est$p[,1])*(1-model.info$est$p[,2])) 
  N = model.info$est$N # estimated population size
  tot.captured = sum(unique.Data$counts)  # total individuals captured
  
  # remove following from the likelihood of the fitted model to get conditional 
  # log likelihood to calculate deviance 
  remove.fitted.log.likelihood = lgamma(N+1)-lgamma(N-tot.captured+1)+
                                  (N-tot.captured)*log(p_00) + tot.captured*log(1-p_00)
  
  ## Likelihood part for the deviance calculation for fitted model(for observed data)
  ## need to add log(n!) to match the deviance part of the saturated model
  obs.fitted.model.conditional.logL = obs.fitted.model.log.likelehood -
                                        remove.fitted.log.likelihood + lgamma(tot.captured+1) 
  
  
  # Deviance for observed data
  Deviance.obs = -2*( obs.fitted.model.conditional.logL) - (-2*( obs.satu.model.logL))
  
  ######### Tukey Statistic value for observed data #####################
  
  tukey.observed.counts = model.info$obs.exp.counts$Observed.Counts
  tukey.expected.counts = model.info$obs.exp.counts$Expected.counts
  
  tukey.statistic.observed.data = sum ((sqrt(tukey.observed.counts) - 
                                          sqrt(tukey.expected.counts))^2)
  
  ############ Deviance for bootstrap data #########################
  
  cats = unique.Data$category
  
  cap.histories = c(paste(cats,"0",sep=""),paste(cats,cats,sep=""),
                    paste("0",cats,sep=""))
  
  boots.Data = NULL
  # all possible observable capture histories("00" is unobservable)
  boots.Data$History = c("U0", "UU", cap.histories,"0U") 
  boots.Data$category = cats
  
  hist = c(boots.Data$History,"00")# add "00" capture history to end 
                                  # of the vector of capture history of expected.data
  
  # the function "log.prob.history" is in the file "neg.log.likelihood"
  # which was used to find the log probabilities
  
  # call the function 'log.prob.history'. 
  expected.log.prob = log.prob.history(model.info$est,hist,cats)  
  # These log probabilities and the counts are in the exact order 
  # related to each capture history. 
  
  expected.prob = exp(expected.log.prob) - (expected.log.prob==0) # expected probabilities
  # expected probabilities for observable capture historis( without history"00")
  obs.hist.exp.prob  =  expected.prob[1:(length(expected.prob)-1)] 
  
  # conditional probabilities for  observable capture histories
  cond.obs.hist.exp.prob = obs.hist.exp.prob/(1-p_00)  
  
  sample.size = sum(unique.Data$counts) # bootstrap sample size

  # devince for  n.boots number of bootstrap samples
  Deviance.boots = rdply(n.boots,parametric.bootstrap.deviance(cond.obs.hist.exp.prob,
                                                               sample.size ,boots.Data)) 
  names(Deviance.boots) = c("sample","DEV","Tukey")
  
  ####### calculate the gof chi-square p-value #####################
  
  p.value = mean(Deviance.boots$DEV >Deviance.obs)
  
  ####### calculate the gof Tukey statistic p-value ################
  
  tukey.p.value = mean(Deviance.boots$Tukey >tukey.statistic.observed.data)

  ####### create Histogram for deviance #########################################
  
  # mark the deviance for observed data
  cuts1 = data.frame(Values=c("Observed Deviance"), vals=c(Deviance.obs)) 
  # x position to show the p-value in the graph
  x.pos = min(Deviance.boots$DEV) + .6*(max(Deviance.boots$DEV,Deviance.obs) - 
                                          min(Deviance.boots$DEV)) 
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
                                     x = x.pos, hjust = 0, y = Inf, vjust = 2, color = "darkred")
  

  ####### create Histogram for Tukey statistic ################################
  
  # mark the deviance for observed data
  cuts2 = data.frame(Values=c("Observed Tukey Statistic"), vals=c(tukey.statistic.observed.data)) 
  # x position to show the p-value in the graph
  x.pos = min(Deviance.boots$Tukey) + .6*(max(Deviance.boots$Tukey,tukey.statistic.observed.data) - 
                                            min(Deviance.boots$Tukey))
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
                                  x = x.pos, hjust = 0, y = Inf, vjust = 1, color = "darkred")
  
  
  
  ## return gof p-value and the the histogram of the bootstrap sampling distribution ##
  result = NULL
  result$p.value = p.value
  result$boots.histogram.deviance  = para.boots.histogram
  #result$boots.histogram.2 = para.boots.histogram.2 
  result$boots.histogram.tukey  = para.boots.histogram.tukey
  
  return(result)
  
} # end of non.parametric.bootstrap 


###############################################################################
########### find the conditional log likelihood for saturated model ###########
###############################################################################

#Function : saturated.model.log.likelihood
#   input  : counts for unique histories ( might have zero entries for counts)
#   output : log likelihood value for saturated model

saturated.model.log.likelihood = function(counts){
  counts = counts[counts!=0] # remove the entries with zero counts
  total = sum(counts) 
  satu.model.logL = lgamma(total+1) - sum(lgamma(counts+1)) + 
                                      sum( counts* log(counts/total)) 
  
  return(satu.model.logL)
  
} # end of saturated.model.log.likelihood 


###############################################################################
## find the conditional log likelihood for fitted model of a bootstrap sample #
###############################################################################

# Function :  fitted.model.conditional.log.likelihood 
#  Inputs: possible observable  histories and categories of a bootstrap sample 
#          (bootstrap.Data) along with model.id,captureDM,thetaDM,lambdaDM,
#           captureOFFSET,thetaOFFSET,lambdaOFFSET
#  output : conditional log likelihood  for the fitted model of a bootstrap 
#            sample and theTukeys statistic for bootstrap sample

fitted.model.conditional.log.likelihood = function(bootstrap.Data,model.id,
                                                   captureDM,thetaDM,lambdaDM,
                                                   captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  ## Likelihood for fitted model for parametric bootstrap data
  boots.fitted.model=  fit.model(model.id,bootstrap.Data,captureDM,thetaDM,lambdaDM,
                                 captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  # probability of the capture history "00"
  p_00 = sum(boots.fitted.model$est$lambda*(1-boots.fitted.model$est$p[,1])*
               (1-boots.fitted.model$est$p[,2])) 
  N = boots.fitted.model$est$N # estimated population size
  tot.captured = sum( bootstrap.Data$counts)  # total individuals captured
  
  # remove following from the likelihood of the fitted model to calculate deviance
  remove.fitted.log.likelihood = lgamma(N+1) - lgamma(N-tot.captured+1) +
                                   (N-tot.captured)*log(p_00) + tot.captured*log(1-p_00)
  
  ## conditional likelihood part for the deviance calculation for
  ## fitted model(for bootstrap data)
  boots.fitted.model.conditional.logL = -(boots.fitted.model$NLL) - 
                                  remove.fitted.log.likelihood +lgamma(tot.captured+1)
  
  ######### Tukey Statistic for bootstrap data #####################
  
  tukey.boots.observed.counts = boots.fitted.model$obs.exp.counts$Observed.Counts
  tukey.boots.expected.counts = boots.fitted.model$obs.exp.counts$Expected.counts
  
  tukey.statistic.boots.data = sum ((sqrt(tukey.boots.observed.counts) - 
                                       sqrt(tukey.boots.expected.counts))^2)
  
  return(c(boots.fitted.model.conditional.logL,tukey.statistic.boots.data) )
  
} # end of fitted.model.log.likelehood


###############################################################################
##########  bootstrap deviance function for one replicate #####################
###############################################################################

# Function : parametric.bootstrap.deviance
# input :  observed data for unique capture histaries
# output : deviance for a bootstrap sample, and Tukey Statistic for bootstrap sample

parametric.bootstrap.deviance = function(expected.prob,sample.size, boots.Data){
  
  bootstrap.Data = boots.Data
  
  # generate bootstrap sample for observable capture histories ( without "00")
  bootstrap.Data$counts = rmultinom(1, size = sample.size, prob = expected.prob) 
  
  ## Likelihood for saturated model for bootstrap data
  boots.satu.model.log.likelihood = saturated.model.log.likelihood(bootstrap.Data$counts)
  
  ## conditional likelihood for fitted model for bootstrap data to calculate deviance
  result = fitted.model.conditional.log.likelihood(bootstrap.Data,model.id=model.id,
                                                   captureDM=captureDM,
                                                   thetaDM=thetaDM,
                                                   lambdaDM=lambdaDM,
                                                   captureOFFSET=captureOFFSET,
                                                   thetaOFFSET=thetaOFFSET,
                                                   lambdaOFFSET=lambdaOFFSET)
  ## conditional likelihood for fitted model for bootstrap data to calculate deviance
  boots.fitted.model.conditional.LogL = result[1] 
  
  tukey.statistic.boots.data = result[2] # Tukey Statistic for bootstrap data 
  
  # devince for  parametric bootstrap sample
  Dev.boots = -2*boots.fitted.model.conditional.LogL - 
                                      (-2* boots.satu.model.log.likelihood) 
  
  return(c(Dev.boots,tukey.statistic.boots.data))
}

###############################################################################
###############################################################################
########## Testing parametric.bootstrap function ##############################
###############################################################################

# # source("load.R")  # load required functions and packages
# # 
#   test.parametric.bootstrap = function(Data){
#     
#   #model identification : unrestricted model
#   model.id = paste("( p(c*t), theta(t), lambda(c) )")
#   
#   # give  the required design matrices
#   captureDM = create.DM(c(1,2,3,4)) # Design matrix for capture prob
#   thetaDM   = create.DM(c(1,2)) # Design matrix for theta
#   lambdaDM  = create.DM(c(1))  # Design matrix for lambda
#   
#   #give the offset vectors(vectors of zeros should be given since no restriction)
#   captureOFFSET = c(0,0,0,0) 
#   thetaOFFSET   = c(0,0)
#   lambdaOFFSET  = c(0)
#   
#   n.boots =10 # enter the number of bootstrap samples needed 
#   
#   test.model = fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
#                          captureOFFSET,thetaOFFSET,lambdaOFFSET)
#    
#    res = parametric.bootstrap(test.model,n.boots)
#    return(res)
#   }
#  # 
#   Data = get.data("Mille_Lacs_Walleye.csv") # get  Mille Lacs Walleye dataset 
#  test.output =  test.parametric.bootstrap(Data)
#  test.output$boots.histogram.deviance
#  test.output$boots.histogram.tukey
###############################################################################
###############################################################################

