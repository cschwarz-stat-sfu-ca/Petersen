########################################################################################################
######### non parametric bootstrap goodness of fit using deviance ######################################
########################################################################################################
#########################################################################################################

# Function: non.parametric.bootstrap 
# inputs : number of bootstrap samples needed(n.boots),
# and model information after fiting the model
# output : Histogram of the bootstrap deviance and gof p-value

non.parametric.bootstrap = function(model.info,n.boots){
  
   #model.info =test.model
   model.id = model.info$model.id
  
   captureDM = model.info$captureDM
     thetaDM = model.info$thetaDM
    lambdaDM = model.info$lambdaDM
  
   captureOFFSET = model.info$captureOFFSET
     thetaOFFSET = model.info$thetaOFFSET
    lambdaOFFSET = model.info$lambdaOFFSET
  
   ############ Deviance for observed data ########################
  
   uniq.history = unique(model.info$rawdata$History)  # get unique capture histories from all the histories in the data set
   obs.counts = outer(uniq.history, model.info$rawdata$History, '==') %*% model.info$rawdata$count #sum of the counts for each of history in Data$History
  
            unique.Data = NULL
    unique.Data$History = uniq.history  #  unique histories
     unique.Data$counts = obs.counts    # total counts for related unique histories
   unique.Data$category = model.info$rawdata$category
  
   ## Likelihood calculation for observed saturated model (for observed data), this is conditioning on "n", 
   # the total individuals captured
   obs.satu.model.logL = saturated.model.log.likelihood(unique.Data$counts) 
  
  
   ## Likelihood for fitted model(for observed data)
   obs.fitted.model.log.likelehood = -model.info$NLL
  
   p_00 = sum(model.info$est$lambda*(1-model.info$est$p[,1])*(1-model.info$est$p[,2])) # probability of the capture history "00"
      N = model.info$est$N # estimated population size
   tot.captured = sum(unique.Data$counts)  # total individuals captured
  
   # remove following from the likelihood of the fitted model to get conditional log likelihood to calculate deviance 
   remove.fitted.log.likelihood = lgamma(N+1)-lgamma(N-tot.captured+1)+(N-tot.captured)*log(p_00) + tot.captured*log(1-p_00)
  
   ## Likelihood part for the deviance calculation for fitted model(for observed data)
   ## need to add log(n!) to match the deviance part of the saturated model
   obs.fitted.model.conditional.logL = obs.fitted.model.log.likelehood - remove.fitted.log.likelihood + lgamma(tot.captured+1) 
  
  
   # Deviance for observed data
   Deviance.obs = -2*( obs.fitted.model.conditional.logL) - (-2*( obs.satu.model.logL))
  
  
  
   ############ Deviance for bootstrap data #########################
   
   # devince for  n.boots number of bootstrap samples
   
   Deviance.boots = rdply(n.boots,bootstrap.deviance(unique.Data)) 
   names(Deviance.boots) = c("sample","DEV")
  
   ####### calculate the gof p-value ################################
   
   p.value = mean(Deviance.boots["DEV"]>Deviance.obs)
   
   ####### create Histogram #########################################
 
   cuts1 = data.frame(Values=c("Observed Deviance"), vals=c(Deviance.obs)) # mark the deviance for observed data
   x.pos = min(Deviance.boots$DEV) + .85*(max(Deviance.boots$DEV) - min(Deviance.boots$DEV)) # x position to show the p-value in the graph
   boots.histogram = ggplot(Deviance.boots, aes(x=DEV,..density.. ))+   
                         geom_histogram(aes(y=..density.. ))+
                            geom_vline(data=cuts1,aes(xintercept=vals,linetype=Values,colour=Values),size=1,show_guide = TRUE)+
                              labs(x = "Deviance",y = "Density",title = paste("Non parametric bootstrap sampling Distribution for Deviance \n model : ", model.id, "\n \n p-value = ", p.value,sep=" "  ))+ 
                                theme(axis.text=element_text(size=10,face="bold"),axis.title=element_text(size=14,face="bold"), plot.title = element_text(size = 14,face="bold"))
   
                                  # annotate("text", label = paste("gof p-value = ", p.value,sep=""), x =x.pos , y = .1) # if need to sow inside the graph
   
   
   #### return gof p-value and the the histogram of the bootstrap sampling distribution ############ 
   result = NULL
   result$p.value = p.value
   result$boots.histogram  = boots.histogram
  
  
   return(result)
  
} # end of non.parametric.bootstrap 


######################################################################################################
########### find the conditional log likelihood for saturated model ##################################
######################################################################################################

#Function : saturated.model.log.likelihood
#   input : counts for unique histories ( might have zero entries for counts)
#   output : log likelihood value for saturated model

saturated.model.log.likelihood = function(counts){
  counts = counts[counts!=0] # remove the entries with zero counts
  total = sum(counts) 
  satu.model.logL = lgamma(total+1) - sum(lgamma(counts+1)) + sum( counts* log(counts/total)) 
  
  return(satu.model.logL)
  
} # end of saturated.model.log.likelihood 


######################################################################################################
########### find the conditional log likelihood for fitted model of a bootstrap sample ###############
######################################################################################################

# Function :  fitted.model.conditional.log.likelihood 
#  Inputs: counts, histories, categories of a bootsttrap sample along with model.id,captureDM,thetaDM,lambdaDM,captureOFFSET,thetaOFFSET,lambdaOFFSET
#  output : conditional log likelihood  for the fitted model of a bootstrap sample

fitted.model.conditional.log.likelihood = function(counts,History, category,model.id,captureDM,thetaDM,lambdaDM,captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  boots.Data = NULL
  boots.Data$History = History
  boots.Data$counts = counts
  boots.Data$category = category
  
  ## Likelihood for fitted model for bootstrap data
  boots.fitted.model=  fit.model(model.id,boots.Data,captureDM,thetaDM,lambdaDM,captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  p_00 = sum(boots.fitted.model$est$lambda*(1-boots.fitted.model$est$p[,1])*(1-boots.fitted.model$est$p[,2])) # probability of the capture history "00"
  N = boots.fitted.model$est$N # estimated population size
  tot.captured = sum( boots.Data$counts)  # total individuals captured
  
  # remove following from the likelihood of the fitted model to calculate deviance
  remove.fitted.log.likelihood = lgamma(N+1) - lgamma(N-tot.captured+1) + (N-tot.captured)*log(p_00) + tot.captured*log(1-p_00)
  
  ## conditional likelihood part for the deviance calculation for fitted model(for bootstrap data)
  boots.fitted.model.conditional.logL = -(boots.fitted.model$NLL) - remove.fitted.log.likelihood +lgamma(tot.captured+1)
  
  return(boots.fitted.model.conditional.logL)
  
} # end of fitted.model.log.likelehood


######################################################################################################
##########  bootstrap deviance function for one replicate ############################################
######################################################################################################

# Function : bootstrap.deviance
# input :  observed data for unique capture histaries
# output : deviance for a bootstrap sample

bootstrap.deviance = function(unique.Data){
  
  tot.captured = sum(unique.Data$counts)  # total number of individuals captured
  prob.counts = unique.Data$counts/tot.captured # probability of individual in unique capture history
  
  boots.sample = rmultinom(1, size = tot.captured, prob = prob.counts) # generate a  bootstrap sample
  
  ## Likelihood for saturated model for bootstrap data
  boots.satu.model.log.likelihood = saturated.model.log.likelihood(boots.sample)
  
  ## conditional likelihood for fitted model for bootstrap data to calculate deviance
  boots.fitted.model.conditional.LogL = fitted.model.conditional.log.likelihood(counts=boots.sample,History=unique.Data$History,
                                                  category = Data$category,model.id=model.id,captureDM=captureDM,thetaDM=thetaDM,
                                                   lambdaDM=lambdaDM,captureOFFSET=captureOFFSET,thetaOFFSET=thetaOFFSET,
                                                    lambdaOFFSET=lambdaOFFSET)
  # devince for  bootstrap sample
  Dev.boots = -2*boots.fitted.model.conditional.LogL - (-2* boots.satu.model.log.likelihood) 
  
  return(Dev.boots)
}

######################################################################################################
######################################################################################################
########## Testing non.parametric.bootstrap function #################################################
######################################################################################################

# source("load.R")  # load required functions and packages
# 
# test.non.parametric.bootstrap = function(Data){
#   
#   #model identification : unrestricted model
#   model.id = paste("( p(c*t), theta(t), lambda(c) )")
#   
#   # give  the required design matrices
#   captureDM = create.DM(c(1,2,3,4)) # Design matrix for capture recapture probabilities
#   thetaDM   = create.DM(c(1,2)) # Design matrix for theta(sampling(sexing) fractions)
#   lambdaDM  = create.DM(c(1))  # Design matrix for lambda(Category proportion)
#   
#   #give the offset vectors(vectors of zeros should be given since no restriction)
#   captureOFFSET = c(0,0,0,0) 
#   thetaOFFSET   = c(0,0)
#   lambdaOFFSET  = c(0)
#   
#   n.boots =100 # enter the number of bootstrap samples needed 
#   
#   test.model = fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,thetaOFFSET,lambdaOFFSET)
#   
#   res = non.parametric.bootstrap(test.model,n.boots)
#   return(res)
# }
# 
# Data = get.data("Mille_Lacs_Walleye.csv") # get  Mille Lacs Walleye dataset 
# test.non.parametric.bootstrap(Data)

######################################################################################################
