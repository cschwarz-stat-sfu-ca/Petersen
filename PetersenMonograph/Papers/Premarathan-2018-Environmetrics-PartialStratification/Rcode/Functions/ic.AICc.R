###############################################################################
############## AIC for the data with individual covariates ####################
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
############################################################################### 
# Conditional likelihood value is used for calculating the AIC value.
# This is okay Since all the models use conditional likelihood values for 
# model comparison
#
# function "ic.AICc" calculates the AIC with a correction for finite sample sizes
#
# AIC = 2*K - 2*log(likelihood) where k is the number of parameters in the model
#    
# AICc = AIC + 2*K*(k+1)/(n-k-1)  ; n denotes the sample size
#
# Input : vector of estimate in logit/log form and Data
#      
# Output: AICc Value. 

ic.AICc <- function(logit.est, Data,indicator,
                    captureDM,thetaDM,lambdaDM,
                    captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  NLL = ic.neg.log.likelihood(logit.est, Data,indicator,
                              captureDM,thetaDM,lambdaDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET)
   np = length(logit.est)  # number of parameters
  AIC = 2*np + 2*NLL  #  Akaike information criterion (AIC)
    n = sum(Data$counts)  # total captured
 AICc = AIC + 2*np*(np+1)/(n-np-1) # after the correction
 
 return(AICc)
  
} # end of ic.AICc 

###############################################################################