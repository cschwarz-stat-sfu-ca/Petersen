###############################################################################
############################################################################### 
# function "AICc" calculates the AIC with a correction for finite sample sizes
#
# AIC = 2*K - 2*log(likelihood) where k is the number of parameters in the model
#    
# AICc = AIC + 2*K*(k+1)/(n-k-1)  ; n denotes the sample size
#
# Input : vector of estimate in logit/log form and Data
#      
# Output: AICc Value. 

AICc <- function(logit.est, Data,captureDM,thetaDM,lambdaDM,
                 captureOFFSET,thetaOFFSET,lambdaOFFSET){
  NLL = neg.log.likelihood(logit.est, Data,captureDM,thetaDM,lambdaDM,
                           captureOFFSET,thetaOFFSET,lambdaOFFSET)
   np = length(logit.est)  # number of parameters
  AIC = 2*np + 2*NLL  #  Akaike information criterion (AIC)
    n = sum(Data$counts)  # total captured
 AICc = AIC + 2*np*(np+1)/(n-np-1) # after the correction
 
 return(AICc)
  
} 

###############################################################################
# testing the above function

# test.AICc = function(Data){
#   temp1 = initial.estimates(Data)
#   temp2 = pack.parm(temp1,Data,captureDM,thetaDM,lambdaDM)
#   temp3 = AICc(temp2,Data,captureDM,thetaDM,lambdaDM,
#                captureOFFSET,thetaOFFSET,lambdaOFFSET)
#   temp3
# }
# 
# test.AICc(Data)
###############################################################################
