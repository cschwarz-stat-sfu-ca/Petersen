###############################################################################
#####         AIC with a correction for finite sample sizes for          ######
#####  Partial Stratification in k-Sample Capture-Recapture Experiments  ###### 
#####                     with known dead removals                       ###### 
###############################################################################


############################################################################### 
# function "PSkSCR.AICc" calculates the AIC with a correction for finite sample sizes
#
# AIC = 2*K - 2*log(likelihood) where k is the number of parameters in the model
#    
# AICc = AIC + 2*K*(k+1)/(n-k-1)  ; n denotes the sample size
#
# Input : vector of estimate in logit/log form and Data
#      
# Output: AICc Value. 

PSkSCR.AICc <- function(logit.est, data,indicator,
                        captureDM,thetaDM,lambdaDM,p_lossDM,
                        captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET){
  
  NLL <- PSkSCR.neg.log.likelihood(logit.est, data,indicator,
                           captureDM,thetaDM,lambdaDM,p_lossDM,
                           captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
  np <- length(logit.est)  # number of parameters
  AIC <- 2*np + 2*NLL  #  Akaike information criterion (AIC)
  n <- sum(abs(data$counts))  # total captured
  AICc <- AIC + 2*np*(np+1)/(n-np-1) # after the correction
  
  return(AICc)
  
} 

###############################################################################
# # testing the above function
# 
# PSkSCR.test.AICc = function(data){
#   temp1 = PSkSCR.initial.estimates(data)$full
#   indicator <- PSkSCR.create.indicator(data)
#   temp2 = PSkSCR.pack.parm(temp1,data,captureDM,thetaDM,lambdaDM)
#   temp3 = PSkSCR.AICc(temp2,data,indicator,
#                       captureDM,thetaDM,lambdaDM,
#                       captureOFFSET,thetaOFFSET,lambdaOFFSET,p_loss)
#   temp3
# }
# 
# PSkSCR.test.AICc(data)
###############################################################################
