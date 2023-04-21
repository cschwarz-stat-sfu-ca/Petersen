###############################################################################
#####                         Bayesian DIC for                           ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################
# function "PSkSCR.BAYESIAN.DIC" calculates the deviance information criterion (DIC)
# for finite sample sizes where the posterior distribution of the model has 
# been obtained by Markov Chain Monte Carlo(MCMC) simulations. 
# 
#### Two different ways of DIC calculation ########
#
#  DEVIANCE=  D(theta)= -2 * log(likelihood) 
#
#   
#        "Dbar" is the posterior mean of the deviance,
#               (i.e. the average of D(theta) over the samples theta)
#        "D(theta_bar)" is the value of the Deviance evaluated at the average of
#                       the samples of theta
#        "pD" is  the effective number of parameters
#         pD= "posterior mean deviance - deviance of posterior means"
#
### Method 1. ####  
#       pD = Dbar - D(theta_bar)     # Spiegelhalter et al. (2002)
#           
### Method 2. ####               
#       pD = pv = 1/2  *  var(D(theta))   # Gelman et al (2004)
#
#  DIC = DIC = Dbar + pD   or  DIC = D(thata_bar) + 2 pD
#
# Function : bayesian.DIC 
#  Input : matrix of posterior values of parameters in logit/log scale (after
#          burn in and after thinning), data, idicator matrices, design matrices 
#          and offset vectors
#  Output : DIC value calculated in two standard ways
#
#

PSkSCR.BAYESIAN.DIC <- function(thinned.beta.matrix , data,indicator,
                                captureDM,thetaDM,lambdaDM,p_lossDM,
                                captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET){
  
  result <- NULL
  
  ####### method 1 ######
  # negative log likelihood for each posterior (thinned) sample
  NLL <- adply(thinned.beta.matrix,1,PSkSCR.neg.log.likelihood,
               data,indicator,
               captureDM,thetaDM,lambdaDM,p_lossDM,
               captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
  
     D <- 2 * NLL[,2] # Deviance for each posterior (thinned) sample
  result$Dbar <- mean(D)     # posterior mean of deviance
      
  # average of posterior samples (logit/log scale)
  avg <- adply(thinned.beta.matrix,2,mean)
  
  # value of the Deviance evaluated at the average of posterior samples
  D_theta_bar <- 2* PSkSCR.neg.log.likelihood(avg[,2], data,indicator,
                                              captureDM,thetaDM,lambdaDM,p_lossDM,
                                              captureOFFSET,thetaOFFSET,
                                              lambdaOFFSET,p_lossOFFSET)
    
  result$pD  <- result$Dbar - D_theta_bar # Spiegelhalter et al. (2002)
  
  result$DIC$pD <- result$Dbar + result$pD 
  
  ####### method 2 ######
  result$pv <- var(D)/2  # Gelman et al (2004)
  
  result$DIC$pv <- result$Dbar + result$pv
  
  
  return(result)
}


###############################################################################