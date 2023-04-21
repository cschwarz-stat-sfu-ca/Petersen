###############################################################################
############ unpack parameters  for the data with individual covariate s#######
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
###############################################################################  

# function "ic.unpack.parm " extract the vector of parameters
# breaks the estimates into 3 groups
#  1. Capture rates ( this is a matrix of n.rows x (length(cats)x2) where 
#     n.rows  = number of rows in the data set
#  2. category proportions (this is vector of size length(cats))
#  3. sampling fractions (this is a vector of size 2) ; theta1 and theta 2 for
#     sampling time 1 and 2

# Input : A vector of reduced estimates in logit scale, Data, 
#         Design Matrices and OFFSET vectors
# Output: matrix for capture rates, a vector for category proportions, and
#         a vector for sub-sampleing proportions as stated above 3 groups

ic.unpack.parm = function(est, Data,captureDM,thetaDM,lambdaDM,
                          captureOFFSET,thetaOFFSET,lambdaOFFSET){
  #browser()
  # extract the parameter estimates in full form from the parameter vector
  
  # est is the parameters in logit form.
  # parameters are in the order
  #   beta parameters for capture probabilities, 
  #   beta parameters for category proportion,
  #   beta parameters for sampling fraction
  
  # cats is the vector of categories.
  cats <- Data$category
  
  parm <- NULL
  
  parm$beta.full.red  <- est
  
  # parm$logit.p  is a matrix of n.rows x (2*length(cats))  where n.rows  is 
  # the number of rows in the data set( number of individual histories)
  
  # if statements extract the reduced estimates. extract what is actually optimised
  
  if(ncol(captureDM)== 0){ # If Capture design matrix has zero columns for each individual 
    parm$logit.p  <- captureOFFSET 
  } else {   
    # extract the reduced estimates to the full estimates using appropriate design matrix
    parm$logit.p.red  <- est[1:ncol(captureDM)]
    temp <- captureDM %*% parm$logit.p.red 
    
    parm$logit.p  <- matrix(temp, ncol=(2*length(cats)), byrow=FALSE)+ captureOFFSET 
  }
  
  
  if(ncol(lambdaDM)==0){  # If lambda design matrix has zero columns
    # this is  vector of size "length(cats)-1"
    parm$logit.lambda <- lambdaOFFSET # this is  vector of size "length(cats)-1"
  } else { 
    # extract the reduced estimates to the full
    #        estimates using appropriate design matrix
    
    # this is  vector of size "length(cats)-1"
    parm$logit.lambda.red <- est[(ncol(captureDM)+1):(ncol(captureDM)+ ncol(lambdaDM))]     
    # this is  vector of size "length(cats)-1"
    parm$logit.lambda <- (lambdaDM %*% parm$logit.lambda.red ) + lambdaOFFSET 
  }
  
  
  if(ncol(thetaDM)== 0){  # If theta design matrix has zero columns
    #parm$logit.theta.red  = rep(NULL, nrow(thetaDM))
    parm$logit.theta  <- thetaOFFSET
  } else { 
    # extract the reduced estimates to the full 
    #       estimates using appropriate design matrix
    parm$logit.theta.red  <- est[(ncol(captureDM)+ ncol(lambdaDM)+1): ( ncol(captureDM)+ 
                                                                         ncol(lambdaDM)+
                                                                         ncol(thetaDM) )]
    parm$logit.theta  <- (thetaDM %*% parm$logit.theta.red) + thetaOFFSET
  }
  
  
  parm$beta.full   = c(est[1:ncol(captureDM)], parm$logit.lambda, parm$logit.theta ) 
  
  ####### convert from logit form to regular form #############
  
  parm$p      <- expit(parm$logit.p)
  
  #   cumulative.logit.category.proportion has length(cats)-1 values since
  #      it preserve the constraint "sum to one"
  # category proportions add to 1.
  # parm$lambda is a vector of size length(cats)
  parm$lambda <- as.vector(cumulative.logit.lambda.to.lambda(parm$logit.lambda)) 
  
  parm$theta  <- as.vector(expit(parm$logit.theta))
  
  ##### give appropriate names ######
  colnames(parm$p) <-c(paste(cats,c("_time_1") ,sep=""), paste(cats,c("_time_2") ,sep=""))
  names( parm$lambda)<- c(paste("lambda_",cats,sep=""))
  names( parm$theta) <-  c("theta_1", "theta_2")
  
    return(parm)
  
} # end of ic.unpack.parm 


###############################################################################
###############################################################################


