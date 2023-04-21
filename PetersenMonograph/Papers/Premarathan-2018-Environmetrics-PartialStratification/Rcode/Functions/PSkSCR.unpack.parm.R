############################################################################### 
##### Function "PSkSCR.unpack.parm " extract the vector of parameters for #####
#####  Partial Stratification in k-Sample Capture-Recapture Experiments   #####
#####                     with known dead removals                        #####
###############################################################################

# function "PSkSCR.unpack.parm " extract the vector of parameters
# breaks the estimates into 4 groups
#  1. Capture rates ( this is a matrix of ( length(cats)x st )  ) where 
#     length(cats) = no. of categories,  st = number of sampling occasions
#  2. category proportions (lambda) (this is vector of size length(cats))
#  3. sampling fractions (this is a vector of size st) ; 
#     theta1, theta 2, ... for sampling time 1, sampling time  2, ...
#  4. population size, N

# Input : A vector of reduced estimates in logit/log scale, data, 
#         Design Matrices and OFFSET vectors
# Output: A vector of probabilities of length (ncats+1)*(st+1) 
#         the order of the probabilities is as stated above 4 groups

PSkSCR.unpack.parm <- function(est, data,
                              captureDM,thetaDM,lambdaDM, p_lossDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET, p_lossOFFSET){
  
  # extract the parameter estimates in full form from the parameter vector
  
  # parameters are in the order
  #   logit.capture.rates, cumulativ.logit.category.proportion,
  #   logit.sampling.fraction, log.population.size.
  # cumulative.logit.category.proportion has length(cats)-1 values since
  #      it preserve the constraint "sum to one"
  # category proportions add to 1.
  # cats is the is the vector of categories.
  # est are in logit  or log form.
  
  cats = data$category
  # number of sampling occasions
  st <-  nchar(data$history[1])
  
  parm <- NULL
  parm$logit.full.red <- est
  
    
  if(ncol(captureDM)== 0){ # If Capture design matrix has zero columns  
    parm$logit.p  <- captureOFFSET 
  } else {   
    # extract the reduced estimates to the full estimates using appropriate
    # design matrix
    parm$logit.p.red  <- est[1:ncol(captureDM)]
    parm$logit.p      <- (captureDM %*% parm$logit.p.red) + captureOFFSET 
  }
  
  
  if(ncol(lambdaDM)==0){  # If lambda design matrix has zero columns
    # this is  vector of size "length(cats)-1"
    parm$logit.lambda <- lambdaOFFSET 
  } else { 
    # extract the reduced estimates to the full
    #        estimates using appropriate design matrix
    
    # this is  vector of size "length(cats)-1"
    parm$logit.lambda.red <- est[(ncol(captureDM)+1):(ncol(captureDM)+ 
                                                       ncol(lambdaDM))]   
    # this is  vector of size "length(cats)-1"
    parm$logit.lambda <- (lambdaDM %*% parm$logit.lambda.red ) + lambdaOFFSET 
  }
  
  
  if(ncol(thetaDM)== 0){  # If theta design matrix has zero columns
    parm$logit.theta  <- thetaOFFSET
  } else { 
    # extract the reduced estimates to the full 
    #       estimates using appropriate design matrix
    parm$logit.theta.red  <- est[( ncol(captureDM)+ ncol(lambdaDM)+1 ): 
                                  ( ncol(captureDM)+ ncol(lambdaDM)+
                                      ncol(thetaDM) ) ]
    parm$logit.theta  <- (thetaDM %*% parm$logit.theta.red) + thetaOFFSET
  }
  
  
  if(ncol(p_lossDM)== 0){  # If theta design matrix has zero columns
    parm$logit.p_loss  <- p_lossOFFSET
  } else { 
    # extract the reduced estimates to the full 
    #       estimates using appropriate design matrix
    parm$logit.p_loss.red  <- est[( ncol(captureDM)+ ncol(lambdaDM)+ 
                                      ncol(thetaDM) + 1 ): 
                                   ( ncol(captureDM)+ ncol(lambdaDM)+ 
                                       ncol(thetaDM)+ ncol(p_lossDM) ) ]
                                       
    parm$logit.p_loss <- (p_lossDM %*% parm$logit.p_loss.red) + p_lossOFFSET
  }
  
  
  parm$log.N  <- est[length(est)]
  parm$logit.full   <- c(parm$logit.p, parm$logit.lambda,
                        parm$logit.theta, parm$logit.p_loss, parm$log.N)   
  
  
  # convert from logit or log form to regular form.
  capture.rates <- as.vector(expit(parm$logit.p)) # in vector form
  names(capture.rates) <- c(paste(rep(cats, st),"_Time_",
                                  rep(c(1:st),each=length(cats)), sep=""))
  
  parm$p      <- matrix(capture.rates, nrow = length(cats),
                       ncol = st,byrow=FALSE) # rows are categories
  
  # parm$lambda is a vector of size length(cats)
  parm$lambda <- as.vector(cumulative.logit.lambda.to.lambda(parm$logit.lambda)) 
  parm$theta  <- as.vector(expit(parm$logit.theta))
  parm$p_loss <- as.vector(expit(parm$logit.p_loss))
  parm$N      <- exp(parm$log.N)
  
  # estimate expected number of fish in each category
  parm$N_lambda <- parm$N * parm$lambda
  
  colnames(parm$p) <- c(paste("Time_",c(1:st), sep=""))
  rownames(parm$p) <- c(paste(cats))
  
  names( parm$lambda) <- c(paste("lambda_",cats,sep=""))
  names( parm$theta) <-  c(paste("theta_",c(1:st), sep=""))
  names(parm$p_loss) <-  c(paste("loss_Time_",c(1:st), sep=""))
  names( parm$N) <-  c("N")
  
  #  parm$full is a vector of length 3*length(cats)+3
  parm$full   <- c(capture.rates, parm$lambda, parm$theta, parm$p_loss, parm$N) 
  
  return(parm)
  
} # end of PSkSCR.unpack.parm

###############################################################################
###############################################################################

###############################################################################
################ Test the functionPSkSCR.unpack.parm ##########################
# 
# #Following code is  to check the above function
# 
# test.PSkSCR.unpack.parm = function(){
# 
#     data <- NULL
#     data$history <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U", "M0","FF","0F")
#     data$counts   <- c(41,9,16,4,23,7,12,25,63,-2,-1,-3)
#     data$category <- c("M","F")
#     str(data)
# 
#     # give  the required design matrices for capture recapture probabilities,
#     # theta(sampling(sexing) fractions ) and lambda(Category proportion)
#     captureDM = create.DM(c(1,2,3,4)) 
#     thetaDM   = create.DM(c(1,2)) 
#     lambdaDM  = create.DM(c(1)) 
#     p_lossDM  = create.DM(c(1,2)) 
#     
#     #give the offset vectors(vectors of zero's should be given since no restriction)
#     captureOFFSET = c(0,0,0,0) 
#     thetaOFFSET   = c(0,0)
#     lambdaOFFSET  = c(0)
#     p_lossOFFSET  = c(0,0)
#    
#     # Set up some values of the parameters for testing
# #     parm <- NULL
# #     parm$p <- matrix(seq(.1,.4,.1),nrow=2)
# #     parm$lambda <- c(0.4, 0.6)
# #     parm$theta <- c(0.3, 0.6)
# #     parm$p_loss <- c(0.05,0.03)
# #     parm$N <- 500
# #     print(parm)
# 
# 
#     indicator <- PSkSCR.create.indicator(data)
#     initial.est <- PSkSCR.initial.estimates(data, indicator)
#     initial.pack.parm.est <- PSkSCR.pack.parm(initial.est$full,data,
#                                               captureDM,thetaDM,lambdaDM, p_lossDM)
# 
#     res <- PSkSCR.unpack.parm(initial.pack.parm.est, data,
#                               captureDM,thetaDM,lambdaDM,p_lossDM,
#                               captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
#     print(res)
# 
# }
#  
# test.PSkSCR.unpack.parm()

################## END of test  function  #####################################
###############################################################################



