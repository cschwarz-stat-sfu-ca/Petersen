###############################################################################
###############################################################################
# function "unpack.parm " extract the vector of parameters
# breaks the estimates into 4 groups
#  1. Capture rates ( this is a length(cats)x2 matrix) where 
#     length(cats) = no. of categories. 2 columns are sampling time 1 and 2
#  2. category proportions (this is vector of size length(cat))
#  3. sampling fractions (this is a vector of size 2) ; theta1 and theta 2 for
#     sampling time 1 and 2
#  4. population size, N

# Input : A vector of reduced estimates in logit/log scale, Data, 
#         Design Matrices and OFFSET vectors
# Output: A vector of probabilities of length 3*length(cats)+3 ; in order as stated above 4 groups

unpack.parm = function(est, Data,captureDM,thetaDM,lambdaDM,
                       captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  # extract the parameter estimates in full form from the parameter vector
  
  # parameters are in the order
  #   logit.capture.rates, cumulativ.logit.category.proportion,
  #   logit.sampling.fraction, log.population.size.
  # cumulative.logit.category.proportion has length(cats)-1 values since
  #      it preserve the constraint "sum to one"
  # category proportions add to 1.
  # cats is the is the vector of categories.
  # est are in logit  or log form.
  
  cats = Data$category
  
  parm = NULL
  parm$logit.full.red   = est
  
  # if statements extract the reduced estimates. extract what is actually optimised
  
  if(ncol(captureDM)== 0){ # If Capture design matrix has zero columns  
    #parm$logit.p.red  = rep(NULL, nrow(captureDM))
    parm$logit.p      = captureOFFSET 
  } else {   
    # extract the reduced estimates to the full estimates using appropriate design matrix
      parm$logit.p.red  = est[1:ncol(captureDM)]
      parm$logit.p      = (captureDM %*% parm$logit.p.red) + captureOFFSET 
    }
  
  
  if(ncol(lambdaDM)==0){  # If lambda design matrix has zero columns
    #parm$logit.lambda.red = rep(NULL, nrow(lambdaDM))  
    # this is  vector of size "length(cats)-1"
    parm$logit.lambda = lambdaOFFSET # this is  vector of size "length(cats)-1"
  } else { 
    # extract the reduced estimates to the full
    #        estimates using appropriate design matrix
    
    # this is  vector of size "length(cats)-1"
      parm$logit.lambda.red = est[(ncol(captureDM)+1):(ncol(captureDM)+ ncol(lambdaDM))]   
      # this is  vector of size "length(cats)-1"
      parm$logit.lambda = (lambdaDM %*% parm$logit.lambda.red ) + lambdaOFFSET 
    }
  
      
  if(ncol(thetaDM)== 0){  # If theta design matrix has zero columns
    #parm$logit.theta.red  = rep(NULL, nrow(thetaDM))
    parm$logit.theta  = thetaOFFSET
  } else { 
    # extract the reduced estimates to the full 
    #       estimates using appropriate design matrix
      parm$logit.theta.red  = est[(ncol(captureDM)+ ncol(lambdaDM)+1): ( ncol(captureDM)+ 
                                                                           ncol(lambdaDM)+
                                                                           ncol(thetaDM) )]
      parm$logit.theta  = (thetaDM %*% parm$logit.theta.red) + thetaOFFSET
    }
  
  
  parm$log.N  = est[length(est)]
  parm$logit.full   = c(parm$logit.p, parm$logit.lambda,
                        parm$logit.theta, parm$log.N)   
  
  
  # convert from logit or log form to regular form.
  capture.rates = as.vector(expit(parm$logit.p)) # in vector form
  
  parm$p      = matrix(capture.rates, nrow = length(cats),
                       ncol = 2,byrow=FALSE) # rows are categories
 
  # parm$lambda is a vector of size length(cats)
  parm$lambda = as.vector(cumulative.logit.lambda.to.lambda(parm$logit.lambda)) 
  parm$theta  = as.vector(expit(parm$logit.theta))
  parm$N      = exp(parm$log.N)
  #  parm$full is a vector of length 3*length(cats)+3
  parm$full   = c(capture.rates, parm$lambda, parm$theta, parm$N) 
  
  # estimate expected number of fish in each category
  parm$N_lambda = parm$N * parm$lambda
  
  colnames(parm$p) = c("time1", "time2")
  rownames(parm$p) = c(paste(cats))
  
  names( parm$lambda)= c(paste("lambda_",cats))
  names( parm$theta) =  c("theta_1", "theta_2")
  names( parm$N) =  c("N")
  
  return(parm)
  
} # end of unpack.parm 

###############################################################################
###############################################################################
###############################################################################
#Following code is  to check the above function

# testing.unpack.parm = function(){
#   
#   data = NULL
#   data$History  = c("U0","UU","M0","MM","F0","FF","0M","0F","0U")
#   data$counts   = c(41,9,16,4,23,7,12,25,63)
#   data$category = c("M","F") # for 2 categories
#   #data$category = c("A","B","C")# for 3 categories
#   str(data)
#   #print(data)
#   
#   # Set up some values of the parameters for testing
# 
#   # testing for 2 categories
#   p      = seq(.1,.4,.1)
#   lambda = c(0.4, 0.6)
#   theta  = c(0.3, 0.6)
#   N      = 500
#   
#   # testing for 3 categories
#   p3 = c(0.2,0.1,0.3,0.2,0.4,0.1)
#   l  = c(0.2,0.5,0.3)  
#   
#   # testparm2 is for 2 categories, and restparm3 is for 3 categories. 
#   #      (probabilities in logit and log form)
#   testparm2 = c(log(p/(1-p)),log((lambda[1]/sum(lambda))/(1-(lambda[1]/sum(lambda)))),
#                    log(theta/(1-theta)),log(N))
#   #testparm3 = c(log(p3/(1-p3)),log((l[1]/sum(l))/(1-(l[1]/sum(l)))),
#                      log((l[2]/sum(l[2:3]))/(1-(l[2]/sum(l[2:3])))), 
#                      log(theta/(1-theta)),log(N))
#   
#   capDM = create.DM(c(1,1,3,4))
#   thetaDM =create.DM(c(1,2))
#   #lambdaDM = create.DM(c(1))
#   lambdaDM  = matrix(, ncol=0,nrow=1)  # Design matrix for lambda
# 
#   captureOFFSET = c(0,0,0,0) 
#   thetaOFFSET   = c(0,0)
#   #lambdaOFFSET  = c(0)
#   lambdaOFFSET  = c(logit(.4))
#  
#   res2 = unpack.parm(testparm2,data, capDM,thetaDM,lambdaDM,
#                      captureOFFSET,thetaOFFSET,lambdaOFFSET)  # for 2 categories
#   #res3 = unpack.parm(testparm3,data) # for 3 categories
#   res2 # for 2 categories
#   #res3 # for 3 categories
# }
# 
# 
# testing.unpack.parm()
#
###############################################################################
###############################################################################
