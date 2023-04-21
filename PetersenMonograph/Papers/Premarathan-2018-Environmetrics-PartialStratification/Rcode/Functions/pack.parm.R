###############################################################################
###############################################################################

#The function "pack.parm" creates initial estimates for the optimization routine.
# Input : initial estimates are in regular form(length is  (3*length(cats)+3)),
#         and raw data
# Output : vector in logit/log form(length is  (3*length(cats)+3)-1)

pack.parm = function(init.est,Data,captureDM,thetaDM,lambdaDM){
  
  cats = Data$category 
  # capture rates for the sample 1
  p1     = init.est[1:length(cats)] 
  # capture rates for the sample 2
  p2     = init.est[(length(cats)+1): (2*length(cats))] 
  # category proportions. This should add to 1.
  lambda = init.est[(2*length(cats)+1):(3*length(cats))] 
  # sampling fraction for sample 1 and sample 2
  theta  = init.est[(3*length(cats)+1):(3*length(cats)+2)]  
  N      = init.est[length(init.est)]
  
  if(ncol(captureDM)==0){
    logit.capture.rate.hat = rep(NULL, nrow(captureDM))
  } else {
      # use least squares to estimate the capture probabilities based on captureDM
    capture.rate = c(p1,p2)
    logit.capture.rate = logit(capture.rate) 
    # transform the estimates using the logit function
    logit.capture.rate.hat   = LS(captureDM, logit.capture.rate) 
    }
  
 
  if(ncol(lambdaDM)==0){
    logit.cats.proprtion.hat = rep(NULL, nrow(lambdaDM))
  } else {
        #initial reduced estimates for lambda
        logit.cats.proprtion = lambda.to.cumulative.logit.lambda(lambda) 
        # use least squares to estimate  based on lambdaDM
        logit.cats.proprtion.hat = LS(lambdaDM,logit.cats.proprtion) 
    }


  if(ncol(thetaDM)==0){
    logit.theta.hat = rep(NULL, nrow(thetaDM))
  } else {   
    # use least squares to estimate the theta based on thetaDM
    logit.theta = logit(theta)
    # transform the estimates using the logit function
    logit.theta.hat = LS(thetaDM,logit.theta)
    }

  
  # transform the initial value for the population size to log scale
  log.N = log(N)

  return(c(logit.capture.rate.hat,logit.cats.proprtion.hat,logit.theta.hat,log.N ))
}

###############################################################################
###############################################################################
#Following code is  to check the above function

# testing.pack.parm = function(){
#   temp1 = initial.estimates(Data)
#   temp2 = pack.parm(temp1,Data,captureDM,thetaDM,lambdaDM)
#   temp  = unpack.parm(temp2,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
#                       thetaOFFSET,lambdaOFFSET)
#   temp
# 
# }
# 
# testing.pack.parm()
#
###############################################################################
###############################################################################
