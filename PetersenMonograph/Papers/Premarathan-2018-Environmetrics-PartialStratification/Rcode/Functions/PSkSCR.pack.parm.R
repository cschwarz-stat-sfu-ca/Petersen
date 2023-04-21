############################################################################### 
#####            "PSkSCR.pack.parm" creates initial estimates using      ######
#####              design matrices for the optimization routine for      ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
############################################################################### 

# Function : PSkSCR.pack.parm
#  Input : initial estimates are in regular form(length is  (ncats+1)*(st+1) )
#          and raw data
#  Output : vector in logit/log form(length is  (ncats+1)*(st+1)-1)
#           -1 is there because of the constraint that the category 
#           proportions equals to 1. so only need to estimate ncats-1
#           parameters for category proportions.
          

PSkSCR.pack.parm <- function(init.est.full,data,
                             captureDM,thetaDM,lambdaDM, p_lossDM){
  
  cats <- data$category 
  
  # number of categories 
  ncats <- length(cats)
  
  # number of sampling occasions
  st <-  nchar(data$history[1])

  # capture probabilities: considered equal capture 
  # probability for all the categories in all capture occasions 
  capture.rate <- init.est.full[1:(ncats*st)] 
  
  # category proportions. This should add to 1.
  lambda <- init.est.full[((ncats*st)+1):((ncats*st)+ ncats)] 
  
  # sampling fraction for each samppling time
  theta  <- init.est.full[((ncats*st)+ ncats+1):((ncats*st)+ ncats +st)]  
  
  # loss on capture probabilities at each sampling time
  p_loss <- init.est.full[((ncats*st)+ ncats + st +1):((ncats*st)+ ncats + (2*st))]
  p_loss[p_loss==0] <- 0.000001 # need to give a small values if loss on capture 
                                # probabilities are zero.(for logit conversion)
  
  N      <- init.est.full[length(init.est.full)]
  
  if(ncol(captureDM)==0){
    logit.capture.rate.hat <- rep(NULL, nrow(captureDM))
  } else {
    # use least squares to estimate the capture probabilities based on captureDM
    logit.capture.rate <- logit(capture.rate) 
    # transform the estimates using the logit function
    logit.capture.rate.hat   <- LS(captureDM, logit.capture.rate) 
  }
  
  
  if(ncol(lambdaDM)==0){
    logit.cats.proportion.hat <- rep(NULL, nrow(lambdaDM))
  } else {
    #initial reduced estimates for lambda
    logit.cats.proportion <- lambda.to.cumulative.logit.lambda(lambda) 
    # use least squares to estimate  based on lambdaDM
    logit.cats.proportion.hat <- LS(lambdaDM,logit.cats.proportion) 
  }
  
  
  if(ncol(thetaDM)==0){
    logit.theta.hat <- rep(NULL, nrow(thetaDM))
  } else {   
    # use least squares to estimate the theta based on thetaDM
    logit.theta <- logit(theta)
    # transform the estimates using the logit function
    logit.theta.hat <- LS(thetaDM,logit.theta)
  }
  
  
  if(ncol(p_lossDM)==0){
    logit.p_loss.hat <- rep(NULL, nrow(p_lossDM))
  } else {   
    # use least squares to estimate the theta based on thetaDM
    logit.p_loss <- logit(p_loss)
    # transform the estimates using the logit function
    logit.p_loss.hat <- LS(p_lossDM,logit.p_loss)
  }
  
  
  # transform the initial value for the population size to log scale
  log.N <- log(N)
  
  pack.ini.est <- c(logit.capture.rate.hat,logit.cats.proportion.hat,
                    logit.theta.hat, logit.p_loss.hat, log.N )
  
  return(pack.ini.est)
  
}# end of function "PSkSCR.pack.parm"

###############################################################################
################ Test the function PSkSCR.pack.parm ###########################
# 
# ## Following code is  to check the above function
# 
# test.PSkSCR.pack.parm  = function(){
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
#    
# #     # Set up some values of the parameters for testing
# #     parm <- NULL
# #     parm$p <- matrix(seq(.1,.4,.1),nrow=2)
# #     parm$lambda <- c(0.4, 0.6)
# #     parm$theta <- c(0.3, 0.6)
# #     parm$p_loss <- c(0.05,0.03)
# #     parm$N <- 500
# #     print(parm)
# 
#     indicator <- PSkSCR.create.indicator(data)
#     initial.est = PSkSCR.initial.estimates(data, indicator)
#     pack.init.est = PSkSCR.pack.parm(initial.est$full,data,
#                                      captureDM, thetaDM,lambdaDM, p_lossDM)
#     res  = PSkSCR.unpack.parm(pack.init.est,data,captureDM,thetaDM,lambdaDM,p_lossDM,
#                               captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
#     print(res)
# 
# }
# 
# test.PSkSCR.pack.parm()

###############################################################################
###############################################################################
