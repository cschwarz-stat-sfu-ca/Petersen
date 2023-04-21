###############################################################################
############# pack parameters  for the data with individual covariate s########
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
###############################################################################  

#The function "ic.pack.parm" creates initial estimates for the optimization routine.
# Input : initial estimates are in regular form(length is  (3*length(cats)+2)),
#         and raw data
# Output : vector in logit form(length is 
#                   ncol(initial.captureDM) + ncol(thetaDM) + ncol(lambdaDM) )

ic.pack.parm = function(init.est,Data,captureDM,thetaDM,lambdaDM){
  
  cats <- Data$category 
  # capture rates for the sample time 1
  p1     <- init.est[1:length(cats)] 
  # capture rates for the sample time 2
  p2     <- init.est[(length(cats)+1): (2*length(cats))] 
  # category proportions. This should add to 1.
  lambda <- init.est[(2*length(cats)+1):(3*length(cats))] 
  # sampling fraction for sample 1 and sample 2
  theta  <- init.est[(3*length(cats)+1):(3*length(cats)+2)]  
  
  
  min_cap_prob <- min(c(p1,p2))
  max_cap_prob <- max(c(p1,p2))
  
  if(min_cap_prob==max_cap_prob){
    min_cap_prob <- max( 0.001, (min_cap_prob-0.05))
    max_cap_prob <- min( 0.9, (min_cap_prob+0.05))
  }
  
  
  if(ncol(captureDM)==0){
    beta.capture.rate.hat <- rep(NULL, nrow(captureDM))
  } else {
    # use least squares to estimate  based on captureDM 
    logit.capture.rate <- logit(runif(nrow(captureDM),min_cap_prob,max_cap_prob))
    beta.capture.rate.hat <- LS(captureDM, logit.capture.rate)[,1]
  }
 
    
  if(ncol(lambdaDM)==0){
    beta.lambda.hat <- rep(NULL, nrow(lambdaDM))
  } else {
        #initial reduced estimates for lambda
        logit.cats.proportion <- lambda.to.cumulative.logit.lambda(lambda) 
        # use least squares to estimate  based on lambdaDM
        beta.lambda.hat <- LS(lambdaDM,logit.cats.proportion) 
    }

  
  if(ncol(thetaDM)==0){
    beta.theta.hat <- rep(NULL, nrow(thetaDM))
  } else {
       # transform the estimates using the logit function
       logit.theta <- logit(theta)
       # use least squares to estimate the theta based on thetaDM
       beta.theta.hat <- LS(thetaDM,logit.theta)
    }

  return(c(beta.capture.rate.hat, beta.lambda.hat, beta.theta.hat))
}

###############################################################################
###############################################################################
#Following code is  to check the above function

# testing.ic.pack.parm = function(){
#   temp1 = initial.estimates(Data)
#   temp2 = ic.pack.parm(temp1$full,Data,initial.captureDM,thetaDM,lambdaDM)
#   temp  = unpack.parm(temp2,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
#                       thetaOFFSET,lambdaOFFSET)
#   temp
# 
# }
# 
# testing.ic.pack.parm()
#
###############################################################################
###############################################################################
# 
# #### for testing #### delete later..
# Data$History  <- Data$History[1:100]
# Data$length   <- Data$length[1:100]
# Data$lengthsq <- Data$length^2
# Data$counts   <- rep(1, length(Data$length))
# str(Data)
# aaa <- logit(runif(400,0.2,0.4))
# aaa <- logit(rep(.2,400))
# df<-data.frame(cbind(captureDM,aaa))
# bbb<- glm(aaa~length+lengthsq+category+time, data=df)
# ddd<-GI(captureDM, aaa)
# summary(bbb)


# ###############################################################################
##   testing capture prbabilities ####
# 
#    Data <- NULL
#    # including the "00" history
#    Data$History <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U", "MM", "FF", "0M")
#    Data$length <- c(19.73, 23.50, 24.69, 22.76, 20.57, 27.32, 22.24, 18.20, 22.49, 26.59, 21.42, 22.98)
#    Data$lengthsq <- Data$length^2
#    Data$category <- c("M","F")
#    str(Data)
#    print(Data)
#    
#    n.rows <- length(Data$History) # number of rows of data
#    
#    # give the capture formula
#    captureformula <- ~length+ lengthsq+ category+time
#    
#    # number of beta parameters related to capture formula
#    #  need add 1, because of the intercept
#    n.beta.param.cap <- length(colnames(attr(terms(captureformula), "factors"))) +1
#    
#    # design Matrices
#    captureDM <- ic.create.DM(Data,captureformula)
#    thetaDM   = create.DM(c(1,2)) 
#    lambdaDM  = create.DM(c(1)) 
#    
#    #give the offset vectors(vectors of zero's should be given since no restriction)
#    captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(Data$category) *2) )
#    thetaOFFSET   = c(0,0)
#    lambdaOFFSET  = c(0)
# 
#    # Set up some values of the parameters for testing
#    parm <- NULL
#    parm$p <- matrix(c(0.1,0.2,0.3,0.4,
#                       0.2,0.3,0.1,0.3,
#                       0.1,0.4,0.3,0.2,
#                       0.3,0.3,0.2,0.2,
#                       0.4,0.2,0.3,0.1,
#                       0.2,0.1,0.1,0.2,
#                       0.2,0.2,0.3,0.1,
#                       0.4,0.2,0.3,0.1,
#                       0.2,0.3,0.3,0.2,
#                       0.3,0.3,0.2,0.2,
#                       0.3,0.1,0.1,0.2,
#                       0.1,0.1,0.1,0.2),ncol=4, byrow=TRUE)
#    parm$lambda <- c(0.4, 0.6)
#    parm$theta <- c(0.8, 0.6)
#    print(parm)
#    
#    indicator <- ic.create.indicator(Data)
#    
#    #### pack capture probabilities #####
#    if(ncol(captureDM)==0){
#      test.beta.cap.rate <- rep(NULL, nrow(captureDM))
#    } else {
#      # use generalized linear models to estimate the  beta parameters of capture 
#      # probabilities based on initial.captureDM
#      test.logit.cap.rate <- logit(as.vector(parm$p))
#      test.logit.beta.cap.rate <- LS(captureDM, test.logit.cap.rate)[,1]
#    }
#    
#    df<- data.frame(captureDM)
#    # give the capture formula
#   cap.fom<- test.logit.cap.rate~length+ lengthsq+ category+time
#   temp<- lm(as.formula(cap.fom), df)
#    summary(temp)$coefficients[,1:2]
# 
#    #### unpack capture probabilities ####
#    cats <- Data$category
#    if(ncol(captureDM)== 0){ # If Capture design matrix has zero columns for each individual 
#      parm$test.logit.p  <- captureOFFSET 
#    } else {   
#      # extract the reduced estimates to the full estimates using appropriate design matrix
#      temp <- captureDM %*% test.logit.beta.cap.rate
#      
#      cap.rate.logit <- matrix(temp, ncol=(2*length(cats)), byrow=FALSE)+ captureOFFSET 
#    }
# 
# new.test.logit.cap.rate<- as.vector(cap.rate.logit)
# new.test.logit.beta.cap.rate <- LS(captureDM, new.test.logit.cap.rate)[,1]
# 
# 
# 
# plot(Data$length ,cap.rate.logit[,4])
#    
# plot(Data$length ,(parm$p- expit(cap.rate.logit))[,1])
# plot(Data$length ,(parm$p- expit(cap.rate.logit))[,1])
   
##########################################################################################