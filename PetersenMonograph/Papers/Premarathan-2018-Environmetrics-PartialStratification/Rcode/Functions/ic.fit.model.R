############################################################################### 
################# fit model  for the data with individual covariates ##########
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
###############################################################################   

# Function "ic.fit.model" 
# Inputs: model id, row data, design matrices and offset vectors for capture 
          # probabilities, category proportions and sub-sample proportions
#  Outputs: MLE , variance, SE of all the parameters
          # ( beta parameters for capture probabilities, 
          # category proportions, sub-sample proportions, and population size),
          # and model information (id, number of parameters,
                                    # negative log-likelihood and AICc values)
          # Plot of Capture Probabilities (logit scale and regular scale)


ic.fit.model <- function(model.id ,Data, true.length,
                         captureformula, thetaDM,lambdaDM, 
                         captureOFFSET,thetaOFFSET,lambdaOFFSET, debug=FALSE){
  
  Result <- NULL
  Result$true_length <- true.length
  
  ncats <- length(Data$category) # number of categories
  Result$cats <- Data$category
  
   
  # Initial Estimates (in regular form) are in following order
  # capture probabilities( M_time1, F_time1, M_time2, F_time2), 
  # lambda( lambda_1, lambda_2) and theta( theta_1, theta_2)
  ic.initial.est <- ic.initial.estimates(Data)
  Result$initial.est <- ic.initial.est
  

  # Design matrix for capture probabilities for all individuals
  ## This is a big matrix with 
  ## number of rows = (length(Data$History) x ncats x 2)
  ## #                                    # 2 for time 1 and time 2
  ## number of columns = (1 + number of individual covariates + 2)
  ## #                 +1 for Intercept and and  2 represent two columns for sex and time
  #
  ## if there are two categories (M and F), then this is a big matrix with 4 blocks
  ## such that each block has length(Data$History) number of rows in the order of 
  ## M_time1, F_time1, M_time2 and F_time2
  
  captureDM <- ic.create.DM(Data,captureformula)
  Result$captureDM <- captureDM
  
  init.est <- ic.initial.est$full
  ic.initial.pack.parm.est <- ic.pack.parm(init.est,Data,
                                           captureDM,thetaDM,lambdaDM) 
  # initial beta estimates
  Result$initial.beta.est <- ic.initial.pack.parm.est
  

  # lower bound for optimization routine
  l.b <-  rep(-10,length(ic.initial.pack.parm.est))
  
  # upper bound for optimization routine
  u.b <-  rep(10,length(ic.initial.pack.parm.est))
 
  # create indicator variables. 
  indicator <- ic.create.indicator(Data)
  
  # optimize the negative log likelihood function(MLEs in logit/log form)
  res <- optim(par=ic.initial.pack.parm.est, fn=ic.neg.log.likelihood, method="L-BFGS-B",
               lower=l.b, upper=u.b , control=list(trace=3,maxit=300), 
               Data=Data, indicator=indicator, 
               captureDM =captureDM, thetaDM=thetaDM,lambdaDM=lambdaDM,
               captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
               lambdaOFFSET=lambdaOFFSET)
  
  Result$res <- res

  #extract the MLEs for the parameter (estimates in regular form)
  Result$est<- ic.unpack.parm(Result$res$par,Data,
                                captureDM,thetaDM,lambdaDM,
                                captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  
  ##############
  # Variance Covariance Matrix of logit/log parameters
  # this is the invers of the hessian matrix
  hess <- ic.inv.fisher.info(Result$res$par,Data, indicator,
                              captureDM, thetaDM,lambdaDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET) 

  
############  Delete later 
#   if(debug) browser()
#   # check the hess above with optimHess() to see if we get the same value
#   hess2 <- optimHess(par=Result$res$par, fn=ic.neg.log.likelihood, 
#                Data=Data, indicator=indicator, 
#                captureDM =captureDM, thetaDM=thetaDM,lambdaDM=lambdaDM,
#                captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
#                lambdaOFFSET=lambdaOFFSET)
############  Delete later 
  
  # reduced Variance Covariance Matrix of logit parameters
  Result$vcv$logit.full.red <- hess   
  #Result$vcv$logit.full.red2<- hess2  # from the optimHess
  
  ###########################
  # MLEs of the beta parameters for capture probabilities
  Result$est$logit.p <- Result$est$logit.p.red 
  
  # variance Covariance Matrix of logit parameters of capture probabilities
  # ( i.e.: beta parameters according to the capture formula)
  Result$vcv$logit.p <- hess[(1:ncol(captureDM)),(1:ncol(captureDM)) ]
  Result$se$logit.p <- sqrt(diag(Result$vcv$logit.p)) 
  
  # Variance Covariance Matrix of logit parameters of lambda and theta 
  Result$vcv$logit.lambda.theta.red <- hess[((ncol(captureDM)+1): length(Result$res$par)),
                                            ((ncol(captureDM)+1): length(Result$res$par))]  
  
  ###############################################
  #find the variance covariance (vcv) of lambda and theta in regular form
  
  #Go from reduced estimates to expanded estimates
  X <- bdiag(lambdaDM,      # lamda parameters in logit form
            thetaDM)       # theta parameters in logit form
  vcv.logit.lambda.theta <- X %*% Result$vcv$logit.lambda.theta.red  %*% t(X)
  Result$vcv$logit.lambda.theta <- vcv.logit.lambda.theta
  Result$se$logit.lambda.theta  <- sqrt(diag(vcv.logit.lambda.theta))
  
  # Go from logit(continuation ratio)  to continuation ratio and then continuation
  # ratio to original p 
  # We use the deltamethod() function in the msm package
  # and so need to create expressions to do the delta method
  # First take anti-logits of all parameters 
  antilogit <- llply(paste("~1/(1+exp(-x",
                          1:(nrow(Result$vcv$logit.lambda.theta)),
                          "))",sep=""),
                    as.formula)
  vcv.cont <- deltamethod(antilogit,
                         Result$est$beta.full[((ncol(captureDM)+1):
                                                            length(Result$est$beta.full))],
                         Result$vcv$logit.lambda.theta,
                         ses=FALSE)
  
  Result$est$logit.lambda.theta.antilogit <-
          expit(Result$est$beta.full[((ncol(captureDM)+1):length(Result$est$beta.full))])
  
  
  anti.cont.formula <- ic.make.anticont.exp(1, ncats) 
  anti.cont.formula <- unlist(anti.cont.formula)
  anti.cont.formula <- c(anti.cont.formula,paste("~x",seq((ncats),(ncats+1)),sep="") )
  anti.cont.formula <- llply(anti.cont.formula,as.formula)
  vcv.lambda.theta <- deltamethod(anti.cont.formula,
                                  Result$est$logit.lambda.theta.antilogit,
                                  vcv.cont, ses=FALSE)
  
  ###############################################
  # extract the individual VCV's  and SE's of the parameters lambda and theta
  
  Result$vcv$lambda <- vcv.lambda.theta[(1:ncats),(1:ncats)]
  Result$se$lambda <- sqrt(diag(Result$vcv$lambda))
  
  Result$vcv$theta <- vcv.lambda.theta[(ncats+1):(ncats+2),(ncats+1):(ncats+2)]
  Result$se$theta <- sqrt(diag(Result$vcv$theta))
  
  ###############################################
  
  #############################################################################
  ## calculate estimate for N ( refer the papers, Huggins 1981 and Huggings 1991,
  ##                            " on the Statistical Analysis of Capture Experiments"
  ##                            "some practical aspect of Conditional likelihood approach...."
  ##                            and Alaho 1990,logistic regression in capture-recapture models)
  phi <- ic.log.prob.history(Result$est, Data$History, Data$category, indicator, phi.flag="YES")$phi
  Result$phi <- phi
  N <- sum(1/phi)
  Result$est$N <- N  # estimate for N 
  
  # unbiased estimator of Var(N)  is S^2
  # note that this under estimate the standard error of N
  Result$VCV$s_sq  <-  sum( (1-phi)/(phi^2)) 
  
  # first derivative  of the log likeligood
  first.derive <- grad(ic.neg.log.likelihood, Result$res$par, method="Richardson", 
                       Data=Data, indicator=indicator,
                       captureDM=captureDM, thetaDM=thetaDM, lambdaDM=lambdaDM,
                       captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
                       lambdaOFFSET=lambdaOFFSET)
  Result$first.derive <- first.derive
  
  D <- sum(1/(phi^2))* first.derive 

# ############## Delete later....   
#   # logit to regular scale - inverse if the hesian matrix
#   antilogit_formula <- llply(c(paste("~1/(1+exp(-x",1:(nrow(hess)),"))",
#                                      sep="")),as.formula)
#   
#   vcv.regular <- deltamethod(antilogit_formula, Result$res$par,
#                          hess, ses=FALSE)
#
#   D <- expit(first.derive)
#
#   Result$VCV$Dt_Iinv_D <- t(D)  %*% vcv.regular %*% (D)
# ##############  
  Result$VCV$Dt_Iinv_D <- t(D)  %*% hess %*% (D) # check later
  

  # Asymptotic variance of N is s^2 +  t(D)* inv(Hessian matrix) *D
  Result$VCV$N <-  Result$VCV$s_sq  + Result$VCV$Dt_Iinv_D
  Result$se$N  <- sqrt( Result$VCV$N)
  
  
  

  ################## calculate the adjusted estimate for N  ###################
  # if the capture probability(phi) for some individual are very very small then
  # the estimate of N is affected ny those values and N will be very big value
  # therefor just consider the phi values larger than 0.0005
  
  adjust.N <- sum(1/phi[phi>0.0005])
  Result$est$adjust.N <- adjust.N  # estimate for adjust.N 
  
  #############################################################################

  # Extract the variance and the SE's of the expected number of individuals
  # in each category.
  
  #estimated category proportions(lambda) and estimated population size(N)
  Result$est$N_lambda <-  Result$est$N *  Result$est$lambda 
  
 
 
 
  ########## Model information #########################
  
  # model identification
  Result$model.id <- model.id
  
  # number of parameters
  Result$np <- length(Result$res$par) 
  
  #Negative log likelihood value
  Result$NLL <- Result$res$value
  
  #corrected AIC information for the model
  Result$AICc <- ic.AICc(Result$res$par,Data,indicator, 
                         captureDM,thetaDM,lambdaDM,
                         captureOFFSET,thetaOFFSET,lambdaOFFSET)
 
  # design matrices and offset vectors for category proportions and 
  # sub-sample proportions

  Result$thetaDM <- thetaDM 
  Result$lambdaDM <- lambdaDM

  Result$thetaOFFSET <- thetaOFFSET
  Result$lambdaOFFSET <- lambdaOFFSET
  #######################################################
  
  ## predictive Plots of Capture Probabilities (logit scale and regular scale) 
  prediction.plots <- ic.plots(Result$res$par,Data, Result$true_length,
                               captureformula, thetaDM,lambdaDM,
                               thetaOFFSET,lambdaOFFSET)
  
  ## Plot of Capture Probabilities (logit scale) with standardized length
  Result$pred.logit.ggplot.p <- prediction.plots$pred.logit.ggplot
  
  ## Plot of Capture Probabilities (regular scale) with standardized length
  Result$pred.ggplot.p <- prediction.plots$pred.ggplot
  
  ## Plot of Capture Probabilities (regular scale) with length in regular scale
  Result$pred.ggplot_reg <- prediction.plots$pred.ggplot_reg
  
  #######################################################
  
  return(Result)
  
} # end of ic.fit.model




###############################################################################
###############################################################################
# "inv.fisher.info" function create the variance covariance matrix using the 
# negative log likelihood.
# Singular Value Decomposition is used to calculate the inverse of the hessian
#
# Input : vector of estimate in logit/log form and Data
#       
# Output: Variance Covariance Matrix (vcv) logit/log form


ic.inv.fisher.info = function(logit.est,Data, indicator,
                              captureDM , thetaDM,lambdaDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  hess.matx = hessian(ic.neg.log.likelihood, x=logit.est, method=("Richardson"),
                      Data=Data, indicator=indicator,
                      captureDM=captureDM, thetaDM=thetaDM, lambdaDM=lambdaDM,
                      captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
                      lambdaOFFSET=lambdaOFFSET)
  hess.svd = svd(hess.matx)  # singular value decomposition
  hess.svd$d = 1/hess.svd$d  # inverse of the singular values
  
  # variance covariance matrix (inverse of the hessian matrix)
  vcv = hess.svd$v %*% diag(hess.svd$d) %*% t(hess.svd$u) 
  
  return(vcv)
}  # end of inv.fisher.info


###############################################################################
ic.make.anticont.exp <- function(x,ncats){
  # make an expression to convert a continuation ratio to the original p 
  #      (This is used to find the category proportion)
  # start = starting index, ncr=number of continuation ratios
  # It will return a list of size ncr
  # x is the starting position
  start = x
  ncr   = ncats-1
  index = seq(start, by=1, length.out=ncr)
  
  myform = c(paste("~x",index,sep=""),"~1")
  for(i in 1:ncr){
    myform = paste(myform,c(rep(" ",i),paste("*(1-x",start:(start+ncr-i),")",sep="")))
  }
  myform
}  # end of make.anticont.exp




