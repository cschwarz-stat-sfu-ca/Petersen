# Incomplete stratification code.
# These functions will not be visible to the user.
# At the moment, these functions only work for the full likelihood and not individual covariates or using
# other variables in the data frame

# I've copied over L's code with minor modifications to make it run here

#'
#' @importFrom reshape2 melt
#' @importFrom plyr llply
#' @importFrom stats optim
#' @importFrom MASS ginv
#' @importFrom rlang .data
#' @noRd

###############################################################################
# Function "IS.fit.model"
# Inputs: model id, row data, design matrices and offset vectors for capture
          # probabilities, category proportions and sub-sample proportions
#  Outputs: MLE , variance, SE of all the parameters( capture probabilities,
          # category proportions, sub-sample proportions, and population size),
          # and model information (id, number of parameters,
                                    # negative log-likelihood and AICc values)
          # table of capture histaries, observed and expected counts
          # residual plot(standardized residual vs expected counts) for model assessment

IS.fit.model = function(model.id ,Data,
                        captureDM,     thetaDM,     lambdaDM,
                        captureOFFSET, thetaOFFSET, lambdaOFFSET,
                        control.optim=list(trace=1,maxit=1000)){

  # deal with R package visible binding. See https://ggplot2.tidyverse.org/articles/ggplot2-in-packages.html

  Result = NULL
  ncats = length(Data$category) # number of categories

  Result$rawdata = Data # raw data set

  # Initial Estimates (in regular form)
  initial.est = IS.initial.estimates(Data)

  # initial.pack.parm.est are in logit/log form
  initial.pack.parm.est = IS.pack.parm(initial.est$full,Data,
                                    captureDM,thetaDM,lambdaDM)
  Result$initial.unpack.parm = IS.unpack.parm(initial.pack.parm.est, Data,
                                           captureDM,    thetaDM,    lambdaDM,
                                           captureOFFSET,thetaOFFSET,lambdaOFFSET)

  # initial value for the log(N)
  log.N.init = initial.pack.parm.est[length(initial.pack.parm.est)]

  # lower bound for optimization routine
  l.b =  c(rep(logit(0.000001),length(initial.pack.parm.est)-1),max(2,log.N.init-3))
  # upper bound for optimization routine
  u.b =  c(rep(logit(.999999),length(initial.pack.parm.est)-1),log.N.init+3)

  # optimize the negative log likelihood function(MLEs in logit/log form)
  res =stats::optim(par=initial.pack.parm.est, fn=IS.neg.log.likelihood, method="L-BFGS-B",
              lower=l.b, upper=u.b , control=control.optim,
              Data=Data,
              captureDM=captureDM, thetaDM=thetaDM,lambdaDM=lambdaDM,
              captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
              lambdaOFFSET=lambdaOFFSET)

  Result$res = res

  #extract the MLEs for the parameter (estimates in regular form)
  Result$est = IS.unpack.parm(Result$res$par,Data,captureDM,thetaDM,lambdaDM,
                           captureOFFSET,thetaOFFSET,lambdaOFFSET)

  # if all the DM have zero columns
  if((ncol(captureDM)== 0)  && (ncol(thetaDM)== 0) && (ncol(lambdaDM)==0)){
     hess = 1/Result$res$par
  } else{
        # Variance Covariance Matrix of logit/log parameters
       hess = IS.inv.fisher.info(Result$res$par,Data,
                              captureDM, thetaDM,lambdaDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET)
    }

  # reduced Variance Covariance Matrix of logit/log parameters
  Result$vcv$logit.full.red = hess

  ###############################################
  #find the variance covariance (vcv) of parameters in regular form

  #Go from reduced estimates to expanded estimates
  #browser()
  X = as.matrix(Matrix::bdiag(captureDM,     # capture recapture parameters in logit form
                    lambdaDM,      # lamda parameters in logit form
                    thetaDM,       # theta patameters in logit form
                    diag(1)))       # the log(N) parameter
  vcv.logit.full = X %*% Result$vcv$logit.full.red %*% t(X)
  Result$vcv$logit.full = vcv.logit.full
  Result$se$logit.full  = sqrt(diag(vcv.logit.full))

  # Go from logit(continution ratio)  to continuation ratio and then continuation
  # ratio to original p and  log(N) -> N
  # We use the deltamethod() function in the msm package
  # and so need to create expressions to do the delta method
  # First take anti-logits of all parameters except the last (which is anti-log)
  antilogit = plyr::llply(c(paste("~1/(1+exp(-x",1:(nrow(Result$vcv$logit.full)-1),"))",sep=""),
                       paste("~exp(x",nrow(Result$vcv$logit.full),")",sep="")),as.formula)
  vcv.cont = msm::deltamethod(antilogit, Result$est$logit.full,
                         Result$vcv$logit.full, ses=FALSE)

  Result$est$logit.full.antilogit = expit(Result$est$logit.full)
  Result$est$logit.full.antilogit[length(Result$est$logit.full.antilogit)] =
                     exp(Result$est$logit.full[length(Result$est$logit.full.antilogit)])


  anti.cont.formula = IS.make.anticont.exp((2*ncats+1), ncats)
  anti.cont.formula = unlist(anti.cont.formula)
  anti.cont.formula = c(paste("~x",seq(1,(2*ncats)),sep=""),anti.cont.formula,
                        paste("~x",seq((3*ncats),(3*ncats+2)),sep="") )
  anti.cont.formula = plyr::llply(anti.cont.formula,as.formula)
  vcv.full = msm::deltamethod(anti.cont.formula,  Result$est$logit.full.antilogit,
                         vcv.cont, ses=FALSE)

  ###############################################
  # extract the individual VCV's  and SE's of the parameters

  Result$vcv$p = vcv.full[1:(2*ncats),1:(2*ncats)]
  Result$se$p  = matrix(sqrt(diag(Result$vcv$p)), nrow=ncats,byrow =FALSE)

  Result$vcv$lambda = vcv.full[(2*ncats+1):(3*ncats),(2*ncats+1):(3*ncats)]
  Result$se$lambda = sqrt(diag(Result$vcv$lambda))

  Result$vcv$theta = vcv.full[(3*ncats+1):(3*ncats+2),(3*ncats+1):(3*ncats+2)]
  Result$se$theta = sqrt(diag(Result$vcv$theta))

  Result$vcv$N = vcv.full[nrow(vcv.full),nrow(vcv.full)]
  Result$se$N  = sqrt(Result$vcv$N)

  Result$se$full = sqrt(diag(vcv.full))

  ###############################################
  # Extract the variance and the SE's of the expected number of individuals
  # in each category.

  #estimated category proportions(lambda) and estimated population size(N)
  lambda.N.hat = Result$est$full[c((2*ncats+1):(3*ncats),nrow(vcv.full))]

  # vcv of estimated category proportions(lambda) and estimated population size(N)
  vcv.lambda.N.hat = vcv.full[c((2*ncats+1):(3*ncats),nrow(vcv.full)),
                              c((2*ncats+1):(3*ncats),nrow(vcv.full))]

  form=c(0)
  for(k in (1:ncats)){
    form[k] = c(paste("~x",k, "*x",length(lambda.N.hat),sep=""))
  }
  form = plyr::llply(form,as.formula)

  # vcv of expected number of fish in each category
  vcv.categoty.total = msm::deltamethod(form, lambda.N.hat, vcv.lambda.N.hat, ses=FALSE)

  # variance of expected number of fish in each category
  Result$variance$N_lambda = diag(vcv.categoty.total)

  # SE of expected number of fish in each category
  Result$se$N_lambda = sqrt(Result$variance$N_lambda)

  ########## Model information #########################


  # model identification
  Result$model.id = model.id

  # number of parameters
  Result$np = length(Result$res$par)

  #Negative log likeligood value
  Result$NLL = IS.neg.log.likelihood(Result$res$par,Data,
                                  captureDM,thetaDM,lambdaDM,
                                  captureOFFSET,thetaOFFSET,lambdaOFFSET)

  #corrected AIC information for the model
  Result$AICc = IS.AICc(Result$res$par,Data,captureDM,thetaDM,lambdaDM,
                     captureOFFSET,thetaOFFSET,lambdaOFFSET)

  # design matrices and offset vectors for capture probabilities, category
  # proportions and sub-sample proportions
  Result$captureDM = captureDM
  Result$thetaDM = thetaDM
  Result$lambdaDM = lambdaDM
  Result$captureOFFSET = captureOFFSET
  Result$thetaOFFSET = thetaOFFSET
  Result$lambdaOFFSET =  lambdaOFFSET

  #########  Model assessment using residual plot #########

  #Find expected counts for model1:
  MLE.results= NULL
  MLE.results$param = Result$est
  MLE.results$category = Data$category

  # exp.data is a list with  expected capture histories,
  # corresponding counts and categories and vatiance for each expected category
  exp.data = IS.expected.counts(MLE.results)

  # Obs.counts
  #  This produces the observed counts as same order of the expected capture
  #  history. If some histories are not present in the
  #  observed data, then their counts will be zero

  #sum of the counts for each of history in Data$History
  obs.counts = outer(exp.data$History, Data$History, '==') %*% Data$count
  residual = obs.counts - exp.data$counts # calculate the residuals

  # calculate the standardized residuals
  standardized.residual = residual/sqrt(exp.data$variance)

  # Output.df gives the Histories, expected and observed counts and residuals
  Output.df =data.frame(History = exp.data$History,Observed.Counts =obs.counts,
                        Expected.counts =exp.data$counts,
                        Residual = residual,
                        Standardized.Residuals = standardized.residual)

  #plot the standardized residuals
  res.plot = ggplot( Output.df,aes(x = .data$Expected.counts, y =.data$Standardized.Residuals,
                                   label=.data$History)) +
                geom_point(size = 3,colour="#990000") +
                xlab("Expected Counts") + ylab("Standardized Residuals")+
                ggtitle(paste("Residual Plot \n  model : ", model.id,sep=" " ))+
                theme(axis.text=element_text(size=10,face="bold"),
                      axis.title=element_text(size=14),
                      plot.title = element_text(size = 12,face="bold")) +
                geom_hline(yintercept=c(-1.96,0,1.96),linewidth=1,
                           linetype=c("dashed","solid","dashed"),
                           colour=c("#660000","black","#660000"))+
                geom_text(aes(label=.data$History),hjust=- 0.2, vjust=0.2)

  # capture Histories and obesred and expected counts
  Result$obs.exp.counts = Output.df
  Result$res.plot = res.plot   # residual plot

  return(Result)

} # end of fit.model

###############################################################################
###############################################################################
# "IS.inv.fisher.info" function create the variance covariance matrix using the
# negative log likelihood.
# Singular Value Decomposition is used to calculate the inverse of the hessian
#
# Input : vector of estimate in logit/log form and Data
#
# Output: Variance Covariance Matrix (vcv) logit/log form

IS.inv.fisher.info = function(logit.est,Data, captureDM , thetaDM,lambdaDM,
                           captureOFFSET,thetaOFFSET,lambdaOFFSET){

    hess.matx = numDeriv::hessian(IS.neg.log.likelihood, x=logit.est, method=("Richardson"),
                        Data=Data,
                        captureDM=captureDM,thetaDM=thetaDM,lambdaDM=lambdaDM,
                        captureOFFSET=captureOFFSET,thetaOFFSET=thetaOFFSET,
                        lambdaOFFSET=lambdaOFFSET)
    hess.svd = svd(hess.matx)  # singular value decomposition
    hess.svd$d = 1/hess.svd$d  # inverse of the singular values

    # variance covariance matrix (inverse of the hessian matrix)
    vcv = hess.svd$v %*% diag(hess.svd$d) %*% t(hess.svd$u)

    return(vcv)
}  # end of IS.inv.fisher.info


###############################################################################
IS.make.anticont.exp <- function(x,ncats){
  # make an expression to convert a continution ratio to the original p
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
}  # end of IS.IS.make.anticont.exp




###############################################################################
############  Test Function ###################################################
###############################################################################
# # Following files are needed to test the function "fit.model"
#
# source(file="get.data.R")
# source(file="create.DM.R")
# source(file="helper.functions.R")
#
# test.fit.model = function(){
#   Data = NULL
#   Data$History   = c("U0", "M0", "MM", "F0" ,"FF", "0M", "0F", "0U", "UU" )
#   Data$counts    = c(41,   16,  4,  23, 7, 12, 25, 63,9)
#   Data$category  = c("M","F")
#
#
#   #model identification"
#   model.id ="No restrictions"
#
#   #get the required desigm matrices
#   captureDM = create.DM(c(1,2,3,4)) # Design matrix for cap prob
#   thetaDM   = create.DM(c(1,2))     # Design matrix for theta
#   lambdaDM  = create.DM(c(1))       # Design matrix for lambda
#
#   #give the offset vectors(vectors of zeros should be given since no restriction)
#   captureOFFSET = c(0,0,0,0)
#   thetaOFFSET   = c(0,0)
#   lambdaOFFSET  = c(0)
#
#
#   test.fit = fit.model(model.id ,Data,captureDM, thetaDM,lambdaDM,
#                        captureOFFSET,thetaOFFSET,lambdaOFFSET)
#   return(test.fit)
# }
# test.fit.model()
#
# ##########  End Test Function ###############################################
# #############################################################################


###############################################################################
# The following function creates initial estimates for the optimization routine.

# Function "IS.initial.estimates"
# Input : Data
# Output : initial values for capture probabilities, category proportions,
#          sub-sample proportions and population size

IS.initial.estimates = function(Data){
    cats = Data$category

    # n1 is the  total number of individuals caught in  sample 1
    # This is calculated counting the histories without a 0 in the first
    #  digit of the history (e.g. U0,UU, M0,MM, ...)
    n1 = sum(Data$counts[ ! (substr(Data$History,1,1)=="0" )])

    # m2 is the number of marked individuals caught in sample 2 ( e.g. UU, MM, FF,...)
    m2 = n1 - sum(Data$counts[ (substr(Data$History,2,2)=="0" )])

    # n2 is the  total number of individuals caught in  sample 2
    # ( m2 + count the histories with a 0 in the first digit of the history)
    n2 = m2 + sum(Data$counts[ (substr(Data$History,1,1)=="0")])

    # n1.star is the number of individuals in sub-sample  in sample 1
    # ( subtract total count of the  historieswith a U in the first digit of
    #   the history from n1)
    n1.star = n1 - sum(Data$counts[substr(Data$History,1,1)=="U" ])


    # n2.star is the number of individuals in sub-sample  in sample 2
    # ( subtract total count of the  histories with  0U from histories with
    #   a 0 in the first digit of the history )
    n2.star =  sum(Data$counts[ (substr(Data$History,1,1)=="0")]) -
                                         sum(Data$counts[Data$History=="0U"])


    ############ Initial values ###################
    # initial value for the population size; this is the simple Lincoln Petersen estimate
    N.init = n1*n2/m2

    theta_1 = n1.star/n1 # initial value for sub-sample proportion in sample 1
    theta_2 = n2.star/(n2-m2) # initial value for sub-sample proportion in sample 2

    # initial values for capture probabilities
    sample.1.cap.prob = rep(n1/N.init , length(cats))  # at time 1
    sample.2.cap.prob = rep(n2/N.init , length(cats))  # at time 2

    # initial values for category proportions:
    #                 use equal category proportions as inital values
    lambda.init  = rep(1/length(cats),length(cats))

    # al initial capture probabilities as a vector
    cap.prob.init = c(sample.1.cap.prob,sample.2.cap.prob)
    theta.init = c(theta_1 ,theta_2) # initial sub-sample proportions

    initial.estimates = NULL
    initial.estimates$p = matrix(cap.prob.init,ncol=2, byrow = FALSE)
    colnames(initial.estimates$p) = c("time1", "time2")
    rownames(initial.estimates$p) = c(paste(cats))

    initial.estimates$lambda = lambda.init
    names( initial.estimates$lambda)= c(paste("lambda_",cats))

    initial.estimates$theta = theta.init
    names(initial.estimates$theta) =  c("theta_1", "theta_2")

    initial.estimates$N = N.init
    names(initial.estimates$N) = c("N")

    initial.estimates$full = c(cap.prob.init,lambda.init,theta.init,N.init)


    return(initial.estimates)

} # end of function "IS.initial.estimates"


###############################################################################
##########  Test Function #####################################################
###############################################################################
# # test the function IS.initial.estimates
#
#  test.IS.initial.estimates = function(){
#   Data = NULL
#   Data$History   = c("u0", "U0", "M0" ,"MM", "F0" ,"F0" ,"FF", "0M",
#                      "0F", "0U", "0U", "0U","UU", "UU" )
#   Data$counts    = c(20, 21,  16 ,  4,  3, 20,  7, 12, 25, 43,  2, 18,  5,  4)
#   Data$category  = c("M","F")
#   ##  here U0 = 41, UU= 9, M0 = 16, MM = 4, F0 = 23, FF= 7, 0M = 12,
#   ##       0F= 25, 0U = 63
#
#   return(IS.initial.estimates(Data))
#  }
#  test.IS.initial.estimates()
#
###############################################################################
########## End  Test Function #################################################
###############################################################################


###############################################################################
###############################################################################

#The function "pack.parm" creates initial estimates for the optimization routine.
# Input : initial estimates are in regular form(length is  (3*length(cats)+3)),
#         and raw data
# Output : vector in logit/log form(length is  (3*length(cats)+3)-1)

IS.pack.parm = function(init.est,Data,captureDM,thetaDM,lambdaDM){

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
    logit.capture.rate.hat   = IS.LS(captureDM, logit.capture.rate)
    }


  if(ncol(lambdaDM)==0){
    logit.cats.proprtion.hat = rep(NULL, nrow(lambdaDM))
  } else {
        #initial reduced estimates for lambda
        logit.cats.proprtion = lambda.to.cumulative.logit.lambda(lambda)
        # use least squares to estimate  based on lambdaDM
        logit.cats.proprtion.hat = IS.LS(lambdaDM,logit.cats.proprtion)
    }


  if(ncol(thetaDM)==0){
    logit.theta.hat = rep(NULL, nrow(thetaDM))
  } else {
    # use least squares to estimate the theta based on thetaDM
    logit.theta = logit(theta)
    # transform the estimates using the logit function
    logit.theta.hat = IS.LS(thetaDM,logit.theta)
    }


  # transform the initial value for the population size to log scale
  log.N = log(N)

  return(c(logit.capture.rate.hat,logit.cats.proprtion.hat,logit.theta.hat,log.N ))
}

###############################################################################
###############################################################################
#Following code is  to check the above function

# testing.pack.parm = function(){
#   temp1 = IS.initial.estimates(Data)
#   temp2 = IS.pack.parm(temp1,Data,captureDM,thetaDM,lambdaDM)
#   temp  = IS.unpack.parm(temp2,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
#                       thetaOFFSET,lambdaOFFSET)
#   temp
#
# }
#
# testing.pack.parm()
#
###############################################################################
###############################################################################


###############################################################################
###############################################################################
# function "IS.unpack.parm " extract the vector of parameters
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

IS.unpack.parm = function(est, Data,captureDM,thetaDM,lambdaDM,
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
  parm$lambda = as.vector(IS.cumulative.logit.lambda.to.lambda(parm$logit.lambda))
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

# testing.IS.unpack.parm = function(){
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
#   res2 = IS.unpack.parm(testparm2,data, capDM,thetaDM,lambdaDM,
#                      captureOFFSET,thetaOFFSET,lambdaOFFSET)  # for 2 categories
#   #res3 = IS.unpack.parm(testparm3,data) # for 3 categories
#   res2 # for 2 categories
#   #res3 # for 3 categories
# }
#
#
# testing.unpack.parm()
#
###############################################################################
###############################################################################


###############################################################################
###############################################################################
# Function IS.neg.log.likelihood
#  Input : parameters are in logit or log form, raw Data, Desigm Matrices and
#          OFFSET vectors
#          (estimates are in the order:logit.capture.rates(p),
#            logit.cumulative.category.proportion(lambda),
#            logit.sampling.fraction(theta), log.population.size(N))
#           (Data is a list with  capture history, counts, category)
#  Output: negative log likelihood velue

IS.neg.log.likelihood = function(logit.est, Data,
                              captureDM,thetaDM,lambdaDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET){

  parm   = IS.unpack.parm(logit.est,Data,
                       captureDM,thetaDM,lambdaDM,
                       captureOFFSET,thetaOFFSET,lambdaOFFSET)

  # N hat
  N = parm$N

  hist = Data$History       # capture history from raw data
  ct   = Data$counts        # counts from raw data
  cats = Data$category      # categories

  History = c(hist,"00")      # add "00" capture history to end of the
                              # vector of capture history from raw data

  counts  = c(ct,(N-sum(ct))) # add total animals did not catch(relating to
                              # capture history "00") to the end of
                              # the vector of counts from the raw data

  # call the function 'IS.log.prob.history'. These probabilities and the counts
  # are in the exact order related to each capture history.
  log.prob = IS.log.prob.history(parm,History,cats)


  # log-ikelihood has 2 parts. Ignore the count! for the observed counts
#  part1 = lgamma(N+1) - sum(lgamma(counts+1)) # factorial part
  part1 = lgamma(N+1) - lgamma(N-sum(ct)+1) # factorial part

  part2 = sum(counts * log.prob)

  # log-likelihhod
  LLH = part1 + part2

  # negative log-likelihood
  NLLH = - LLH  # -(log-likelihood),convert to negative value for optimization purpose

  return(NLLH)

} # end of "IS.neg.log.likelihood"

###############################################################################
###############################################################################
##### CJS functions to compute the probability of each history passed to it. ##

IS.log.prob.history <- function(parm, history, cats){
# This computes the log(probability) of each history based on the parms list
# of parameters and the vector of category labels.
#
# We assume that parms has been unpacked and has elements (at a minimum)
#  p - length(cat) x 2 matrix
#  lambda - proportions of categories in population (must add to 1)
#  theta  - probability that an unmarked fish captured is “categorized”
#
# cats is a vector of category codes. The value of “U” CANNOT be used as a valid
#   code as this indicates an unknown category. Upper and lower case categories
#   are ok, but the user must use caution.
#
# history - a vector of histories. Each history is a character string of length 2
#  e.g. M0 indicates a male captured at time 1 and never seen again
#
# This routine will compute the probability of history (00) if you pass this to this
# routine even though this history is unobservable.
#
# Returns a vector of the log(prob) of each history. If a history has probability of 0
# a 0 is returned rather than a -Inf or NA. [This will allow graceful handling
# of cases when certain parameters are fixed at 1 or 0


# We construct a (length(cat)+2) x (length(cat)+2) matrix representing all of the
# possible capture histories log(probabilites) for (cat, "0" ,”U”) combinations
  prob <- matrix(0, nrow=length(cats)+2, ncol=length(cats)+2)

# Categories are captured at both occasions (e.g. MM)
  prob[matrix(rep(1:length(cats),2),ncol=2)] <- parm$lambda*parm$p[,1]*parm$theta[1]*parm$p[,2]

# Categories captured at first occasion and then never seen (e.g M0)
  prob[matrix(c(1:length(cats),
                rep(length(cats)+1,
                    length(cats))),ncol=2)] <- parm$lambda*parm$p[,1]*parm$theta[1]*(1-parm$p[,2])

# Categories not captured at first occassion, but captured at time 2 (e.g. 0M)
  prob[matrix(c(rep(length(cats)+1,
                    length(cats)),1:length(cats)),
              ncol=2)] <- parm$lambda*(1-parm$p[,1])*parm$p[,2]*parm$theta[2]

# History (00)
  prob[length(cats)+1,length(cats)+1] <- sum(parm$lambda*(1-parm$p[,1])*(1-parm$p[,2]))

# History (U0)
  prob[length(cats)+2,
       length(cats)+1] <- sum(parm$lambda*parm$p[,1]*(1-parm$theta[1])*(1-parm$p[,2]))

# History (0U)
  prob[length(cats)+1,
       length(cats)+2] <- sum(parm$lambda*(1-parm$p[,1])*parm$p[,2]*(1-parm$theta[2]))

# History (UU)
  prob[length(cats)+2,
       length(cats)+2] <- sum(parm$lambda*parm$p[,1]*(1-parm$theta[1])*parm$p[,2])

  if(sum(prob)<.9999999){stop("Error in computing cell probabilities")}

# Now index each history into the matrix of probabilities and return the log(prob)
  index <- matrix(c(match(substr(history,1,1),c(cats,"0","U")),
                    match(substr(history,2,2),c(cats,"0","U"))),ncol=2)
  #print(index)
  log_p_hist <- log(prob[index]+(prob[index]==0))
  #print(log_p_hist)

  return(log_p_hist)
}  # end of IS.log.prob.history


###############################################################################
###############################################################################
# #Following code is to check the above functions
#
# # data and patameter estimates for checking functions
# testing.data <- function(){
#
#   data <- NULL
#   data$History  = c("U0","UU","M0","MM","F0","FF","0M","0F","0U")
#   data$counts    = c(41,9,16,4,23,7,12,25,63)
#   data$category = c("M","F")
#   str(data)
#   #print(data)
#
#   # Set up some values of the parameters for testing
#   parm = NULL
#   parm$p      = matrix(seq(.1,.4,.1),nrow=2)
#   parm$lambda = c(0.4, 0.6)
#   parm$theta  = c(0.3, 0.6)
#   parm$N      = 500
#   str(parm)
#   #print(parm)
#
#   return(list(parm,data))
# }
#
#
# # test  both 'neg.log.like' and 'IS.log.prob.history' functions
# # 'IS.log.prob.history'should produce sum of the probabilities equal to 1
#
# check.both.functions  = function(){
#   # Check the IS.log.prob.history function
#   parm = testing.data()[[1]]
#   data = testing.data()[[2]]
#
#   res1  = neg.lg.ld(parm,data)
#   print(paste("Negative log Likelihhod: " , res1 , sep = ""))
#
#   res2 = IS.log.prob.history(parm,c(data$History,"00"),data$category)
#   print(paste("sum of the probabilities: " , sum(exp(res2)), sep = "")) # to see add to 1
# }
#
# check.both.functions()
###############################################################################

###############################################################################
# ##### CJS function to check probability #####################################
#
# check.IS.log.prob.history <- function(){
# # Check the IS.log.prob.history function
#
#    data <- NULL
#    data$History <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U")
#    data$Count   <- c(41,9,16,4,23,7,12,25,63)
#    data$Category <- c("M","F")
#    str(data)
#    print(data)
#
#    # Set up some values of the parameters for testing
#    parm <- NULL
#    parm$p <- matrix(seq(.1,.4,.1),nrow=2)
#    parm$lambda <- c(0.4, 0.6)
#    parm$theta <- c(0.3, 0.6)
#    print(parm)
#
#    # Test out on all possible capture histories (including the 00 history)
#    res <- IS.log.prob.history(parm, c(data$History,"00"), data$Category)
#    sum(exp(res))  # see if add to 1
# }
#
# check.IS.log.prob.history()
#
#
###############################################################################
###############################################################################


# Numerical version of "logit" and "expit" functions

# p scale to logit
logit = function(prob){
  result =log(prob/(1-prob))
  return(result)
}

#antilogit(i.e. from logit to p scale)
expit = function(logit.prob){
  # if logit value is too small or too big then the regular value is -Inf or Inf
  # Therefore to address that issue, keep the logit value between -10 and 10
  #logit.prob =  max(-10, min(10, logit.prob)) # this dosent work for vector

  for(i in 1:length(logit.prob)){
    logit.prob[i] = max(-10, min(10, logit.prob[i]))
  }
  result = 1/(1+exp(-logit.prob))
  return(result)
}

##############################################################################
# To pack lambda parameters

# regular probabilities to cumulative logit probabilities
# Input : lambda vector of length(prob)
# Output: Vector of logit cumulatives of lambda. lenght is length(prob)-1
lambda.to.cumulative.logit.lambda = function(prob){
  result = regular.to.cumulative.logit.prob(prob)
  return(result)
}

# The following transformation preserves the "sum to one" constraint.
# lambda.to.cumulative.logit.lambda ;
# Input: A vector of categary proportions  of length = length(prob)
# Ouput: A vector of length = length(prob)-1

# Takes a vector of length(prob) Probabilities and returns a transformed vector with length(prob)-1 terms.
# Where the terms are caculated as follows;
#
# P_i  = logit( P_i/(P_i + ... P_k) ) ,    i = 1..k-1
#
# Assumptions: The vector prob sums to 1

regular.to.cumulative.logit.prob = function(prob){
  return(logit(cum.prob(prob)))
}

# The following function transform  the vector of size length(prob)  to form of a vactor of size length(prob)-1 vector
# for converting the lambda vector to cumulative form

cum.prob = function(prob){
  cum.prob = rev(cumsum(rev(prob)))
  cp  = prob/cum.prob
  cp[is.na(cp)]= 1
  return(cp[-length(cp)])
}

##############################################################################
# To unpack lambda parameters

# Inverts logit  probabilities back to regular probabiliries
IS.cumulative.logit.lambda.to.lambda = function(logit.prob){
  result = IS.inv.logit.lambda.to.lambda(logit.prob)
  #result[result < .00001] = 0  # force all small entries to be zero
  return(result)
}

# The following transformation preserves the "sum to one" constraint.
# IS.inv.logit.lambda.to.lambda;
# Input: A vector of logit.prob probabilities with length(logit.prob) elements;
# Ouput: A vector of probabilities "p" of length(logit.prob)+1; this vector sums to 1.
#
# The input vector elements are in the following form
# p_i = logit( P_i/(P_i + ... P_k) ) where, i = 1,...,(k-1)

IS.inv.logit.lambda.to.lambda = function(logit.prob){
  prob = expit(logit.prob)
  p = IS.inv.cum.prob(prob)
  return(p)
}

IS.inv.cum.prob = function(prob){
  temp = cumprod(1-prob)
  temp = c(1,temp[-length(temp)])*prob
  return(c(temp,1-sum(temp)))
}


# using  Least Squares method
#Input : a Matrix X(Design Matrix) and Vector y
#Output: Least Squares Solution
IS.LS = function(x,y){
  b = solve(t(x) %*% x) %*% t(x) %*% y
  return(b)
}


# using Generalized inverse method
#Input : a Matrix X(Design Matrix) and Vector y
#Output: generalized linear model Solution
GI = function(x,y){
  b = MASS::ginv(t(x) %*% x) %*% t(x) %*% y
  return(b)
}


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

IS.AICc <- function(logit.est, Data,captureDM,thetaDM,lambdaDM,
                 captureOFFSET,thetaOFFSET,lambdaOFFSET){
  NLL = IS.neg.log.likelihood(logit.est, Data,captureDM,thetaDM,lambdaDM,
                           captureOFFSET,thetaOFFSET,lambdaOFFSET)
   np = length(logit.est)  # number of parameters
  AIC = 2*np + 2*NLL  #  Akaike information criterion (AIC)
    n = sum(Data$counts)  # total captured
 AICc = AIC + 2*np*(np+1)/(n-np-1) # after the correction

 return(AICc)

}

###############################################################################
# testing the above function

# test.IS.AICc = function(Data){
#   temp1 = initial.estimates(Data)
#   temp2 = pack.parm(temp1,Data,captureDM,thetaDM,lambdaDM)
#   temp3 = IS.AICc(temp2,Data,captureDM,thetaDM,lambdaDM,
#                captureOFFSET,thetaOFFSET,lambdaOFFSET)
#   temp3
# }
#
# test.IS.AICc(Data)
###############################################################################


###############################################################################
###############################################################################
# "IS.expected.counts" generate the expected counts give the parameter set.
# Input : capture probabilities, category proportions,sub-sample proportions
#         and categotires in the data set
# output : expected counts for the all possible observable capture histories
#           (except the capture history "00"- animals not captured)
#        : variance for each expected count
IS.expected.counts = function(result){

   param  = result$param
   cats  = result$category

   expected.Data = NULL
   # Create capture histories for Categories captured at first occasion and then
   # never seen (e.g. M0)and  Categories are captured at both occasions (e.g. MM) and
   # Categories not captured at first occasion, but captured at time 2 (e.g. 0M)
   cap.hist = c(paste(cats,"0",sep=""),paste(cats,cats,sep=""),paste("0",cats,sep=""))

   # all possible observable capture histories("00" is unobservable)
   expected.Data$History = c("U0", "UU", cap.hist,"0U")
   expected.Data$category = cats

   # add "00" capture history to end of the vector of capture history of expected.data
   hist = c(expected.Data$History,"00")

   # the function "log.prob.history" is in the file "neg.log.likelihood" which was
   # used to find the log probabilities
   exp.cats = expected.Data$category
   expected.log.prob = IS.log.prob.history(param,hist,exp.cats)#call the function 'log.prob.history'.
   # These log probabilities and the counts are in the exact
   # order related to each capture history.

   # expected probabilities
   expected.prob = exp(expected.log.prob) - (expected.log.prob==0)
   # have to use "-(expected.log.prob==0)" because of  "prob[index]==0" used
   # inside the neg.log.likelihood
   ## sum(expected.prob) # this is just to check. Sum of the probabilities should equal to 1

   N =   param$N

   # expected counts for all the capture histories including history "00"
   expected.counts.all = N * expected.prob
   # variance for all the  capture histories including history "00"
   expected.variance.all = N*expected.prob*(1-expected.prob )
   # variance for  the capture histories without history "00"
   expected.Data$variance =  expected.variance.all[1:length(expected.variance.all)-1]
   #print.expected.counts =cat(paste("Expected count for the capture History",hist[1:length(hist)],
   #                          sep = " ", ":",expected.counts.all[1:length(expected.counts.all)],'\n'))

   # expected couts for all the capture
   expected.count = expected.counts.all[1:length(expected.counts.all)-1]
   ## histories without history "00"
   ## data frame for expected histories and counts
   ## df.expected.counts = data.frame(expected.Data$History , expected.count)
   ## names(df.expected.counts) = c("expected.History","expected.count")

   expected.Data$counts = expected.count

   return(expected.Data)

}


###############################################################################
########### Test function #####################################################
###############################################################################
#
# test.IS.expected.counts= function(){
#    data <- NULL
#    data$History <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U")
#    data$count   <- c(41,9,16,4,23,7,12,25,63)
#    data$category <- c("M","F")
#    str(data)
#    print(data)
#
#    # Set up some values of the parameters for testing
#    param <- NULL
#    param$p <- matrix(c( 0.1848103,0.1597991 ,0.1829310,0.2142439),nrow=2)
#    param$lambda <- c(0.3666383,(1-0.3666383))
#    param$theta <- c(0.5, 0.37)
#    param$N = 591
#
#    print(param)
#
#    test=NULL
#    test$param = param
#    test$category = data$category
#    return(IS.expected.counts(test))
# }
#
# test.IS.expected.counts()

###############################################################################
########### End Test function #################################################
###############################################################################

# print the results

IS.print.output = function(x){

  cats = x$rawdata$category
  ########## model information ####################################
  cat("\n")
  cat("Model information: \n")
  cat("\n")
  cat("Model Name: ", x$model.id,"\n")
  cat("Neg Log-Likelihood: ",x$NLL,"\n")
  cat("Number of Parameters: ",x$np,"\n")
  cat("AICc value:" ,x$AICc,"\n")
  cat("\n \n")

  ########### Raw Data ############################################
  cat("\n")
  cat("Raw data: \n")
  cat("\n")
  print(x$rawdata)
  cat("\n \n")
  ncats = length(x$rawdata$category)
  cat("\n \n")

  ########### initial values for optimization routine ##############
  cat("Initial values used for optimization routine: \n")
  cat("\n")
  cat("Initial capture probabilities: \n")
  print(x$initial.unpack.parm$p)
  cat("\n")
  cat("Initial category proportions: \n")
  print(x$initial.unpack.parm$lambda)
  cat("\n")
  cat("Initial sub-sample  proportions: \n")
  print(x$initial.unpack.parm$theta)
  cat("\n")
  cat("Initial population size for optimization(simple Lincoln Petersen estimator is used) :" , prettyNum(round(x$initial.unpack.parm$N), big.mark = ",") , "\n")
  cat("\n \n")

  ############################################################################################################
  ## Design matrices and offset vectors for capture probabilities, category proportions and sub-sample proportions
  cats =x$rawdata$category

  cat("Design matrix  and OFFSET vector for capture probabilities: \n")
  cap.DM = x$captureDM
  rownames(cap.DM)= c((paste("p1", cats,sep="")),(paste("p2", cats,sep="")))
  print(data.frame(Beta=cap.DM, OFFSET.vector=x$captureOFFSET))
  cat("\n")

  cat("Design matrix and OFFSET vector for sub-sample proportions (theta): \n")
  rownames(x$thetaDM)= c("theta_1","theta_2")
  print(data.frame(Beta=x$thetaDM, OFFSET.vector=x$thetaOFFSET))
  cat("\n")

  cat("Design matrix and OFFSET vector for category proportions (lambda):\n")
  rownames(x$lambdaDM)= c(paste("lambda_", cats[1:(length(cats)-1)],sep=""))
  print(data.frame(Beta=x$lambdaDM, OFFSET.vector=x$lambdaOFFSET))
  cat("\n")


  ############ Find the MLEs ##########################################################
  cat("\n")
  cat("Find MLEs: \n")
  cat("\n")
  cat("MLEs for capture probabilities: \n")
  print(x$est$p)
  cat("\n")
  cat("MLEs for category proportions: \n")
  print(x$est$lambda)
  cat("\n")
  cat("MLEs for sub-sample  proportions: \n")
  print(x$est$theta)
  cat("\n")
  cat("MLE for Population size : ", prettyNum(round(x$est$N), big.mark = ",")) #print number with commas

  category.total = prettyNum(round(x$est$N_lambda), big.mark = ",")
  cat(paste("\nMLE for Population size of category", cats,sep = " ", ":",category.total))



  cat("\n \n")

  ############# SE's of the above MLEs ##################################################
  cat("SE's of the MLEs \n")
  cat("\n")
  cat("SE's of the MLEs of capture probabilities: \n")
  se_p = x$se$p
  colnames(se_p) = c("time1", "time2")
  rownames(se_p) = c(paste(cats))
  print(se_p)

  cat("\n")
  cat("SE's of the MLEs of the category proportions: \n")
  se_lambad = x$se$lambda
  names( se_lambad)= c(paste("lambda_",cats))
  print(se_lambad)
  cat("\n")

  cat("SE's of the MLEs of the sub-sample  proportions: \n")
  se_theta = x$se$theta
  names(se_theta) =  c("theta_1", "theta_2")
  print(se_theta)
  cat("\n")

  cat("SE of the MLE of the population size: ", prettyNum(round(x$se$N), big.mark = ","))
  cat(paste("\nSE for Population size of category", cats, sep = " ", ":",prettyNum(round(x$se$N_lambda), big.mark = ",")))
  cat("\n")

  ######### Table of Histories, observed and expected counts and residuals ################
  cat("\n")
  cat(paste("Observed and Expected counts for capture histories for the model ", x$model.id,sep=" " ,'\n'))
  print(x$obs.exp.counts)

  ###########  Residual plot ##############################################################
  cat("\n")
  cat(paste("Standardized residual plot for the model", x$model.id,sep=" " ,'\n'))
  plot(x$res.plot)


} # end of "IS.print.output"






