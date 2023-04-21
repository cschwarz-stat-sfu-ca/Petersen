###############################################################################
# Function "fit.model" 
# Inputs: model id, row data, design matrices and offset vectors for capture 
          # probabilities, category proportions and sub-sample proportions
#  Outputs: MLE , variance, SE of all the parameters( capture probabilities, 
          # category proportions, sub-sample proportions, and population size),
          # and model information (id, number of parameters,
                                    # negative log-likelihood and AICc values)
          # table of capture histaries, observed and expected counts
          # residual plot(standardized residual vs expected counts) for model assessment

fit.model = function(model.id ,Data,captureDM, thetaDM,lambdaDM,
                     captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  Result = NULL
  ncats = length(Data$category) # number of categories
  
  Result$rawdata = Data # raw data set
  
  # Initial Estimates (in regular form)
  initial.est = initial.estimates(Data)
  
  # initial.pack.parm.est are in logit/log form
  initial.pack.parm.est = pack.parm(initial.est$full,Data,
                                    captureDM,thetaDM,lambdaDM)  
  Result$initial.unpack.parm = unpack.parm(initial.pack.parm.est, Data, 
                                           captureDM, thetaDM,lambdaDM,
                                           captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  # initial value for the log(N)
  log.N.init = initial.pack.parm.est[length(initial.pack.parm.est)]  
  
  # lower bound for optimization routine
  l.b =  c(rep(logit(0.000001),length(initial.pack.parm.est)-1),max(2,log.N.init-3))
  # upper bound for optimization routine
  u.b =  c(rep(logit(.999999),length(initial.pack.parm.est)-1),log.N.init+3) 

  # optimize the negative log likelihood function(MLEs in logit/log form)
  res =optim(par=initial.pack.parm.est, fn=neg.log.likelihood, method="L-BFGS-B",
              lower=l.b, upper=u.b , control=list(trace=1,maxit=1000),
              Data=Data, 
              captureDM=captureDM, thetaDM=thetaDM,lambdaDM=lambdaDM,
              captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
              lambdaOFFSET=lambdaOFFSET)
  
  Result$res = res
  
  #extract the MLEs for the parameter (estimates in regular form)
  Result$est = unpack.parm(Result$res$par,Data,captureDM,thetaDM,lambdaDM,
                           captureOFFSET,thetaOFFSET,lambdaOFFSET)  
  
  # if all the DM have zero columns
  if((ncol(captureDM)== 0)  && (ncol(thetaDM)== 0) && (ncol(lambdaDM)==0)){ 
     hess = 1/Result$res$par
  } else{
        # Variance Covariance Matrix of logit/log parameters
       hess = inv.fisher.info(Result$res$par,Data, 
                              captureDM, thetaDM,lambdaDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET) 
    }
  
  # reduced Variance Covariance Matrix of logit/log parameters
  Result$vcv$logit.full.red = hess   
  
  ###############################################
  #find the variance covariance (vcv) of parameters in regular form
  
  #Go from reduced estimates to expanded estimates
  X = bdiag(captureDM,     # capture recapture parameters in logit form
            lambdaDM,      # lamda parameters in logit form
            thetaDM,       # theta patameters in logit form
            diag(1))       # the log(N) parameter
  vcv.logit.full = X %*% Result$vcv$logit.full.red %*% t(X)
  Result$vcv$logit.full = vcv.logit.full
  Result$se$logit.full  = sqrt(diag(vcv.logit.full))
 
  # Go from logit(continution ratio)  to continuation ratio and then continuation
  # ratio to original p and  log(N) -> N
  # We use the deltamethod() function in the msm package
  # and so need to create expressions to do the delta method
  # First take anti-logits of all parameters except the last (which is anti-log)
  antilogit = llply(c(paste("~1/(1+exp(-x",1:(nrow(Result$vcv$logit.full)-1),"))",sep=""),
                       paste("~exp(x",nrow(Result$vcv$logit.full),")",sep="")),as.formula)
  vcv.cont = deltamethod(antilogit, Result$est$logit.full,
                         Result$vcv$logit.full, ses=FALSE)
  
  Result$est$logit.full.antilogit = expit(Result$est$logit.full)
  Result$est$logit.full.antilogit[length(Result$est$logit.full.antilogit)] = 
                     exp(Result$est$logit.full[length(Result$est$logit.full.antilogit)])
  
  
  anti.cont.formula = make.anticont.exp((2*ncats+1), ncats) 
  anti.cont.formula = unlist(anti.cont.formula)
  anti.cont.formula = c(paste("~x",seq(1,(2*ncats)),sep=""),anti.cont.formula,
                        paste("~x",seq((3*ncats),(3*ncats+2)),sep="") )
  anti.cont.formula = llply(anti.cont.formula,as.formula)
  vcv.full = deltamethod(anti.cont.formula,  Result$est$logit.full.antilogit,
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
  form = llply(form,as.formula)
  
  # vcv of expected number of fish in each category 
  vcv.categoty.total = deltamethod(form, lambda.N.hat, vcv.lambda.N.hat, ses=FALSE)
  
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
  Result$NLL = neg.log.likelihood(Result$res$par,Data,
                                  captureDM,thetaDM,lambdaDM,
                                  captureOFFSET,thetaOFFSET,lambdaOFFSET)

  #corrected AIC information for the model
  Result$AICc = AICc(Result$res$par,Data,captureDM,thetaDM,lambdaDM,
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
  exp.data = expected.counts(MLE.results)  
  
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
  res.plot = ggplot( Output.df,aes(x = Expected.counts, y =Standardized.Residuals,
                                   label=History)) +
                geom_point(size = 3,colour="#990000") +
                xlab("Expected Counts") + ylab("Standardized Residuals")+
                ggtitle(paste("Residual Plot \n  model : ", model.id,sep=" " ))+
                theme(axis.text=element_text(size=10,face="bold"),
                      axis.title=element_text(size=14),
                      plot.title = element_text(size = 12,face="bold")) + 
                geom_hline(yintercept=c(-1.96,0,1.96),size=1,
                           linetype=c("dashed","solid","dashed"),
                           colour=c("#660000","black","#660000"))+ 
                geom_text(aes(label=History),hjust=- 0.2, vjust=0.2)
  
  # capture Histories and obesred and expected counts
  Result$obs.exp.counts = Output.df 
  Result$res.plot = res.plot   # residual plot
  
  return(Result)
  
} # end of fit.model

###############################################################################
###############################################################################
# "inv.fisher.info" function create the variance covariance matrix using the 
# negative log likelihood.
# Singular Value Decomposition is used to calculate the inverse of the hessian
#
# Input : vector of estimate in logit/log form and Data
#       
# Output: Variance Covariance Matrix (vcv) logit/log form

inv.fisher.info = function(logit.est,Data, captureDM , thetaDM,lambdaDM,
                           captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
    hess.matx = hessian(neg.log.likelihood, x=logit.est, method=("Richardson"),
                        Data=Data, 
                        captureDM=captureDM,thetaDM=thetaDM,lambdaDM=lambdaDM,
                        captureOFFSET=captureOFFSET,thetaOFFSET=thetaOFFSET,
                        lambdaOFFSET=lambdaOFFSET)
    hess.svd = svd(hess.matx)  # singular value decomposition
    hess.svd$d = 1/hess.svd$d  # inverse of the singular values
    
    # variance covariance matrix (inverse of the hessian matrix)
    vcv = hess.svd$v %*% diag(hess.svd$d) %*% t(hess.svd$u) 
    
    return(vcv)
}  # end of inv.fisher.info
 

###############################################################################
make.anticont.exp <- function(x,ncats){
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
}  # end of make.anticont.exp




###############################################################################
############  Test Function ###################################################
###############################################################################
# # Following files are needed to test the function "fit.model" 
# 
# source(file="AICc.R")
# source(file="get.data.R")
# source(file="create.DM.R")
# source(file="expected.counts.R")
# source(file="helper.functions.R")
# source(file="initial.estimates.R")
# source(file="neg.log.likelihood.R")
# source(file="unpack.parm.R")
# source(file="pack.parm.R")
# source(file="print.output.R")
# source(file="unpack.parm.R")
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
