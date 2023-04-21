###############################################################################
#####                          Fit Model  for                            ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################

# Function : PSkSCR.fit.model
#  Inputs: model id, row data, design matrices and offset vectors for capture 
#          probabilities, category proportions and sub-sample proportions
#  Outputs: MLE , variance, SE of all the parameters( capture probabilities, 
#           category proportions, sub-sample proportions, and population size),
#           and model information (id, number of parameters,
#           negative log-likelihood and AICc values)
#           table of capture histories, observed and expected counts
#           residual plot(standardized residual vs expected counts) for model assessment



PSkSCR.fit.model <- function(model.id ,data,captureDM, thetaDM,lambdaDM,p_lossDM,
                             captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET,
                             dataprepare="NO"){
  
  #######################################
  ##### Preparing data ##################
  
  # each histories can be duplicated many times in the given data with different
  # for counts ( both negative and positive values). Data is prepare such that each 
  # history is duplicated maximum two times(if there is negative value for the 
  # counts) eg. if the history "0U0" is there many times in the original data file 
  # with both negative and positive values for corresponding counts, after preparing
  # the history "0U0" appears only two times( to represent the negative and 
  # positive counts)
  # The data preparation is needed so that it is easier to get the expected counts 
  # for the observed capture histories. Expected counts are required to create 
  # residual plots
  
  if(dataprepare=="NO"){
    # histories of the given data
    hist <- data$history 
    # corresponding counts for the  histories of the given data
    counts <- data$counts 
    # replicate each history to tally the corresponding counts
    individual.hist <- rep(hist,abs(counts))
    # create -1 an +1 to tally the corresponding counts
    sign.counts <- rep( sign(counts),abs(counts))
    
    # create a two-way table with unique histories and capture indicators
    # (eg. Capture histories might duplicate two times if there are both positive
    #      and negative  counts for each history )
    summary <-  table(individual.hist,sign.counts)
    
    data.summary <- as.data.frame(summary) # create a data frame
    names(data.summary) <- c("history", "cap.ind", "frequency")
    
    new.cap.ind <- rep(0, nrow(data.summary))
    new.cap.ind[ data.summary$cap.ind==-1] <- -1
    new.cap.ind[ data.summary$cap.ind!=-1] <- 1
    
    data.summary$counts <- new.cap.ind* data.summary$frequency
    
    obs.summary <- data.summary[data.summary$counts != 0, ]
    
    data$history <- as.character(obs.summary$history)
    data$counts <- obs.summary$counts
  }
  
  
  ###### End Preparing data ##############
  ########################################
  
  
  Result <- NULL
  ncats <- length(data$category) # number of categories
  st <-  nchar(data$history[1]) # number of sampling occasions
  Result$st <- st
  
  Result$rawdata <- data # raw data set
  
  # create required indicator variables in matrix form. These are required for 
  # calculating the likelihood
  indicator <- PSkSCR.create.indicator(data)
  Result$indicator <- indicator
  
  # Initial Estimates (in regular form)
  initial.est <- PSkSCR.initial.estimates(data,indicator)
  
  # initial.pack.parm.est are in logit/log form
  initial.pack.parm.est <- PSkSCR.pack.parm(initial.est$full,data,
                                            captureDM,thetaDM,lambdaDM, p_lossDM)  
  Result$initial.unpack.parm <- PSkSCR.unpack.parm(initial.pack.parm.est, data, 
                                                   captureDM, thetaDM,lambdaDM,p_lossDM,
                                                   captureOFFSET,thetaOFFSET,
                                                   lambdaOFFSET, p_lossOFFSET )
  

  # initial value for the log(N)
  log.N.init <- initial.pack.parm.est[length(initial.pack.parm.est)]  
  
  # lower bound for optimization routine
  l.b <- c(rep(logit(0.000001),length(initial.pack.parm.est)-1),max(2,log.N.init-2))
  # upper bound for optimization routine
  u.b <- c(rep(logit(.99999),length(initial.pack.parm.est)-1),log.N.init+2) 
  
  # optimize the negative log likelihood function(MLEs in logit/log form)
  res <- optim(par=initial.pack.parm.est, fn=PSkSCR.neg.log.likelihood, method="L-BFGS-B",
               lower=l.b, upper=u.b , control=list(trace=3,maxit=1000),
               data=data, indicator=indicator, 
               captureDM=captureDM, thetaDM=thetaDM,lambdaDM=lambdaDM,p_lossDM=p_lossDM,
               captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
               lambdaOFFSET=lambdaOFFSET,p_lossOFFSET=p_lossOFFSET)

  
  Result$res <- res
  
  #extract the MLEs for the parameter (estimates in regular form)
  Result$est <- PSkSCR.unpack.parm(Result$res$par,data,
                                   captureDM,thetaDM,lambdaDM,p_lossDM,
                                   captureOFFSET,thetaOFFSET,
                                   lambdaOFFSET,p_lossOFFSET)  
  
  # if all the DM have zero columns
  if((ncol(captureDM)== 0)  && (ncol(thetaDM)== 0) && (ncol(lambdaDM)==0) && (ncol(p_lossDM)== 0)){ 
    hess = 1/Result$res$par
  } else{
    # Variance Covariance Matrix of logit/log parameters
    hess = PSkSCR.inv.fisher.info(Result$res$par,data, indicator,
                                  captureDM, thetaDM,lambdaDM,p_lossDM,
                                  captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
                                  
        }
  
  # reduced Variance Covariance Matrix of logit/log parameters
  Result$vcv$logit.full.red = hess   
  
  ###############################################
  #find the variance covariance (vcv) of parameters in regular form
  
  #Go from reduced estimates to expanded estimates
  X = bdiag(captureDM,     # capture recapture parameters in logit form
            lambdaDM,      # lamda parameters in logit form
            thetaDM,       # theta patameters in logit form
            p_lossDM,      # loss on capture parameters in logit form
            diag(1))       # the log(N) parameter
  vcv.logit.full = X %*% Result$vcv$logit.full.red %*% t(X)
  Result$vcv$logit.full = vcv.logit.full
  Result$se$logit.full  = sqrt(diag(vcv.logit.full))
  
  # Go from logit(continuation ratio)  to continuation ratio and then continuation
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
  
  
  anti.cont.formula = PSkSCR.make.anticont.exp(((ncats*st)+1), ncats) 
  anti.cont.formula = unlist(anti.cont.formula)
  anti.cont.formula = c(paste("~x",seq(1,(ncats*st)),sep=""),anti.cont.formula,
                        paste("~x",seq(((ncats*st)+ncats),((ncats*st)+ncats+st+st)),sep="") )
  anti.cont.formula = llply(anti.cont.formula,as.formula)
  vcv.full = deltamethod(anti.cont.formula,  Result$est$logit.full.antilogit,
                         vcv.cont, ses=FALSE)
  
  ###############################################
  # extract the individual VCV's  and SE's of the parameters
  
  #VCV's  and SE's of capture probabilities
  Result$vcv$p = vcv.full[1:(ncats*st),1:(ncats*st)]
  Result$se$p  = matrix(sqrt(diag(Result$vcv$p)), nrow=ncats,byrow =FALSE)
  
  #VCV's  and SE's of category proportions
  Result$vcv$lambda = vcv.full[((ncats*st)+1):((ncats*st)+ncats),
                               ((ncats*st)+1):((ncats*st)+ncats)]
  Result$se$lambda = sqrt(diag(Result$vcv$lambda))
  
  #VCV's  and SE's of sab-sample proportions
  Result$vcv$theta = vcv.full[((ncats*st)+ncats+1):((ncats*st)+ncats+st),
                              ((ncats*st)+ncats+1):((ncats*st)+ncats+st)]
  Result$se$theta = sqrt(diag(Result$vcv$theta))
  
  #VCV's  and SE's of loss on capture probabilities
  Result$vcv$p_loss = vcv.full[((ncats*st)+ncats+st+1):((ncats*st)+ncats+st+st),
                               ((ncats*st)+ncats+st+1):((ncats*st)+ncats+st+st)]
  Result$se$p_loss = sqrt(diag(Result$vcv$p_loss))
  
  #VCV's  and SE's of population size
  Result$vcv$N = vcv.full[nrow(vcv.full),nrow(vcv.full)]
  Result$se$N  = sqrt(Result$vcv$N)
  
  # SE's of all the parameters
  Result$se$full = sqrt(diag(vcv.full))
  
  ###############################################
  # Extract the variance and the SE's of the expected number of individuals
  # in each category.
  
  #estimated category proportions(lambda) and estimated population size(N)
  lambda.N.hat = Result$est$full[c(((ncats*st)+1):((ncats*st)+ncats),nrow(vcv.full))]
  
  # vcv of estimated category proportions(lambda) and estimated population size(N) 
  vcv.lambda.N.hat = vcv.full[c(((ncats*st)+1):((ncats*st)+ncats),nrow(vcv.full)),
                              c(((ncats*st)+1):((ncats*st)+ncats),nrow(vcv.full))]
  
  form=c(0)
  for(k in (1:ncats)){
    form[k] = c(paste("~x",k, "*x",length(lambda.N.hat),sep=""))
  }
  form = llply(form,as.formula)
  
  # vcv of expected number of fish in each category 
  vcv.category.total = deltamethod(form, lambda.N.hat, vcv.lambda.N.hat, ses=FALSE)
  
  # variance of expected number of fish in each category 
  Result$variance$N_lambda = diag(vcv.category.total)
  
  # SE of expected number of fish in each category 
  Result$se$N_lambda = sqrt(Result$variance$N_lambda)
  
  ########## Model information #########################
  
  
  # model identification
  Result$model.id = model.id
  
  # number of parameters
  Result$np = length(Result$res$par) 
  
  #Negative log likelihood value
  Result$NLL = PSkSCR.neg.log.likelihood(Result$res$par,data,indicator,
                                         captureDM,thetaDM,lambdaDM,p_lossDM,
                                         captureOFFSET,thetaOFFSET,
                                         lambdaOFFSET,p_lossOFFSET)
  
  #corrected AIC information for the model
  Result$AICc = PSkSCR.AICc(Result$res$par,data,indicator,
                            captureDM,thetaDM,lambdaDM,p_lossDM,
                            captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
  
  # design matrices and offset vectors for capture probabilities, category 
  # proportions and sub-sample proportions
  Result$captureDM = captureDM
  Result$thetaDM = thetaDM 
  Result$lambdaDM = lambdaDM
  Result$p_lossDM = p_lossDM
  
  Result$captureOFFSET = captureOFFSET
  Result$thetaOFFSET = thetaOFFSET
  Result$lambdaOFFSET = lambdaOFFSET
  Result$p_lossOFFSET = p_lossOFFSET
  
  #########  Model assessment using residual plot #########
  
  #Find expected counts for model1:
  
  parm <- Result$est
  
  # exp.data is a list with  expected capture histories, 
  # corresponding counts and categories and variance for each expected category
  # This produces the expected counts as same order as observed counts, because 
  # data is prepared as in the top of this function  in order to make it  easier. 
  # all the observable histories those are not observed are grouped into "OTHER"
  # history.
  expected.data <- PSkSCR.expected.counts(parm, data,indicator)

 
  # Obs.counts
  # add 0 to data$ counts to tally with the history "OTHER"
  obs.counts = c(data$counts,0)
  residual = abs(obs.counts) - abs(expected.data$counts) # calculate the residuals
  
  # calculate the standardized residuals
  standardized.residual = residual/sqrt(expected.data$variance) 
  
  # group gives whether the coount is a Positive or Negative value
  group <- rep(c("alive"),length(expected.data$history))
  group[expected.data$counts < 0] <- c("dead")
  
  # Output.df gives the Histories, expected and observed counts and residuals
  Output.df <- data.frame(History = expected.data$history,Observed.counts =obs.counts,
                          Expected.counts =round(expected.data$counts,2), 
                          Residual = round(residual,2),
                          Standardized.Residuals = round(standardized.residual,2),
                          log.absolute.expected.counts = round(log10(abs(expected.data$counts)),2),
                          Group = group)
  
   #plot the standardized residuals
  res.plot <- ggplot( Output.df,aes(x = log.absolute.expected.counts, y =Standardized.Residuals,
                                    label=History)) +
    #geom_point(size = 3,colour="#990000") +
    geom_point(aes(shape = as.factor(Group)), show_guide = FALSE,size = 4,colour="#990000") +
    xlab("Absolute values of Expected Counts in log Scale") + ylab("Standardized Residuals")+
    ggtitle(paste("PSkSCR Analysis: Residual Plot \n  model : ", model.id,sep=" " ))+
    theme(axis.text=element_text(size=10,face="bold"),
          axis.title=element_text(size=16),
          plot.title = element_text(size = 14,face="bold")) + 
    geom_hline(yintercept=c(-1.96,0,1.96),size=1,
               linetype=c("dashed","solid","dashed"),
               colour=c("#660000","black","#660000"))+ 
    geom_text(aes(label=History),hjust= 0.5, vjust=-0.5)
  
  # capture Histories and obesred and expected counts
  Result$obs.exp.counts = Output.df 
  Result$res.plot = res.plot   # residual plot
  
  return(Result)
  
} # end of PSkSCR.fit.model

###############################################################################
###############################################################################
# "PSkSCR.inv.fisher.info" function create the variance covariance matrix using the 
# negative log likelihood.
# Singular Value Decomposition is used to calculate the inverse of the hessian
#
# Input : vector of estimate in logit/log form and Data
#       
# Output: Variance Covariance Matrix (vcv) logit/log form

PSkSCR.inv.fisher.info = function(logit.est,data, indicator,
                                  captureDM , thetaDM,lambdaDM,p_lossDM,
                                  captureOFFSET,thetaOFFSET,
                                  lambdaOFFSET,p_lossOFFSET){
  
  hess.matx = hessian(PSkSCR.neg.log.likelihood, x=logit.est, method=("Richardson"),
                      data=data, indicator=indicator,
                      captureDM=captureDM,thetaDM=thetaDM,
                      lambdaDM=lambdaDM,p_lossDM=p_lossDM,
                      captureOFFSET=captureOFFSET,thetaOFFSET=thetaOFFSET,
                      lambdaOFFSET=lambdaOFFSET,p_lossOFFSET=p_lossOFFSET)
  hess.svd = svd(hess.matx)  # singular value decomposition
  hess.svd$d = 1/hess.svd$d  # inverse of the singular values
  
  # variance covariance matrix (inverse of the hessian matrix)
  vcv = hess.svd$v %*% diag(hess.svd$d) %*% t(hess.svd$u) 
  
  return(vcv)
}  # end of PSkSCR.inv.fisher.info


###############################################################################
###############################################################################
PSkSCR.make.anticont.exp <- function(x,ncats){
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


###############################################################################
###############################################################################