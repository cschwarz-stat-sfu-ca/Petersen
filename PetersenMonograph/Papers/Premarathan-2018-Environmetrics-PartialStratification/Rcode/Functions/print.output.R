# print the results
 
print.output = function(x){
  
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
  x$res.plot  
    
  
} # end of "print.output"


