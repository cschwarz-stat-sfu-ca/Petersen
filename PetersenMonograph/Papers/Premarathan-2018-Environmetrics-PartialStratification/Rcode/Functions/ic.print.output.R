###############################################################################
######### Print output  for the data with individual covariates ###############
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
############################################################################### 

# print the results

ic.print.output = function(x){

  cats <-  x$cats 
  ########## model information ###############################################
  cat("\n")
  cat("Model information: \n")
  cat("\n")
  cat("Model Name: ", x$model.id,"\n")
  cat("Neg Log-Likelihood: ",x$NLL,"\n")
  cat("Number of Parameters: ",x$np,"\n")
  cat("AICc value:" ,x$AICc,"\n")
  cat("\n \n")
  
  
  ########### initial values for optimization routine #########################
  cat("Initial values used for optimization routine (regular form): \n")
  cat("\n")
  cat("Initial capture probabilities: \n")
  print(x$initial.est$p)
  cat("\n")
  cat("Initial category proportions (lambda): \n")
  print(x$initial.est$lambda)
  cat("\n")
  cat("Initial sub-sample  proportions(theta): \n")
  print(x$initial.est$theta)
  cat("\n \n")
  
  
  #############################################################################

  cat("Design matrix for category proportions (lambda):\n")
  rownames(x$lambdaDM)= c(paste("lambda_", cats[1:(length(cats)-1)],sep=""))
  print(data.frame(Beta=x$lambdaDM))
  cat("\n")
  
  cat("Design matrix  for sub-sample proportions (theta): \n")
  rownames(x$thetaDM)= c("theta_1","theta_2")
  print(data.frame(Beta=x$thetaDM))
  cat("\n")
  
  cat("Initial Beta parameters for optimization (ie. parameters are in logit scale)
      (First", ncol(x$captureDM), "are for capture probabilities and last", ncol(x$thetaDM), "for thetas 
      and others are for lambdas( since ", ncol(x$lambdaDM)+1, "categories", ncol(x$lambdaDM), "beta parameters for lambdas) \n")
  print(x$initial.beta.est)
  cat("\n")
  
  ############ Find the MLEs and SE ###########################################

  cat("Find MLEs and SE: \n")
  cat("\n")
  cat("MLEs for Beta parameters for capture probabilities: \n")

  beta.p <- data.frame( cbind(paste("p:",c(colnames(x$captureDM)),sep=""),
                        round(x$est$logit.p,4), round(x$se$logit.p,5)))

  names(beta.p) <- c( "Beta_Pamrameter", "Estimate", "se" )
  df.beta.p <- format(beta.p,  justify = "left")
  print(df.beta.p)
  cat("\n \n")

  cat("MLEs forReal parameters for category proportions(lambda): \n")
  df.real.lambda <- data.frame(cbind(x$est$lambda,x$se$lambda))
  names(df.real.lambda) <- c( "Estimate", "se")
  print(df.real.lambda)
  cat("\n \n")

  cat("MLEs for Real parameters for sub-sample proportions(theta): \n")
  df.real.theta <- data.frame(cbind(x$est$theta,x$se$theta))
  names(df.real.theta) <- c( "Estimate", "se")
  print(df.real.theta)
  cat("\n \n")

  cat("Estimates for Population size and category populations : \n ")
  cat("\n")
  names.N <- c("Population size (N)", paste("- Population size of category", cats,sep = " "))
  df.real.N <- data.frame(cbind(names.N , prettyNum(round(c(x$est$N,x$est$N_lambda)), big.mark = ",")))
  names(df.real.N) <- c( "Description", "Estimate" )
  print(df.real.N)
  cat("\n \n")
  
  cat("SE of the estimate of population size(N) is sqrt(part1 + part2)\n ")
  cat("SE of the estimate of population size(N): ", prettyNum(round(x$se$N), big.mark = ","))
  cat("\n \n ")
  cat("varince part1 (S-sq) -  population size(N): ", prettyNum(round(x$VCV$s_sq ), big.mark = ","))
  cat("\n ")
  cat("varince part2 ( t(D)*inv(I) * D ) -  population size(N): ", prettyNum(round(x$VCV$Dt_Iinv_D ), big.mark = ","))
  cat("\n \n \n")
  
  
  ########### WITH  adjust.N ################## 
  cat("MLE for Population size and category populations with adjust N :
      Considered capture probability (phi) >0.0005\n ")
  cat("\n")
  names.N <- c("Population size (N)", paste("- Population size of category", cats,sep = " "))
  df.real.N <- data.frame(cbind(names.N , prettyNum(round(c(x$est$adjust.N ,x$est$adjust.N*x$est$lambda )), big.mark = ",")))
  names(df.real.N) <- c( "Description", "Estimate" )
  print(df.real.N)
  cat("\n \n")
  
  #############################################
    
  
  ## Predictive Plot of Capture Probabilities (logit scale and regular scale) ##
  cat("\n \n ")
  cat("Predictive Plot of Capture Probabilities (logit scale and regular scale) : \n")
  print(x$pred.logit.ggplot.p)

  print(x$pred.ggplot.p)
  
  print(x$pred.ggplot_reg) 
  
 
} # end of "print.output"

###############################################################################
