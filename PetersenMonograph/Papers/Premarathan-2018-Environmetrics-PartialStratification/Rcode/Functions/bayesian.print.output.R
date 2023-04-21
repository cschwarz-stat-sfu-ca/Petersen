
# Print results

bayesian.print.output <- function(x){
 
  ########## model information ################################################
  cat("\n")
  cat("Model information: \n")
  cat("\n")
  cat("Model Name: ", x$model.id,"\n")
  cat("Number of beta Parameters: ",x$np,"\n")
  cat("DIC value:" ,x$IDC$DIC,"\n")
  cat("\n \n")
  
  ########### Raw Data ########################################################
  cat("\n")
  cat("Raw data: \n")
  cat("\n")
  print(x$Data) 
  cat("\n \n")
  ncats = length(x$Data$category)
  cat("\n")
  
  ##############################################################################
  ## Design matrices and offset vectors for capture probabilities, 
  ## category proportions and sub-sample proportions
  cats =x$Data$category
  
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
  
  ########## MCMC information #################################################
  cat("\n")
  cat("MCMC information using Metropolis Hastings Algorithm: \n")
  cat("\n")
  cat("Number of post burn in iterations = ", x$n.post.burnin,"\n")
  cat("Number of burn in iterations = ", x$n.burn.in,"\n")
  cat("Thin = ", x$n.thin,"\n")
  cat("Number of Parameters = ",x$np,"\n")
  cat("Number of chains = ",x$n.chains, "\n")
  cat("Target acceptance Rate of the proposal for each chain = ",
      x$target.acc.rate, "\n" )
  cat("Number of Acceptance rate calculation steps for each cahin = ",
      x$nstops, "\n" )
  cat("(Each step has ", ((x$n.post.burnin + x$n.burn.in)/x$nstops) ,
      "iterations in acceptance rate calculation)", "\n")
  cat("\n")
  
  cat("Acceptance rate at the last step(adaptively change the acceptance","\n")
  cat("ratio upto burn in iterations) for beta parameters", "\n")
  last.step.rate<-x$MCMC$summary.tables$rate[,x$nstops,]
  chains <- paste("chain_", 1:x$n.chains,sep="")
  colnames(last.step.rate)<- chains
  rownames(last.step.rate) <- paste("beta[", 1:x$np, "]",sep="")
  
  print(last.step.rate)
  cat("\n \n")
  
  ############ Summary statistics of the posterior distributions ##############

  cat("Posterior summary values","\n")
  print(x$plots$posterior.summary)
  cat("\n \n")
  
  cat("Summary statistics of the posterior distributions
(after thinning and after discarding the burn in ) \n")
  cat("\n")
  print(x$MCMC$summary.tables$summary_stat_thin )
  cat("\n \n")
  
    
  ############# trace plots, Posterior density plots, ACF plots, ##############
  #################  and Gelman Rubin Diagnose plots ##########################
  
  ##### plots for beta parameters  ##### 
  cat("Trace plots for the beta parmeters", "\n")
  print(x$plots$beta.trace.plots)
  cat("\n")
  cat("Posterior plots for the beta parmeters", "\n")
  print(x$plots$beta.posterior.plots)
  cat("\n")
  cat("Autocorrelation plots for the beta parmeters", "\n")
  print(x$plots$beta.acf.plots)
  cat("\n")
  cat("Gelman Rubin Diagnose plots", "\n")
  print(x$plots$BGR.diag.plots)
  cat("\n \n")
  
  ##### plots for capture probabilities  ##### 
  cat("Trace plots for the capture probabilities", "\n")
  print(x$plots$capture.trace.plots)
  cat("\n")
  cat("Posterior plots for the capture probabilities", "\n")
  print(x$plots$capture.posterior.plots)
  cat("\n")
  cat("Autocorrelation plots for the capture probabilities", "\n")
  print(x$plots$capture.acf.plots)
  cat("\n \n")
  
  ##### plots for category proportions(lambda)  ##### 
  cat("Trace plots for the category proportions", "\n")
  print(x$plots$lambda.trace.plots)
  cat("\n")
  cat("Posterior plots for the category proportions", "\n")
  print(x$plots$lambda.posterior.plots)
  cat("Autocorrelation plots for the category proportions", "\n")
  print(x$plots$lambda.acf.plots)
  cat("\n \n")
 
  ##### plots for sub-sample proportions(theta)  ##### 
  cat("Trace plots for the sub-sample proportions", "\n")
  print(x$plots$theta.trace.plots)
  cat("\n")
  cat("Posterior plots for the sub-sample proportions", "\n")
  print(x$plots$theta.posterior.plots)
  cat("\n")
  cat("Autocorrelation plots for the sub-sample proportions", "\n")
  print(x$plots$theta.acf.plots)
  cat("\n \n")
  
  ##### plots for population size (N) and category total  ##### 
  cat("Trace plots for the population size", "\n")
  print(x$plots$pop.trace.plots)
  cat("\n")
  cat("Posterior plots for the population size", "\n")
  print(x$plots$pop.posterior.plots)
  cat("\n")
  cat("Autocorrelation plots for the population size", "\n")
  print(x$pop.acf.plots)
  cat("\n \n")
  
  
  ########## Predictive plots and Bayesian p-value ############################
  cat("Bayesian p-value and Predictive plots", "\n")
  cat("Bayesian p-value using Deviance = ", 
      x$bayesian.predictive.plots$deviance.p.value, "\n") 
  cat("Bayesian p-value using Freeman-Tukey statistic = ",
      x$bayesian.predictive.plots$FT.p.value, "\n")
  print(x$bayesian.predictive.plots$predictive.plots)
  
  
} # end of bayesian.print.output


