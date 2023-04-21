###############################################################################
#####                           Print results                            ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################
###############################################################################

PSkSCR.BAYESIAN.print.output <- function(res){
 
  ########## model information ################################################
  cat("\n")
  cat("Model information: \n")
  cat("\n")
  cat("Model Name: ", res$model.id,"\n")
  cat("Number of beta Parameters: ",res$np,"\n")
  cat("DIC value (using $pD ):" ,res$DIC$DIC$pD,"\n")
  cat("DIC value ((using $pv):" ,res$DIC$DIC$pv,"\n")
  cat("\n \n")
  
  ########### Raw Data ########################################################
  cat("\n")
  cat("Raw data: \n")
  cat("\n")
  print.data <- data.frame(cbind(res$data$history, res$data$counts))
  names(print.data) <- c("history", "counts")
  print(print.data)
  cat("\n")
  cat("Categories: ", res$data$category,"\n")
  cat("\n \n")
  ncats = length(res$data$category)
  cat("\n")
  
  ##############################################################################
  ## Design matrices and offset vectors for capture probabilities, 
  ## category proportions, sub-sample proportions and loss on capture proportion
  cats = res$data$category
  st = nchar(res$data$history[1]) # number of sampling occasions
  
  cat("Design matrix  and OFFSET vector for capture probabilities: \n")
  cap.DM = res$captureDM
  rownames(cap.DM)= c(paste(rep( c(paste("p_", c(1:st),sep="")), 
                                 each = length(cats)),rep(cats, st),sep=""))
  print(data.frame(Beta=cap.DM, OFFSET.vector=res$captureOFFSET))
  cat("\n")
  
  cat("Design matrix and OFFSET vector for sub-sample proportions (theta): \n")
  rownames(res$thetaDM)= c(paste("theta_",c(1:st), sep=""))
  print(data.frame(Beta=res$thetaDM, OFFSET.vector=res$thetaOFFSET))
  cat("\n")
  
  cat("Design matrix and OFFSET vector for loss on capture probabilities: \n")
  rownames(res$p_lossDM)= c(paste("loss_Time_",c(1:st), sep=""))
  print(data.frame(Beta=res$p_lossDM, OFFSET.vector=res$p_lossOFFSET))
  cat("\n")
  
  cat("Design matrix and OFFSET vector for category proportions (lambda):\n")
  rownames(res$lambdaDM)= c(paste("lambda_", cats[1:(length(cats)-1)],sep=""))
  print(data.frame(Beta=res$lambdaDM, OFFSET.vector=res$lambdaOFFSET))
  cat("\n")
  
  
  ########## MCMC information #################################################
  cat("\n")
  cat("MCMC information using Metropolis Hastings Algorithm: \n")
  cat("\n")
  cat("Number of post burn in iterations = ", res$n.post.burnin,"\n")
  cat("Number of burn in iterations = ", res$n.burn.in,"\n")
  cat("Thin = ", res$n.thin,"\n")
  cat("Number of Parameters = ",res$np,"\n")
  cat("Number of chains = ",res$n.chains, "\n")
  cat("Target acceptance Rate of the proposal for each chain = ",
      res$target.acc.rate, "\n" )
  cat("Number of Acceptance rate calculation steps for each cahin = ",
      res$nstops, "\n" )
  cat("(Each step has ", ((res$n.post.burnin + res$n.burn.in)/res$nstops) ,
      "iterations in acceptance rate calculation)", "\n")
  cat("\n")
  
  cat("Acceptance rate at the last step(adaptively change the acceptance","\n")
  cat("ratio upto burn in iterations) for beta parameters", "\n")
  last.step.rate<-res$MCMC$summary.tables$rate[,res$nstops,]
  chains <- paste("chain_", 1:res$n.chains,sep="")
  colnames(last.step.rate)<- chains
  rownames(last.step.rate) <- paste("beta[", 1:res$np, "]",sep="")
  
  print(last.step.rate)
  cat("\n \n")
  
  ############ Summary statistics of the posterior distributions ##############

  cat("Posterior summary values","\n")
  print(res$plots$posterior.summary)
  cat("\n \n")
  
  cat("Summary statistics of the posterior distributions
  (after thinning and after discarding the burn in ) \n")
  cat("\n")
  print(res$MCMC$summary.tables$summary_stat_thin )
  cat("\n \n")
  
    
  ############# trace plots, Posterior density plots, ACF plots, ##############
  #################  and Gelman Rubin Diagnose plots ##########################
  
  ##### plots for beta parameters  ##### 
  cat("Trace plots for the beta parmeters", "\n")
  print(res$plots$beta.trace.plots)
  cat("\n")
  cat("Posterior plots for the beta parmeters", "\n")
  print(res$plots$beta.posterior.plots)
  cat("\n")
  cat("Autocorrelation plots for the beta parmeters", "\n")
  print(res$plots$beta.acf.plots)
  cat("\n")
  cat("Gelman Rubin Diagnose plots", "\n")
  print(res$plots$BGR.diag.plots)
  cat("\n \n")
  
  ##### plots for capture probabilities  ##### 
  cat("Trace plots for the capture probabilities", "\n")
  print(res$plots$capture.trace.plots)
  cat("\n")
  cat("Posterior plots for the capture probabilities", "\n")
  print(res$plots$capture.posterior.plots)
  cat("\n")
  cat("Autocorrelation plots for the capture probabilities", "\n")
  print(res$plots$capture.acf.plots)
  cat("\n \n")
  
  ##### plots for category proportions(lambda)  ##### 
  cat("Trace plots for the category proportions", "\n")
  print(res$plots$lambda.trace.plots)
  cat("\n")
  cat("Posterior plots for the category proportions", "\n")
  print(res$plots$lambda.posterior.plots)
  cat("Autocorrelation plots for the category proportions", "\n")
  print(res$plots$lambda.acf.plots)
  cat("\n \n")
 
  ##### plots for sub-sample proportions(theta)  ##### 
  cat("Trace plots for the sub-sample proportions", "\n")
  print(res$plots$theta.trace.plots)
  cat("\n")
  cat("Posterior plots for the sub-sample proportions", "\n")
  print(res$plots$theta.posterior.plots)
  cat("\n")
  cat("Autocorrelation plots for the sub-sample proportions", "\n")
  print(res$plots$theta.acf.plots)
  cat("\n \n")
  
  ##### plots for loss on capture proportions  ##### 
  cat("Trace plots for the loss on capture proportions", "\n")
  print(res$plots$p_loss.trace.plots)
  cat("\n")
  cat("Posterior plots for the loss on capture proportions", "\n")
  print(res$plots$p_loss.posterior.plots)
  cat("\n")
  cat("Autocorrelation plots for loss on capture proportions", "\n")
  print(res$plots$p_loss.acf.plots)
  cat("\n \n")
  
  ##### plots for population size (N) and category total  ##### 
  cat("Trace plots for the population size", "\n")
  print(res$plots$pop.trace.plots)
  cat("\n")
  cat("Posterior plots for the population size", "\n")
  print(res$plots$pop.posterior.plots)
  cat("\n")
  cat("Autocorrelation plots for the population size", "\n")
  print(res$plots$pop.acf.plots)
  cat("\n \n")
  
  
} # end of PSkSCR.BAYESIAN.print.output

###############################################################################
###############################################################################

