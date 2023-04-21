############################################################################### 
#####                 "PSkSCR.BAYESIAN.summary.table "  for               #####
#####  Partial Stratification in k-Sample Capture-Recapture Experiments   #####
#####                     with known dead removals                        #####
###############################################################################

###############################################################################
###############################################################################
# Function :  PSkSCR.BAYESIAN.summary.table 
# inputs : output from the MCMC.MH function(mcmc_out), data, 
#          design matrices and offset vectors
# output : object result that contains mean and sd for the parameters in regular
#          form after after burn in and after thining

PSkSCR.BAYESIAN.summary.table <- function(mcmc_out ,data,
                                          captureDM,thetaDM,lambdaDM,p_lossDM,
                                          captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET){
  
  thinned.matrix <- matrix( ncol=length(mcmc_out$thinned.array[1,1,]), nrow =0)
  
  #combine chains (after burn in and after thining)
  for( i in 1:mcmc_out$n.chains ){
    thinned.matrix <- rbind(thinned.matrix, mcmc_out$thinned.array[ , i, ])
  }
  
  
  cats <- data$category
  ncats <- length(data$category)
  result <- NULL
  

  # thinned posterior samples of all chains together in a one matrix
  result$param_names <- colnames(thinned.matrix) 
  
  result$thinned.matrix <- thinned.matrix
  
  result$n.post.burnin <-  mcmc_out$n.post.burnin # number of MCMC post burn in
  result$np <- nrow(mcmc_out$beta.array)
  result$rate <- mcmc_out$rate # acceptance rate
  result$n.burn.in <- mcmc_out$n.burn.in # number of MCMC burn in
  result$n.thin <- mcmc_out$n.thin
  result$n.chains <- mcmc_out$n.chains
  
  ###################################################
  # compute summary statistics ( after burn in and after thining )
  # calculate means,sd, and percentiles of the posterior distribution
  # mean,sd, percentiles are data frames
  
  mean_thin_df <- adply(thinned.matrix,2,mean)
  
  sd_thin_df <- adply(thinned.matrix,2,sd)
  
  percent_vals <- c(0.01, 0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.99)
  percentiles_thin <-adply(thinned.matrix,2,quantile,probs = percent_vals)
    
  summary_stat_thin <- data.frame(result$param_names,mean_thin_df[,2], sd_thin_df[,2],
                                  percentiles_thin[,2:length(percentiles_thin)])
  colnames(summary_stat_thin) <- c("Parameter", "Mean", "SD",
                                   paste(percent_vals *100, "%", sep=""))
  
  result$summary_stat_thin <- summary_stat_thin
  ###############################################


  return(result)
  
} # end of PSkSCR.BAYESIAN.summary.table




