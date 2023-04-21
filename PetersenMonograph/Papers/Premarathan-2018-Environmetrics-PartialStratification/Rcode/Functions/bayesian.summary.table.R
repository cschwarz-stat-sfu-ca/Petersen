
# this file contains two functions
#  "unpack.parm.full" and "bayesian.summary.table"
 
###############################################################################
###############################################################################

# # Function "unpack.parm.full" unpack the logit/log sclae parameters to
# # regular parameters
# #   Input : estimates in logit/log scales , Data, design matrices and offset vectors
# #  Output : estimates in regular form ( capture probabilities, 
# #           category proportions(lambda), and sub-sample proportions(theta) )
# unpack.parm.full <- function(logit.est, Data,captureDM,thetaDM,lambdaDM,
#                              captureOFFSET,thetaOFFSET,lambdaOFFSET){
#   temp <- unpack.parm(logit.est, Data,captureDM,thetaDM,lambdaDM,
#                       captureOFFSET,thetaOFFSET,lambdaOFFSET)
#   return(as.vector(temp$full))
# }

###############################################################################
###############################################################################
# Function :  bayesian.output 
# inputs : output from the MCMC.MH function(mcmc_out), Data, 
#          design matrices and offset vectors
# output : object result that contains mean and sd for the parameters in regular
#          form after discarding the first half of the chain, trace plots, density 
#          plots for the parameters, # acceptance rate, etc..

bayesian.summary.table <- function(mcmc_out ,Data,captureDM,thetaDM,lambdaDM,
                                  captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  thinned.matrix <- matrix( ncol=length(mcmc_out$thinned.array[1,1,]), nrow =0)
  
  for( i in 1:mcmc_out$n.chains ){
    thinned.matrix <- rbind(thinned.matrix, mcmc_out$thinned.array[ , i, ])
  }
  
  # apply to each column of thin.beta.array and output a Data frame with regular scale
  # this data frame combines( row bind) results in each chain
#   regular_thin_df <- adply(mcmc_out$thinned.array,c(1,2),unpack.parm.full, Data=Data,
#                            captureDM=captureDM,thetaDM=thetaDM,
#                            lambdaDM=lambdaDM,captureOFFSET=captureOFFSET,
#                            thetaOFFSET=thetaOFFSET,lambdaOFFSET=lambdaOFFSET)
#     
  
  cats <- Data$category
  ncats <- length(Data$category)
  result <- NULL
  
#   ###############################################
#   # naming the parameters
#   
#   #names for beta parameters
#   betanames <- paste("beta[",1:mcmc_out$n.beta,"]",sep="")
#   
#   
#   #Names for capture probabilities
#   capnames <- c()
#   for(i in 1:ncats){
#     temp <-    c(paste("p_",i,cats, sep=""))
#     capnames <- c(capnames,temp)
#   }
#   
#   #Names for category proportions
#   lambdanames <-  c(paste("lambda_",cats, sep=""))
#   
#   
#   #Names for sab-sample proportions
#   thetanames <-   c("theta_1", "theta_2")
#   
#   #Name for population size
#   popname <-  c("N")
#   
#   
#   #names for category totals
#   
#   cat.tot.names <- paste("N_", cats,sep="")
#   
#   param_names <- c(betanames,capnames,lambdanames,thetanames,popname,cat.tot.names)
#   result$param_names <- param_names
#   ####################################################
  
  # thinned posterior samples of all chains together in a one matrix
  result$param_names <- colnames(thinned.matrix) 
  
  result$thinned.matrix <- thinned.matrix
  
  result$n.post.burnin <-  mcmc_out$n.post.burnin # number of MCMC post burn in
  result$np <- nrow(mcmc_out$beta.array)
  result$rate <- mcmc_out$rate # acceptance rate
  result$n.burn.in <- mcmc_out$n.burn.in # number of MCMC burn in
  result$n.thin <- mcmc_out$n.thin
  result$n.chains <- mcmc_out$n.chains
  
  #############################################################################
  # compute summary statistics ( after bun in and after thining )
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
  
} # end of bayesian.summary.table




