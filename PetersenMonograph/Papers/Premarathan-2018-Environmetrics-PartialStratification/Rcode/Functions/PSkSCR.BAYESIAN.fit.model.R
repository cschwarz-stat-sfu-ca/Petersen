###############################################################################
#####                      Bayesian Fit Model  for                       ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################

# Function "PSkSCR.BAYESIAN.fit.model" 
#   Inputs: model id, row data, design matrices and offset vectors for capture 
#           probabilities,category proportions, sub-sample proportions, 
#           n.updates, n.burn.in, n.thin,nstops,n.chains,alpha,priors,
#           initial.proposal.sigma,target.acc.rate


PSkSCR.BAYESIAN.fit.model <- function(model.id,data,
                               captureDM, thetaDM,lambdaDM,p_lossDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET,
                               n.post.burnin=60000, n.burn.in=40000, n.thin=50,
                               nstops=20,n.chains=3,alpha=0.05,
                               priors,
                               initial.proposal.sigma,target.acc.rate=0.3){
  
      
  Result = NULL
     
  # fit the MLE model and get the results
  MLE_model_res <- PSkSCR.fit.model(model.id ,data,
                                    captureDM, thetaDM,lambdaDM,p_lossDM,
                                    captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
  Result$MLE_model_res <- MLE_model_res
  
  # "MLE_model_res$rawdata" is the original data prepared for analysis
  # see the function "PSkSCR.fit.model" for more detail about data preparation
  data <- MLE_model_res$rawdata

  Result$data <- data
  
  #required indicator variables in matrix form
  indicator <- MLE_model_res$indicator
  
  #vector of beta parameters(as in MARK). these are the MLEs in logit/log scale
  beta <- MLE_model_res$est$logit.full.red  
                                        
  n.beta <- length(beta) # number of beta parameters in logit/log scale
  
                  update <- NULL # update list have elements as follows
  #sigma for normal proposal distribution for all beta parameters
   update$proposal.sigma <- matrix(initial.proposal.sigma,n.beta,1)
         update$proposed <- matrix(0,n.beta,1)  # proposal generated
         update$accepted <- matrix(0,n.beta,1)  # number of acceptance
  # target acceptance rate for all beta parameters
  update$target.acc.rate <- matrix(target.acc.rate,n.beta,1) 
  
  # object with design matrices for capture probabilities, theta, lambda and p_loss
         design <- NULL 
       design$p <- captureDM
   design$theta <- thetaDM
  design$lambda <- lambdaDM
  design$p_loss <- p_lossDM
  
  ####################################################
  # MCMC output using Metropolis Hastings method
  # function "PSkSKR.BAYESIAN.MCMC.MH" is in the file "PSkSKR.BAYESIAN.MCMC.MH.R" file
  mcmc_out <- PSkSCR.BAYESIAN.MCMC.MH(n.post.burnin, n.burn.in, n.thin, nstops, n.chains,
                                      beta,design, priors, update, data, indicator,
                                      captureDM,thetaDM,lambdaDM,p_lossDM,
                                      captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
  Result$thinned.array <- mcmc_out$thinned.array
  
  # summary table for all parameters( betaparameters, regular parameters and category totals)
  Result$MCMC$summary.tables  <- PSkSCR.BAYESIAN.summary.table(mcmc_out ,data,
                                                               captureDM,thetaDM,
                                                               lambdaDM,p_lossDM,
                                                               captureOFFSET,thetaOFFSET,
                                                               lambdaOFFSET,p_lossOFFSET) 
  
  #(thinned)samples from the posterior for beta parameters (for all chains)
  Result$thinned.beta.matrix <- Result$MCMC$summary.tables$thinned.matrix[ , (1:n.beta)]

  
   ########## Model information #########################
  
  # model identification
  Result$model.id = model.id
  
  # number of parameters
  Result$np = n.beta 
  
  # DIC value
  Result$DIC <- PSkSCR.BAYESIAN.DIC( Result$thinned.beta.matrix, data, indicator,
                                     captureDM,thetaDM,lambdaDM,p_lossDM,
                                     captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
  
  Result$indicator <- indicator
  
  #design matrices and offset vectors for capture probabilities, 
  # category proportions and sub-sample proportions
  Result$captureDM <- captureDM
  Result$thetaDM <- thetaDM 
  Result$lambdaDM <- lambdaDM
  Result$p_lossDM <- p_lossDM
  
  Result$captureOFFSET <- captureOFFSET
  Result$thetaOFFSET <- thetaOFFSET
  Result$lambdaOFFSET <- lambdaOFFSET
  Result$p_lossOFFSET <- p_lossOFFSET
  
  ## MCMC info ###
  Result$n.post.burnin <- n.post.burnin
  Result$n.burn.in <-n.burn.in
  Result$n.thin <- n.thin
  Result$nstops <- nstops
  Result$n.chains <- n.chains
  Result$target.acc.rate <- target.acc.rate
  Result$priors <- priors
  Result$alpha <- alpha
  
  
  #############################
  # trace plots, Posterior density plots, ACF plots, and Gelman Rubin Diagnose
  # plots with Rhat values
  Result$plots <- PSkSCR.BAYESIAN.plots(Result)

  #############################
  
  return(Result)
  
} # end of "PSkSCR.BAYESIAN.fit.model"

###############################################################################
###############################################################################