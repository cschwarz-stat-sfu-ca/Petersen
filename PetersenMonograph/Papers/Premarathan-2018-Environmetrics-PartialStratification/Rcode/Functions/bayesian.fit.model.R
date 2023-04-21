# Function "bayesian.fit.model" 
#   Inputs: model id, row data, design matrices and offset vectors for capture 
#           probabilities,category proportions, sub-sample proportions, 
#           n.updates, n.burn.in, n.thin,nstops,n.chains,alpha,priors,
#           initial.proposal.sigma,target.acc.rate


bayesian.fit.model <- function(model.id,Data,captureDM,thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET,
                               n.post.burnin=60000, n.burn.in=40000, n.thin=50,
                               nstops=20,n.chains=3,alpha=0.95,
                               priors,
                               initial.proposal.sigma,target.acc.rate=0.3){
  
  Result = NULL
  Result$Data <- Data
   
  # fit the MLE model and get the results
  MLE_model_res <- fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
                             captureOFFSET,thetaOFFSET,lambdaOFFSET)
  Result$MLE_model_res <- MLE_model_res
  
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
  
  # object with design matrices for capture probabilities, theta and lambda
         design <- NULL 
       design$p <- captureDM
   design$theta <- thetaDM
  design$lambda <- lambdaDM
  
  ####################################################
  # MCMC output using Metropolis Hastings method
  # function "MCMC.MH" is in the file "bayesian.MCMC.MH.R" file
  mcmc_out <- MCMC.MH(n.post.burnin, n.burn.in,n.thin,nstops,n.chains,
                      beta,design,priors,update,Data,
                      captureDM,thetaDM,lambdaDM,
                      captureOFFSET,thetaOFFSET,lambdaOFFSET)
  Result$thinned.array <- mcmc_out$thinned.array
  
  # summary table for all parameters( betaparameters, regular parameters and category totals)
  Result$MCMC$summary.tables  <- bayesian.summary.table(mcmc_out ,Data,
                                      captureDM,thetaDM,lambdaDM,
                                      captureOFFSET,thetaOFFSET,lambdaOFFSET) 
  
  #(thinned)samples from the posterior for beta parameters (for all chains)
  Result$thinned.beta.matrix <- Result$MCMC$summary.tables$thinned.matrix[ , (1:n.beta)]

  
   ########## Model information #########################
  
  # model identification
  Result$model.id = model.id
  
  # number of parameters
  Result$np = n.beta 
  
  # DIC value
  Result$DIC <- bayesian.DIC( Result$thinned.beta.matrix,
                             Data,captureDM,thetaDM,lambdaDM,
                             captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  #design matrices and offset vectors for capture probabilities, 
  # category proportions and sub-sample proportions
  Result$captureDM = captureDM
  Result$thetaDM = thetaDM 
  Result$lambdaDM = lambdaDM
  Result$captureOFFSET = captureOFFSET
  Result$thetaOFFSET = thetaOFFSET
  Result$lambdaOFFSET =  lambdaOFFSET
  
  ## MCMC info ###
  Result$n.post.burnin <- n.post.burnin
  Result$n.burn.in <-n.burn.in
  Result$n.thin <- n.thin
  Result$nstops <- nstops
  Result$n.chains <- n.chains
  Result$target.acc.rate <- target.acc.rate
  Result$priors <- priors
  Result$alpha <- alpha
  
  
  #####################################################
  # trace plots, Posterior density plots, ACF plots, and Gelman Rubin Diagnose
  # plots with Rhat values
  Result$plots <- bayesian.plots(Result)
  
  ######################################################
  # Bayesian p-value scatter plots using the Discrepancy functions 
  # (a) deviance  (b) Freeman-Tukey (FT) statistic
  Result$bayesian.predictive.plots <- bayesian.p.value(Result)
  
  ######################################################
  
  
  return(Result)
  
}