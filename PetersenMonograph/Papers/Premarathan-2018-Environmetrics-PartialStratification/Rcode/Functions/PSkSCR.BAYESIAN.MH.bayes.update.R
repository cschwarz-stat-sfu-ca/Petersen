###############################################################################
#####   Bayesian - Partial Stratification in k-Sample Capture-Recapture  ######
#####                  Experiments with known dead removals              ######
###############################################################################

###############################################################################
###### Metropolis Hastings on each parameter in terms of using a symmetric ####
#######################   proposal distribution   #############################

# this file contains functions for calculations of prior densities and the
# function "PSkSCR.BAYESIAN.MH.bayes.update" for Metropolis Hastings algorithm.


###############################################################################
###############################################################################
# calculation of prior densities 

## normal distribution on the logit scale can be made to fit very closely to the
## chosen beta distribution. if x~N(mu,sd) in logit scale then distribution 
## function of exp(x)/(1 + exp(x)) in regular scale has a Beta(alpha,beta) distribution
#  

# this function returns the log(prior density) in logit/log scale
# Input x is in logit/log scale
log.prior.density <- function(x,mu,sigma){
    log.prior.d <- dnorm(x,mean=mu, sd = sigma,log=TRUE) # log(prior density)
    return(log.prior.d)
}



###############################################################################
###############################################################################

# Function : PSkSCR.BAYESIAN.MH.bayes.update
# inputs  :  beta - vector of parameters (on the logit/log scale)
#            design - design matrices for capture probabilities, lambda(category
#                     proportions) and theta(sub sample proportions)
#            priors - prior distributions
#            update - a list with proposal.sigma, generated proposed value, 
#                     number of acceptance and target acceptance rate
# outputs : list of the updated "beta" vector(logit/log scale) and the object "update"

PSkSCR.BAYESIAN.MH.bayes.update <- function(beta,design,priors,update,data,indicator){
  
  parmindex <- 1 # keep track of beta parameters
  
  #########################################
  # update the capture probabilities
  if(ncol(design$p)> 0) {
    for(i in seq_len(ncol(design$p))){
   
      oldbeta <- beta[parmindex]
    
      proposed_beta <- beta
      proposed_beta[parmindex] <- beta[parmindex] + 
                                    rnorm(1,0, update$proposal.sigma[parmindex])
        
      # proposal distribution is symmetric.  so Q_{ij}=Q_{ji}
    
      # numerator of the Metropolis Hasting ratio
      log_r_top <-  -(PSkSCR.neg.log.likelihood( proposed_beta, data, indicator,
                                                 captureDM,thetaDM,
                                                 lambdaDM, p_lossDM, 
                                                 captureOFFSET,thetaOFFSET,
                                                 lambdaOFFSET, p_lossOFFSET)) + 
                      log.prior.density( proposed_beta[parmindex],
                                            priors$p.mu[i],priors$p.sigma[i]) 
    
      # denominator of the Metropolis Hasting ratio
      log_r_bot <- -(PSkSCR.neg.log.likelihood( beta, data, indicator,
                                                captureDM,thetaDM,
                                                lambdaDM, p_lossDM, 
                                                captureOFFSET, thetaOFFSET,
                                                lambdaOFFSET, p_lossOFFSET)) +
                    log.prior.density(oldbeta,priors$p.mu[i],priors$p.sigma[i])
    
      MH_r <- exp(log_r_top - log_r_bot) # Metropolis Hastings ratio
    
      # make a decision
      if(runif(1) < min(MH_r,1)){
        # accept the move
        beta[parmindex] <- proposed_beta[parmindex]
        update$accepted[parmindex] <- update$accepted[parmindex] + 1
      }
    
      update$proposed[parmindex] <- update$proposed[parmindex] + 1
      parmindex <- parmindex + 1
    
    } 
  } # end of update the capture rates
  
  ############################################
  
  # update the category proportions (lambda)
  if(ncol(design$lambda) > 0) {
    for(i in seq_len(ncol(design$lambda))){
    
      oldbeta <- beta[parmindex]
    
      # proposal distribution is symmetric.  so Q_{ij}=Q_{ji}
      proposed_beta <- beta
      proposed_beta[parmindex] <- beta[parmindex] + 
                                   rnorm(1,0, update$proposal.sigma[parmindex])
      
    # numerator of the Metropolis Hasting ratio
      log_r_top <- -(PSkSCR.neg.log.likelihood( proposed_beta, data, indicator,
                                                captureDM,thetaDM,
                                                lambdaDM, p_lossDM, 
                                                captureOFFSET,thetaOFFSET,
                                                lambdaOFFSET, p_lossOFFSET)) + 
                                log.prior.density( proposed_beta[parmindex],
                                             priors$lambda.mu[i],priors$lambda.sigma[i]) 
    
      # denominator of the Metropolis Hasting ratio
      log_r_bot <- -(PSkSCR.neg.log.likelihood( beta, data, indicator,
                                                captureDM,thetaDM,
                                                lambdaDM, p_lossDM, 
                                                captureOFFSET, thetaOFFSET,
                                                lambdaOFFSET, p_lossOFFSET)) +
                                log.prior.density(oldbeta,priors$lambda.mu[i],
                                                        priors$lambda.sigma[i])
    
      MH_r <- exp(log_r_top - log_r_bot) # Metropolis Hastings ratio
    
      # make a decision
      if(runif(1) < min(MH_r,1)){
        # accept the move
        beta[parmindex] <- proposed_beta[parmindex]
        update$accepted[parmindex] <- update$accepted[parmindex] + 1
      }
    
      update$proposed[parmindex] <- update$proposed[parmindex] + 1
      parmindex <- parmindex + 1
    
    }
  }# end of  update the category proportions
  ############################################
  
  # update the sub-sample proportions (theta)
  if(ncol(design$theta) > 0){
    for(i in seq_len(ncol(design$theta))){
    
      oldbeta <- beta[parmindex]
    
      # proposal distribution is symmetric.  so Q_{ij}=Q_{ji}
      proposed_beta <- beta
      proposed_beta[parmindex] <- beta[parmindex] + 
                                  rnorm(1,0, update$proposal.sigma[parmindex])
    
      # numerator of the Metropolis Hasting ratio
      log_r_top <- -(PSkSCR.neg.log.likelihood( proposed_beta, data, indicator,
                                                captureDM,thetaDM,
                                                lambdaDM, p_lossDM, 
                                                captureOFFSET,thetaOFFSET,
                                                lambdaOFFSET, p_lossOFFSET)) + 
                              log.prior.density( proposed_beta[parmindex],
                                           priors$theta.mu[i],priors$theta.sigma[i] ) 
    
      # denominator of the Metropolis Hasting ratio
      log_r_bot <- -(PSkSCR.neg.log.likelihood( beta, data, indicator,
                                                captureDM,thetaDM,
                                                lambdaDM, p_lossDM, 
                                                captureOFFSET, thetaOFFSET,
                                                lambdaOFFSET, p_lossOFFSET)) +
                               log.prior.density(oldbeta,priors$theta.mu[i],
                                                        priors$theta.sigma[i])
    
      MH_r <- exp(log_r_top - log_r_bot) # Metropolis Hastings ratio
    
      # make a decision
      if(runif(1) < min(MH_r,1)){
        # accept the move
        beta[parmindex] <- proposed_beta[parmindex]
        update$accepted[parmindex] <- update$accepted[parmindex] + 1
      }
    
      update$proposed[parmindex] <- update$proposed[parmindex] + 1
      parmindex <- parmindex + 1
    
    } 
  } # end of update the category proportions
  ####################################################
  
  # update the loss on capture prbabilities (p_loss (nu))
  if(ncol(design$p_loss) > 0){
    for(i in seq_len(ncol(design$p_loss))){
      
      oldbeta <- beta[parmindex]
      
      # proposal distribution is symmetric.  so Q_{ij}=Q_{ji}
      proposed_beta <- beta
      proposed_beta[parmindex] <- beta[parmindex] + 
                                  rnorm(1,0, update$proposal.sigma[parmindex])
      
      # numerator of the Metropolis Hasting ratio
      log_r_top <- -(PSkSCR.neg.log.likelihood( proposed_beta, data, indicator,
                                                captureDM,thetaDM,
                                                lambdaDM, p_lossDM, 
                                                captureOFFSET,thetaOFFSET,
                                                lambdaOFFSET, p_lossOFFSET)) + 
                            log.prior.density( proposed_beta[parmindex],
                                      priors$p_loss.mu[i],priors$p_loss.sigma[i] ) 
      
      # denominator of the Metropolis Hasting ratio
      log_r_bot <- -(PSkSCR.neg.log.likelihood( beta, data, indicator,
                                                captureDM,thetaDM,
                                                lambdaDM, p_lossDM, 
                                                captureOFFSET, thetaOFFSET,
                                                lambdaOFFSET, p_lossOFFSET)) +
                            log.prior.density(oldbeta, priors$p_loss.mu[i],
                                                          priors$p_loss.sigma[i])
      
      MH_r <- exp(log_r_top - log_r_bot) # Metropolis Hastings ratio
      
      # make a decision
      if(runif(1) < min(MH_r,1)){
        # accept the move
        beta[parmindex] <- proposed_beta[parmindex]
        update$accepted[parmindex] <- update$accepted[parmindex] + 1
      }
      
      update$proposed[parmindex] <- update$proposed[parmindex] + 1
      parmindex <- parmindex + 1
      
    } 
  } # end of update the loss on capture probabilities
  ####################################################
    
  # update the population size(N)
    oldbeta <- beta[parmindex]
    
  # proposal distribution is symmetric.  so Q_{ij}=Q_{ji}
    proposed_beta <- beta
    proposed_beta[parmindex] <- beta[parmindex] + 
                                  rnorm(1,0, update$proposal.sigma[parmindex])
    
    # numerator of the Metropolis Hasting ratio
    log_r_top <- -(PSkSCR.neg.log.likelihood( proposed_beta, data, indicator,
                                              captureDM,thetaDM,
                                              lambdaDM, p_lossDM, 
                                              captureOFFSET,thetaOFFSET,
                                              lambdaOFFSET, p_lossOFFSET)) + 
                         log.prior.density( proposed_beta[parmindex],
                                            priors$N.mu,priors$N.sigma) 
    
    # denominator of the Metropolis Hasting ratio
    log_r_bot <- -(PSkSCR.neg.log.likelihood( beta, data, indicator,
                                              captureDM,thetaDM,
                                              lambdaDM, p_lossDM, 
                                              captureOFFSET, thetaOFFSET,
                                              lambdaOFFSET, p_lossOFFSET)) +
                         log.prior.density(oldbeta,priors$N.mu,priors$N.sigma)
    
    MH_r <- exp(log_r_top - log_r_bot) # Metropolis Hastings ratio
    
    # make a decision
    if(runif(1) < min(MH_r,1)){
      # accept the move
      beta[parmindex] <- proposed_beta[parmindex]
      update$accepted[parmindex] <- update$accepted[parmindex] + 1
    }
    
    update$proposed[parmindex] <- update$proposed[parmindex] + 1

    
  # end of update the population size(N)
  ####################################################
  
  return(list(beta=beta, update=update))
  
} # end of PSkSCR.BAYESIAN.MH.bayes.update

#############################################################################
#############################################################################
#############################################################################

#############################################################################
############ Testing the function "bayes.update"  ###########################
#############################################################################
# 
# # setwd("U:\\Lasantha/Research/Rcode")
# source("load.R")  # load required functions and packages
# 
# 
# test.MH.bayes.update <- function(){
# 
#            Data <- NULL
#    Data$History <- c("U0", "UU", "M0", "MM", "F0", "FF", "0M", "0F", "0U")
#     Data$counts <- c(40, 1, 5067,   40, 1551,   33,   41,  237, 3075)
#   Data$category <- c( "M", "F")
#
#   #model identification : unrestricted model
#   model.id <- paste("{ p(c*t), theta(t), lambda(c) }")
#   
#   # give  the required design matrices for capture recapture probabilities,
#   # theta(sampling(sexing) fractions ) and lambda(Category proportion)
#   captureDM <- create.DM(c(1,2,3,4)) 
#   thetaDM   <- create.DM(c(1,2)) 
#   lambdaDM  <- create.DM(c(1)) 
#   
#   #give the offset vectors(vectors of zero's should be given since no restriction)
#   captureOFFSET <- c(0,0,0,0) 
#   thetaOFFSET   <- c(0,0)
#   lambdaOFFSET  <- c(0)
#   
#   MLE_test <-  fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
#                             thetaOFFSET,lambdaOFFSET)
#   
#   beta <- MLE_test$est$logit.full  #vector of parameters(on the logit/log scale)
#                                        # these are the MLEs
#   n.beta <- length(beta) # number of parameters(on the logit/log scale)
#   
#                   update <- list() # update list have elements as follows
#    update$proposal.sigma <- matrix(0.10,n.beta,1) #sigma for normal proposal distribution
#          update$proposed <- matrix(0,n.beta,1)  # proposal generated
#          update$accepted <- matrix(0,n.beta,1)  # number of acceptance
#   update$target.acc.rate <- matrix(0.30,n.beta,1) # target acceptance rate is 30%
#   
#          design <- list() # object with design matrices
#        design$p <- captureDM
#    design$theta <- thetaDM
#   design$lambda <- lambdaDM
#   
#   # test the function "MH.bayes.update" 
#   out <- MH.bayes.update(beta,design,update,Data)
#   
#   return(list(MLE_test$est$logit.full,out))
#     
# }
# 
#
# test.MH.bayes.update()
#
###############################################################################
####################  end of testing  #########################################
###############################################################################

############ Following are used to find the values for the priors  ############   
############ in logit and log scale                                ############

# rdirichlet(1, c(1,2))
# 
# ########################
# # logit value of regular probablity
# logit(0.3298)
# 
# ## normal distribution on the logit scale can be made to fit very closely
# # to the chosen beta distribution.
# logitnorm = data.frame(rnorm(100000,mean=-0.7, sd=0.2))
# names(logitnorm) <- c("val")
# regular <-  data.frame(1/( 1+ exp(-logitnorm)))
# colMeans(regular)
# sd(regular[,1])
# names(regular) <- c("val")
# ggplot( regular,aes(x = val))+ geom_density()


#########################
# 
# # log value of the regular population size
# log(205000)
# 
# ## normal distribution in log scale
# 
# lognorm = data.frame(rnorm(100000,mean=12.25, sd=0.4))
# names(logitnorm) <- c("val")
# regular <-  data.frame(exp(lognorm))
# mean(regular[,1])
# sd(regular[,1])
# names(regular) <- c("val")
# 
# ggplot( regular,aes(x = val))+ geom_density()
# #########################
###############################################################################

