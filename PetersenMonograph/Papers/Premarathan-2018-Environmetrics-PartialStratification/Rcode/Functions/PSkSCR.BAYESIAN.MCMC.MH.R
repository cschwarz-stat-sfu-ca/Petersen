###############################################################################
#####             MCMC using Metropolis Hastings algorithm for           ######
#####                   "n.updates" number of iterations                 ######
#####   Bayesian - Partial Stratification in k-Sample Capture-Recapture  ######
#####                  Experiments with known dead removals              ######
###############################################################################

###############################################################################
###############################################################################
# Function : PSkSCR.BAYESIAN.MCMC.MH  - MCMC with Metropolis Hastings algorithm
#  Inputs : number of MCMC updates, n.burn.in,n.thin, parameters(in logit/log scale), 
#           design matrices(for capture probabilities, lambda and theta),
#           update list(sigma for normal proposal distribution, number of 
#           proposal generated,number of acceptance,target acceptance rate),
#           and Data
# Outputs : a list with a matrix of parameter values in logit/log form(posterior) 
#            of size n.beta x n.updates and a matrix of acceptance rate

PSkSCR.BAYESIAN.MCMC.MH <- function(n.post.burnin ,n.burn.in, n.thin, nstops, n.chains,
                                    beta, design, priors, update, data, indicator,
                                    captureDM, thetaDM, lambdaDM, p_lossDM,
                                    captureOFFSET, thetaOFFSET, lambdaOFFSET, p_lossOFFSET){
  
  ncats <- length(data$category)
  cat <- data$category
  st <-  nchar(data$history[1]) # number of sampling occasions
  
  n.updates <- n.burn.in + n.post.burnin # total number of iterations
  
  initial.beta <- beta
  n.beta <- length(beta) # number of beta parameters
    
  ### beta.array contains the MCMC samples for all the iterations ###
  ## since thinned samples are interested from n.post.burnin, no need to 
  ## have this huge array
  #beta.array <- array(0,c(n.updates,n.chains,n.beta)) 
   
  # if n.thin > 1 part of the output is saved out of n.post.burnin
  # if n.thin = 1 its an array with eliments c(n.post.burnin, n.chains, #parameters)
  thinned.array <- array( 0,  c(((n.post.burnin)/n.thin), n.chains, 
                                ( n.beta+((ncats*st)+ncats+st+st+1) +length(data$category) ) ))
                                            
  #acceptance rates in (n.updates/nstops) iterations
  rate <- array(0,c(n.beta,nstops,n.chains)) 
   
  for(ch in 1:n.chains ) {
    if(ch==1){
      beta <- initial.beta  # initial values for the parameters for first chain
                            # in logit/log scale
    }else{ beta <- initial.beta  +rnorm(1,0,0.2) } #initial values for the parameters
                                                   #for other chains in logit/log scale
    
    #beta.array[1,ch , ] <- beta
        
    kr <- 1 # keep track of column of rate matrix
    kt <- 1 # keep track of column of thin.beat.matrix 
  
    for( iter in 2:n.updates){
      out <- PSkSCR.BAYESIAN.MH.bayes.update(beta,design,priors,update,data, indicator)
      beta <- out$beta
      #beta.array[iter, ch ,] <- beta  # update beta.array
      update <- out$update
    
       
      # part of the output is saved
      if(iter == (n.burn.in + (kt * n.thin))){
        # unpack beta parameters to regular parameters
        reg.parm <- PSkSCR.unpack.parm(beta, data,
                                       captureDM,thetaDM,lambdaDM, p_lossDM,
                                       captureOFFSET,thetaOFFSET,lambdaOFFSET, p_lossOFFSET)$full
        # get the categoty totals ( i.e:  category proportion * Population size)
        N.lambda <- reg.parm[(ncats*st+1):(ncats*st +ncats) ] * reg.parm[length(reg.parm)]
        
        thinned.array[ kt,ch , ] <- c(beta,reg.parm,N.lambda)
        kt <- kt + 1 
      }
        
      if(iter== kr*n.updates/nstops){ 
        for(i in 1:n.beta){
          rate[i,kr,ch] <-  update$accepted[i]/ update$proposed[i]
          update$accepted[i] <- 0
          update$proposed[i] <- 0  
        
          if( iter <= n.burn.in ){
            if((rate[i,kr,ch] > (update$target.acc.rate[i]+ 0.05)) || 
                    (rate[i,kr,ch] < (update$target.acc.rate[i]- 0.05))){             
              # set up acceptance rate between 25% to 35% 
              # adjust the sigma for normal proposal distribution if we are in 
              # the burn in  iterations
              update$proposal.sigma[i] <- update$proposal.sigma[i]*
                                           rate[i,kr,ch]/update$target.acc.rate[i]
            }
          }
        
        } # end of for loop
        kr <-  kr + 1
      }
    }# end of for loop for n.updates
   
  }# end of for loop for n.chains
  
  #################################################
  # naming the parameters
  
  #names for beta parameters
  betanames <- paste("beta[",1:n.beta,"]",sep="")
  
  #Names for capture probabilities
  capnames <- c()
  for(i in 1:st){
    temp <-    c(paste("p_",i,cat, sep=""))
    capnames <- c(capnames,temp)
  }
  
  #Names for category proportions
  lambdanames <-  c(paste("lambda_",cat, sep=""))
  
  #Names for sab-sample proportions
  thetanames <-   paste("theta_",c(1:st),sep="")
  
  #Names for loss on capture proportions
  p_lossnanames <-  paste("loss_Time_",c(1:st), sep="")
  
  #Name for population size
  popname <-  c("N")
  
  #names for category totals
    cat.tot.names <- paste("N_", cat,sep="")
  
  param_names <- c(betanames,capnames,lambdanames,thetanames,p_lossnanames,popname,cat.tot.names)
  ####################################################
  #naming the chains
  chain.names <- paste("chain_", (1:n.chains), sep="")
  
  #names for the thinned.array
  dimnames(thinned.array) <- list( NULL,c(chain.names), c(param_names)) 

  MCMC.MH.out <- NULL
  MCMC.MH.out$n.beta <- n.beta # number of beta parameters
  #MCMC.MH.out$beta.array <- beta.array  # beta parameters for all iterations
  MCMC.MH.out$thinned.array <- thinned.array # beta parameters, regular parameters
                                             # and category totals for thinned samples
  MCMC.MH.out$rate <- rate  # acceptance rate for each (n.updates/nstops) iterations
  MCMC.MH.out$n.post.burnin <- n.post.burnin
  MCMC.MH.out$n.burn.in <-n.burn.in
  MCMC.MH.out$n.thin <- n.thin
  MCMC.MH.out$n.chains <- n.chains
  
  return(MCMC.MH.out)

} #end of MCMC.MH
 
###############################################################################
###############################################################################

