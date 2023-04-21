###############################################################################
########### Bayesian Analysis of the Mille Lacs Walleye dataset ###############
##############################################################################
####  Metropolis Hastings on each parameter in terms of using a symmetric #####
####  proposal distributions                                              #####            
###############################################################################
# We can create different models and  the results are stored in
# bayes_model_1, bayes_model_2...etc.
# Run Each Model, one at a time

# setwd("U:\\Lasantha/Research/Rcode") 
# setwd("D:\\SFU/Research/Rcode/Rcode")
# setwd("C:\\Users/wpremara/Dropbox/Research/Rcode")

source("load.R")  # load required functions and packages

options(width=350)
set.seed(123)

##############################################################################
Data = get.data("Mille_Lacs_Walleye.csv") # get  Mille Lacs Walleye dataset 
str(Data)

##############################################################################
# following notation used for model id
#   t - time 
#   c - category
#  e.g : p(t*c) means capture probabilities  vary by time and category
#      : p(t)  means capture probability  vary by time but not category
#      : lambda(c) means category proportions vary by category
#      : p(.) means capture probabilities  not vary by time and category
#      : lambda(MLE) means category proportions fixed at MLE values  
#
# the integer after the p, lambda and theta indicate about the prior
# distributions for p, lambda and theta. if the integer is same for same
# parameters in two or more models, that indicate the prior distributios are same
# for those marametes
#
# Ex: The models "{ p(c*t)1, theta(t)1, lambda(c)1 } and 
# { p(c*t)1, theta(t)1, lambda(c)2 } are both unrestricted models. But they have 
# different prior distributions for lambda
#
#
# If the parameter values are fixed, then pass them through the offset vectors. 
# Required fixed values should be given as a logit value
# { eg. If capture probabilities at time 1 for male and female are 0.2 and 0.3
# and time 2 are 0.3 and 0.2 then  the capture offset is
# c( logit(0.2), logit(0.3), logit(0.3),logit(0.2) ) }


###############################################################################
######## visulizer functions are in the file " bayesian.visualizer.R" #########
###############################################################################

################
################ regular values to logit/log scale

## Following help to figure out the mean and sd for the  normal prior 
## distributions in logit/log form.  

#### enter the mean and sd for prior for probability in regular form ##########
visualizer.prob.to.logit.prob(priormean=0.2,priorsd=0.05)

#### enter the mean and sd for prior for population size in regular form ######
visualizer.N.to.log.N(priormean=205000,priorsd=40000)

#### enter the  the vector of parameters for Dirichlet distribution in  #######
#### regular form and vector of categories for prior distribution #############      
visualizer.lambda.to.logit.lambda(c(20,40),Data$category)

###############
###############  logit/log scale to regular values

#### enter the mean and sd for prior for probability in logit form ##########
visualizer.logit.to.regular.prob(logit_priormean=-2,logit_priorsd=1)

#### enter the mean and sd for prior for population size in log form ######
visualizer.log.N.to.N(log_priormean=12.25,log_priorsd=1)


##############################################################################
###############################################################################
###############################################################################

# Metropolis Hastings analysis
# User has to specify the following
#   n.updates : give the number of iterations
#   n.burn.in : give the number of iterations to be discarded as burn in
#   n.thin : e.g. if n.thin=50, then every 50th iteration will be saved
#   n.chains : number of chains
#   alpha : what size of credible intervals needed 
#           (eg : c(0.05) means 95%  credible intervals)
#   nstops : adapting the step size . acceptance rate are calculated for
#            each n.updates/nstops iterations. Adaptively change the sd of
#            proposal for the n.burn.in iterations
#   priors : prior distributions for beta parameters(as in MARK)
#   initial.proposal.sigma : initial sd of the proposal distributions
#   target.acc.rate : target acceptance rate of the proposal

###############################################################################
###############################################################################
# unrestricted model with prior for category proportion Dirichlet(2,4)

#model identification : unrestricted model
model.id <- paste("{ p(c*t)(0,1), theta(t)(0,1), lambda(c)(2,4) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)


######## priors - prior distributions for beta parameters(as in MARK) ########

priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(0,0,0,0 ) #vector of length ncol(design$p) 
priors$p.sigma <-c(1.78,1.78,1.78,1.78) #vector of length ncol(design$p)

# priors for category proportion(lambda) are normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda)- Dirichlet(2,4) in regular scale
priors$lambda.mu <- c(-0.69) #vector of length ncol(design$lambda) 
priors$lambda.sigma <-c(0.96) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(0, 0) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1.78,1.78) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log  scale
priors$N.mu <- c(12.25) 
priors$N.sigma <-c(5) 

# fit the bayesian model 
bayes_model_1 <- bayesian.fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
                                    captureOFFSET,thetaOFFSET,lambdaOFFSET,
                                    n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                    nstops=20,n.chains=3,alpha=c(0.05),
                                    priors,
                                    initial.proposal.sigma = 0.1,
                                    target.acc.rate =0.3)

bayesian.print.output(bayes_model_1)

###############################################################################
###############################################################################
#model identification: Category proportions are equal (lambda_m = lambda_f=0.5)
model.id = paste("{ p(c*t)(0,1), theta(t)(1)(0,1), lambda(0.5) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM = create.DM(c(1,2,3,4))    
thetaDM   = create.DM(c(1,2))       
lambdaDM  = matrix(, ncol=0,nrow=1) 

#give the offset vectors(vectors of zeros should be given since no restriction)
captureOFFSET = c(0,0,0,0) 
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(logit(0.5))


# priors - prior distributions for beta parameters(as in MARK)
priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(0, 0,0,0 ) #vector of length ncol(design$p) 
priors$p.sigma <-c(1.78,1.78,1.78,1.78) #vector of length ncol(design$p)

## since category proportions are fixed, no need to give the prior for lambda

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(0, 0) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1.78,1.78) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution in log  scale
priors$N.mu <- c(12.25) 
priors$N.sigma <-c(5) 


# fit the bayesian model 
bayes_model_2 <- bayesian.fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
                                    captureOFFSET,thetaOFFSET,lambdaOFFSET,
                                    n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                    nstops=20,n.chains=3,alpha=c(0.05),
                                    priors,
                                    initial.proposal.sigma = 0.1,
                                    target.acc.rate =0.3)

bayesian.print.output(bayes_model_2)


###############################################################################
###############################################################################
#model identification:  Capture Probabilities vary by category but not time. 
model.id = paste("{ p(c)(0,1), theta(t)(0.7,1)(0,0.3), lambda(c)(2,4)}")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM = create.DM(c(1,2,1,2)) 
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = c(0,0,0,0) 
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

# priors - prior distributions for beta parameters(as in MARK)
priors <- NULL 

# priors for capture historyare normal distribution in logit scale
priors$p.mu <- c(0, 0) #vector of length ncol(design$p) 
priors$p.sigma <-c (1.78,1.78) #vector of length ncol(design$p)

# priors for category proportion(lambda) are normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda)
priors$lambda.mu <- c(-0.69) #vector of length ncol(design$lambda) 
priors$lambda.sigma <-c(0.96) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(2, -2) #vector of length ncol(design$theta) 
priors$theta.sigma <- c(1,1) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log scale
priors$N.mu <- c(12.25) 
priors$N.sigma <-c(5) 


# fit the bayesian model 
bayes_model_3 <- bayesian.fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
                                    captureOFFSET,thetaOFFSET,lambdaOFFSET,
                                    n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                    nstops=20,n.chains=3,alpha=c(0.05),
                                    priors,
                                    initial.proposal.sigma = 0.1,
                                    target.acc.rate =0.3)

bayesian.print.output(bayes_model_3)

###############################################################################
##############################################################################
# unrestricted model With diffrent prior for lambda, Dirichlet(2,4)
#model identification : unrestricted model
model.id <- paste("{ p(c*t)(0,0.3), theta(t)(0.7,1)(0,0.3), lambda(c)(2,4) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)


######## priors - prior distributions for beta parameters(as in MARK) ########

priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(-2, -2,-2,-2 ) #vector of length ncol(design$p) 
priors$p.sigma <-c(1,1,1,1) #vector of length ncol(design$p)

# priors for category proportion(lambda) are normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda), Dirichlet(2,4) in regular scale
priors$lambda.mu <- c(-0.69) #vector of length ncol(design$lambda) 
priors$lambda.sigma <-c(0.96) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(2, -2) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1,1) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log  scale
priors$N.mu <- c(12.25) 
priors$N.sigma <-c(5) 

# fit the bayesian model 
bayes_model_4 <- bayesian.fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
                                    captureOFFSET,thetaOFFSET,lambdaOFFSET,
                                    n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                    nstops=20,n.chains=3,alpha=c(0.05),
                                    priors,
                                    initial.proposal.sigma = 0.1,
                                    target.acc.rate =0.3)

bayesian.print.output(bayes_model_4)

###############################################################################

# unrestricted model With diffrent prior for lambda, Dirichlet(3,4)
#model identification : unrestricted model
model.id <- paste("{ p(c*t)(0,0.3), theta(t)(0.7,1)(0,0.3), lambda(c)(3,4) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)


######## priors - prior distributions for beta parameters(as in MARK) ########

priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(-2, -2,-2,-2 ) #vector of length ncol(design$p) 
priors$p.sigma <-c(1,1,1,1) #vector of length ncol(design$p)

# priors for category proportion(lambda) are normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda), Dirichlet(3,4) in regular scale
priors$lambda.mu <- c(-0.29) #vector of length ncol(design$lambda) 
priors$lambda.sigma <-c(0.83) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(2, -2) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1,1) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log  scale
priors$N.mu <- c(12.25) 
priors$N.sigma <-c(5) 

# fit the bayesian model 
bayes_model_5 <- bayesian.fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
                                    captureOFFSET,thetaOFFSET,lambdaOFFSET,
                                    n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                    nstops=20,n.chains=3,alpha=c(0.05),
                                    priors,
                                    initial.proposal.sigma = 0.1,
                                    target.acc.rate =0.3)

bayesian.print.output(bayes_model_5)

###############################################################################
# unrestricted model With diffrent prior for lambda, Dirichlet(20,40)
#model identification : unrestricted model
model.id <- paste("{ p(c*t)(0,0.3), theta(t)(0.7,1)(0,0.3), lambda(c)(20,40) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)


######## priors - prior distributions for beta parameters(as in MARK) ########

priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(-2, -2,-2,-2 ) #vector of length ncol(design$p) 
priors$p.sigma <-c(1,1,1,1) #vector of length ncol(design$p)

# priors for category proportion(lambda) are normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda), Dirichlet(20,40) in regular scale
priors$lambda.mu <- c(-0.69) #vector of length ncol(design$lambda) 
priors$lambda.sigma <-c(0.27) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(2, -2) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1,1) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log  scale
priors$N.mu <- c(12.25) 
priors$N.sigma <-c(5) 

# fit the bayesian model 
bayes_model_6 <- bayesian.fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
                                    captureOFFSET,thetaOFFSET,lambdaOFFSET,
                                    n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                    nstops=20,n.chains=3,alpha=c(0.05),
                                    priors,
                                    initial.proposal.sigma = 0.1,
                                    target.acc.rate =0.3)

bayesian.print.output(bayes_model_6)

###############################################################################
######## Bayesian Model comparison table.######################################

# Model comparison table
# enter which method have to be used in DIC calculation 
# enter  "Spiegelhalter" or "Gelman" for two different methods of DIC calculation
#   method = "Spiegelhalter" to use the "pD" for DIC calculation
#   method = "Gelman" to use the "pv" for DIC calculation


#### Two different ways of DIC calculation ########
#
#  DEVIANCE=  D(theta)= -2 * log(likelihood) 
#
#   
#        "Dbar" is the posterior mean of the deviance,
#               (i.e. the average of D(theta) over the samples theta)
#        "D(theta_bar)" is the value of the Deviance evaluated at the average of
#                       the samples of theta
#        "pD" is  the effective number of parameters
#         pD= "posterior mean deviance - deviance of posterior means"
#
### Method 1. ####  
#       pD = Dbar - D(theta_bar)     # Spiegelhalter et al. (2002)
#           
### Method 2. ####               
#       pD = pv = 1/2  *  var(D(theta))   # Gelman et al (2004)
#
#  DIC = DIC = Dbar + pD   or  DIC = D(thata_bar) + 2 pD

## 


bayesian.model.comparison.table(method="Spiegelhalter", 
                                models = list(bayes_model_1, bayes_model_2,
                                              bayes_model_3, bayes_model_4,
                                              bayes_model_5, bayes_model_6)) 

bayesian.model.comparison.table(method="Gelman", 
                                models = list(bayes_model_1, bayes_model_2,
                                              bayes_model_3, bayes_model_4,
                                              bayes_model_5, bayes_model_6)) 

###############################################################################
###############################################################################
######################### END #################################################




