############################################################################### 
#####                       Bayesian Analysis of                         ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################
####  Metropolis Hastings on each parameter in terms of using a symmetric #####
####  proposal distributions                                              #####            
###############################################################################
# We can create different models and  the results are stored in
# BAYESIAN_PSkSCR_model_1, BAYESIAN_PSkSCR_model_2, ...etc..
# Run Each Model, one at a time

#############  path below should points to the Functions directory  ###########
# setwd("U:\\Lasantha/Research/Rcode")
# setwd("C:\\Lasantha/SFU/Research/Rcode")
# setwd("D:\\SFU/Research/Rcode")

source("load.R")  # load required functions and packages
set.seed(123)


options(width=350)

##############################################################################
# get generated data to analysis
data <- PSkSCR.get.data("PSkSCR_Test_data_1.csv")
str(data)

##############################################################################

# following notation are used for model identification

### parameters ######################
## p - capture probabilities
## theta - sub-sample proportions
## lambda - category proportions
## nu - loss on capture probabilities
#####################################

#   t - time 
#   c - category
#  e.g : p(t*c) means capture probabilities  vary by time and category
#      : p(t)  means capture probability  vary by time but not category
#      : lambda(c) means category proportions vary by category
#      : p(.) means capture probabilities  not vary by time and category
#      : lambda(MLE) means category proportions fixed at MLE values  

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
visualizer.prob.to.logit.prob(priormean=0.01  ,priorsd=0.034)

#### enter the mean and sd for prior for population size in regular form ######
visualizer.N.to.log.N(priormean=381000,priorsd=50000)

#### enter the  the vector of parameters for Dirichlet distribution in  #######
#### regular form and vector of categories for prior distribution #############      
visualizer.lambda.to.logit.lambda(c(3,1),data$category)

###############
###############  logit/log scale to regular values

#### enter the mean and sd for prior for probability in logit form ##########
visualizer.logit.to.regular.prob(logit_priormean=-4,logit_priorsd=1)
visualizer.logit.to.regular.prob(logit_priormean=1.1,logit_priorsd=1.4)

#### enter the mean and sd for prior for population size in log form ######
visualizer.log.N.to.N(log_priormean=12.5,log_priorsd=5)


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


###############################################################################
###############################################################################
#### Model - 1 ####

#model identification : unrestricted model
model.id <- paste("{ p(c*t)(0,1), theta(t)(0,1), nu(t)(0,1), lambda(c)(3,1) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4,5,6)) 
thetaDM   <- create.DM(c(1,2,3)) 
lambdaDM  <- create.DM(c(1)) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(0)
p_lossOFFSET  <- c(0,0,0)


######## priors - prior distributions for beta parameters(as in MARK) ########


priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(0,0,0,0,0,0) #vector of length ncol(design$p) 
priors$p.sigma <-c(1.78, 1.78, 1.78, 1.78, 1.78, 1.78) #vector of length ncol(design$p)

# priors for category proportion(lambda) are normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda)
priors$lambda.mu <- c(1.1) #vector of length ncol(design$lambda) 
priors$lambda.sigma <- c(1.4) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(0,0,0) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1.78, 1.78, 1.78) #vector of length ncol(design$theta)

# priors for loss on capture(nu) are normal distribution in logit scale
priors$p_loss.mu <- c(0,0,0) #vector of length ncol(design$theta) 
priors$p_loss.sigma <-c(1.78, 1.78, 1.78) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log  scale
priors$N.mu <- c(12.8) 
priors$N.sigma <- c(5) 



# fit the bayesian model 
BAYESIAN_PSkSCR_model_1 <- PSkSCR.BAYESIAN.fit.model(model.id,data,
                                                     captureDM,thetaDM,lambdaDM,p_lossDM,
                                                     captureOFFSET,thetaOFFSET,
                                                     lambdaOFFSET,p_lossOFFSET,
                                                     n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                                     nstops=40,n.chains=3,alpha=c(0.05),
                                                     priors,
                                                     initial.proposal.sigma = 0.1,
                                                     target.acc.rate =0.4)
# print output
PSkSCR.BAYESIAN.print.output(BAYESIAN_PSkSCR_model_1)

# Bayesian p-value scatter plots using the Discrepancy functions 
# (a) deviance  (b) Freeman-Tukey (FT) statistic
PSkSCR_bayesian_predictive_plots_model_1 <- PSkSCR.BAYESIAN.p.value(BAYESIAN_PSkSCR_model_1)
PSkSCR_bayesian_predictive_plots_model_1

## Run a model one at a time and save the work space
save.image("U:\\Lasantha/Research/Rcode/output/Environment/PSkSCR_BAYESIAN_Model_1.RData")

###############################################################################
###############################################################################
#### Model - 2 ####
#model identification : unrestricted model
model.id <- paste("{ p(c*t)(0,1), theta(t)(0,1), nu(t)(0,1), lambda(c)(30,10) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4,5,6)) 
thetaDM   <- create.DM(c(1,2,3)) 
lambdaDM  <- create.DM(c(1)) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(0)
p_lossOFFSET  <- c(0,0,0)


######## priors - prior distributions for beta parameters(as in MARK) ########

priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(0,0,0,0,0,0) #vector of length ncol(design$p) 
priors$p.sigma <-c(1.78, 1.78, 1.78, 1.78, 1.78, 1.78) #vector of length ncol(design$p)

# priors for category proportion(lambda) are normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda)
priors$lambda.mu <- c(1.1) #vector of length ncol(design$lambda) 
priors$lambda.sigma <- c(0.4) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(0,0,0) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1.78, 1.78, 1.78) #vector of length ncol(design$theta)

# priors for loss on capture(nu) are normal distribution in logit scale
priors$p_loss.mu <- c(0,0,0) #vector of length ncol(design$theta) 
priors$p_loss.sigma <-c(1.78, 1.78, 1.78) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log  scale
priors$N.mu <- c(12.8) 
priors$N.sigma <- c(5) 



# fit the bayesian model 
BAYESIAN_PSkSCR_model_2 <- PSkSCR.BAYESIAN.fit.model(model.id,data,
                                                     captureDM,thetaDM,lambdaDM,p_lossDM,
                                                     captureOFFSET,thetaOFFSET,
                                                     lambdaOFFSET,p_lossOFFSET,
                                                     n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                                     nstops=40,n.chains=3,alpha=c(0.05),
                                                     priors,
                                                     initial.proposal.sigma = 0.1,
                                                     target.acc.rate =0.4)
# print output
PSkSCR.BAYESIAN.print.output(BAYESIAN_PSkSCR_model_2)

# Bayesian p-value scatter plots using the Discrepancy functions 
# (a) deviance  (b) Freeman-Tukey (FT) statistic
PSkSCR_bayesian_predictive_plots_model_2 <- PSkSCR.BAYESIAN.p.value(BAYESIAN_PSkSCR_model_2)
PSkSCR_bayesian_predictive_plots_model_2


## Run a model one at a time and save the work space
save.image("U:\\Lasantha/Research/Rcode/output/Environment/PSkSCR_BAYESIAN_Model_2.RData")


###############################################################################
###############################################################################
#### Model - 3 ####

#model identification : unrestricted model
model.id <- paste("{ p(c*t)(0,1), theta(t)(0,1), nu(t)(0,1), lambda(0.5) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4,5,6)) 
thetaDM   <- create.DM(c(1,2,3)) 
lambdaDM  <- matrix(, ncol=0,nrow=1) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(logit(0.5))
p_lossOFFSET  <- c(0,0,0)


######## priors - prior distributions for beta parameters(as in MARK) ########


priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(0,0,0,0,0,0) #vector of length ncol(design$p) 
priors$p.sigma <-c(1.78, 1.78, 1.78, 1.78, 1.78, 1.78) #vector of length ncol(design$p)

## since category proportions are fixed, no need to give the prior for lambda

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(0,0,0) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1.78, 1.78, 1.78) #vector of length ncol(design$theta)

# priors for loss on capture(nu) are normal distribution in logit scale
priors$p_loss.mu <- c(0,0,0) #vector of length ncol(design$theta) 
priors$p_loss.sigma <-c(1.78, 1.78, 1.78) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log  scale
priors$N.mu <- c(12.8) 
priors$N.sigma <- c(5) 



# fit the bayesian model 
BAYESIAN_PSkSCR_model_3 <- PSkSCR.BAYESIAN.fit.model(model.id,data,
                                                     captureDM,thetaDM,lambdaDM,p_lossDM,
                                                     captureOFFSET,thetaOFFSET,
                                                     lambdaOFFSET,p_lossOFFSET,
                                                     n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                                     nstops=40,n.chains=3,alpha=c(0.05),
                                                     priors,
                                                     initial.proposal.sigma = 0.1,
                                                     target.acc.rate =0.4)
# print output
PSkSCR.BAYESIAN.print.output(BAYESIAN_PSkSCR_model_3)

# Bayesian p-value scatter plots using the Discrepancy functions 
# (a) deviance  (b) Freeman-Tukey (FT) statistic
PSkSCR_bayesian_predictive_plots_model_3 <- PSkSCR.BAYESIAN.p.value(BAYESIAN_PSkSCR_model_3)
PSkSCR_bayesian_predictive_plots_model_3


## Run a model one at a time and save the work space
save.image("U:\\Lasantha/Research/Rcode/output/Environment/PSkSCR_BAYESIAN_Model_3.RData")


###############################################################################
###############################################################################
#### Model - 4 ####

#model identification : unrestricted model
model.id <- paste("{ p(c*t)(0,0.1), theta(t)(0,1), nu(t)(0,0.1), lambda(c)(3,1) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4,5,6)) 
thetaDM   <- create.DM(c(1,2,3)) 
lambdaDM  <- create.DM(c(1)) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(0)
p_lossOFFSET  <- c(0,0,0)


######## priors - prior distributions for beta parameters(as in MARK) ########


priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(-4,-4,-4,-4,-4,-4) #vector of length ncol(design$p) 
priors$p.sigma <-c(1,1,1,1,1,1) #vector of length ncol(design$p)

# priors for category proportion(lambda) are normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda)
priors$lambda.mu <- c(1.1) #vector of length ncol(design$lambda) 
priors$lambda.sigma <- c(1.4) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(0,0,0) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1.78, 1.78, 1.78) #vector of length ncol(design$theta)

# priors for loss on capture(nu) are normal distribution in logit scale
priors$p_loss.mu <- c(-4,-4,-4) #vector of length ncol(design$theta) 
priors$p_loss.sigma <-c(1,1,1) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log  scale
priors$N.mu <- c(12.8) 
priors$N.sigma <- c(5) 



# fit the bayesian model 
BAYESIAN_PSkSCR_model_4 <- PSkSCR.BAYESIAN.fit.model(model.id,data,
                                                     captureDM,thetaDM,lambdaDM,p_lossDM,
                                                     captureOFFSET,thetaOFFSET,
                                                     lambdaOFFSET,p_lossOFFSET,
                                                     n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                                     nstops=40,n.chains=3,alpha=c(0.05),
                                                     priors,
                                                     initial.proposal.sigma = 0.1,
                                                     target.acc.rate =0.4)
# print output
PSkSCR.BAYESIAN.print.output(BAYESIAN_PSkSCR_model_4)

# Bayesian p-value scatter plots using the Discrepancy functions 
# (a) deviance  (b) Freeman-Tukey (FT) statistic
PSkSCR_bayesian_predictive_plots_model_4 <- PSkSCR.BAYESIAN.p.value(BAYESIAN_PSkSCR_model_4)
PSkSCR_bayesian_predictive_plots_model_4


## Run a model one at a time and save the work space
save.image("U:\\Lasantha/Research/Rcode/output/Environment/PSkSCR_BAYESIAN_Model_4.RData")


###############################################################################
###############################################################################
#### Model - 5 ####

#model identification : unrestricted model
model.id <- paste("{ p(c*t)(0,0.1), theta(t)(0,1), nu(t)(0,0.1), lambda(c)(30,10) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4,5,6)) 
thetaDM   <- create.DM(c(1,2,3)) 
lambdaDM  <- create.DM(c(1)) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(0)
p_lossOFFSET  <- c(0,0,0)


######## priors - prior distributions for beta parameters(as in MARK) ########

priors <- NULL 

# priors for capture history are normal distribution in logit scale
priors$p.mu <- c(-4,-4,-4,-4,-4,-4) #vector of length ncol(design$p) 
priors$p.sigma <-c(1,1,1,1,1,1) #vector of length ncol(design$p)

# priors for category proportion(lambda) are normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda)
priors$lambda.mu <- c(1.1) #vector of length ncol(design$lambda) 
priors$lambda.sigma <- c(0.4) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) are normal distribution in logit scale
priors$theta.mu <- c(0,0,0) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(1.78, 1.78, 1.78) #vector of length ncol(design$theta)

# priors for loss on capture(nu) are normal distribution in logit scale
priors$p_loss.mu <- c(-4,-4,-4) #vector of length ncol(design$theta) 
priors$p_loss.sigma <-c(1,1,1) #vector of length ncol(design$theta)

# Prior distribution for the population size(N) is normal distribution  in log  scale
priors$N.mu <- c(12.8) 
priors$N.sigma <- c(5) 



# fit the bayesian model 
BAYESIAN_PSkSCR_model_5 <- PSkSCR.BAYESIAN.fit.model(model.id,data,
                                                     captureDM,thetaDM,lambdaDM,p_lossDM,
                                                     captureOFFSET,thetaOFFSET,
                                                     lambdaOFFSET,p_lossOFFSET,
                                                     n.post.burnin=100000, n.burn.in =60000,n.thin=50,
                                                     nstops=40,n.chains=3,alpha=c(0.05),
                                                     priors,
                                                     initial.proposal.sigma = 0.1,
                                                     target.acc.rate =0.4)
# print output
PSkSCR.BAYESIAN.print.output(BAYESIAN_PSkSCR_model_5)

# Bayesian p-value scatter plots using the Discrepancy functions 
# (a) deviance  (b) Freeman-Tukey (FT) statistic
PSkSCR_bayesian_predictive_plots_model_5 <- PSkSCR.BAYESIAN.p.value(BAYESIAN_PSkSCR_model_5)
PSkSCR_bayesian_predictive_plots_model_5


## Run a model one at a time and save the work space
save.image("U:\\Lasantha/Research/Rcode/output/Environment/PSkSCR_BAYESIAN_Model_5.RData")


###############################################################################
###############################################################################

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


PSkSCR.BAYESIAN.model.comparison.table(method="Spiegelhalter", 
                                       models = list(BAYESIAN_PSkSCR_model_1,
                                                     BAYESIAN_PSkSCR_model_2,
                                                     BAYESIAN_PSkSCR_model_3,
                                                     BAYESIAN_PSkSCR_model_4,
                                                     BAYESIAN_PSkSCR_model_5)) 

PSkSCR.BAYESIAN.model.comparison.table(method="Gelman", 
                                       models = list(BAYESIAN_PSkSCR_model_1,
                                                     BAYESIAN_PSkSCR_model_2,
                                                     BAYESIAN_PSkSCR_model_3,
                                                     BAYESIAN_PSkSCR_model_4,
                                                     BAYESIAN_PSkSCR_model_5)) 
###############################################################################
###############################################################################
######################### END #################################################

save.image("U:\\Lasantha/Research/Rcode/output/Environment/PSkSCR_BAYESIAN_Analysis.RData")
