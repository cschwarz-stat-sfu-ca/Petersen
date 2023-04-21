###############################################################################
#########  Bayesian Analysis of the Mille Lacs Walleye dataset  ###############
####################  With individual covariate  ##############################
###############################################################################
####  Metropolis Hastings on each parameter in terms of using a symmetric #####
####  proposal distributions                                              #####            
###############################################################################
# We can create different models and  the results are stored in
# IC.BAYESIAN_model_1, IC.BAYESIAN_model_2...etc.
# Run Each Model, one at a time

#setwd("D:\\SFU/Research/Rcode/Rcode")
#setwd("U:\\Lasantha/Research/Rcode")

source("load.R")  # load required functions and packages

options(width=350)

#Data <- ic.get.data("ic.test.data.csv") # get  Mille Lacs Walleye dataset 
Data <- ic.get.data("ic.millelacs-2013.csv")
#Data <- ic.get.data("ic.test.data.csv")
str(Data)

##############################################################################
######## take a sample of 1000 rows from the millelacs-2013 data set #########
## This is to run the code and check 
df <- data.frame(Data[1:4])
sample_df <-df[sample(nrow(df), 400),]

Data$History <- as.character(sample_df$History)
Data$counts   <- sample_df$counts 
Data$length   <- sample_df$length
Data$lengthsq <- sample_df$lengthsq
str(Data)
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
######## priors - prior distributions for beta parameters(as in MARK) #########
###############################################################################
###############################################################################
######## visulizer functions are in the file " bayesian.visualizer.R" #########

################
################ regular values to logit/log scale

## Following help to figure out the mean and sd for the  normal prior 
## distributions in logit form. 

#### enter the mean and sd for prior for probability in regular form ##########
visualizer.prob.to.logit.prob(priormean=0.2,priorsd=0.04)

#### enter the mean and sd for prior for population size in regular form ######
visualizer.N.to.log.N(priormean=205000,priorsd=40000)

#### enter the  the vector of parameters for Dirichlet distribution in  #######
#### regular form and vector of categories for prior distribution #############      
visualizer.lambda.to.logit.lambda(c(20,40),Data$category)

###############
###############  logit scale to regular values

#### enter the mean and sd for prior for probability in logit form ##########
visualizer.logit.to.regular.prob(logit_priormean=0,logit_priorsd=1.78)

visualizer.logit.to.regular.prob(logit_priormean=0,logit_priorsd=10)

visualizer.logit.to.regular.prob(logit_priormean=3,logit_priorsd=.2)

###############################################################################
###############################################################################
###############################################################################
# Bayesian Model 1

# give the capture formula
captureformula <- ~length+ lengthsq+ category + time

# number of beta parameters related to capture formula
#  need add 1, because of the intercept
n.beta.param.cap <- length(colnames(attr(terms(captureformula), "factors"))) +1

# design Matrices
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = create.DM(c(1)) 

n.rows <- length(Data$History) # number of rows of data
#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(Data$category) *2) )
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

#model identification :All capture probabilities are equal
model.id = paste("{ p(",captureformula,"), theta(t), lambda(c) }",sep="")[2]


# priors - prior distributions for beta parameters(as in MARK)
priors <- NULL 

# priors for the parameters in capture formula in logit scale
# capture probabilities depend on the capture formula
# according to the capture formula we have to give the priors for the parameters for 
# Constant, length, lengthsq, category and time
# have to specify 5 prior distributions according to capture formula
## Give a uniform prior(in regular form) for the Constant( ie. N(0,1.78) in logit form)
## Give plat priors for length, lengthsq, category and time coefficients
##    ( Normal priors with mean 0  and large SD. ie. N( 0, 10) in logit form)

priors$beta.logit.p.mu <- c(0, 0, 0, 0, 0) #vector of length of  n.beta.param.cap
priors$beta.logit.p.sigma <-c(1.78,10,10,10,10) #vector of length of  n.beta.param.cap

# priors for category proportion(lambda) is normal distribution in logit scale
# use Dirichlet prior on the category proportion(lambda)
priors$lambda.mu <- c(2) #vector of length ncol(design$lambda) 
priors$lambda.sigma <-c(0.2) #vector of length ncol(design$lambda)

# priors for sub-sample proportion(theta) is normal distribution in logit scale
priors$theta.mu <- c(1, 0.9) #vector of length ncol(design$theta) 
priors$theta.sigma <-c(0.2,0.2) #vector of length ncol(design$theta)

# ################### to be deleted
# n.post.burnin=140
# n.burn.in =140
# n.thin=2
# nstops=10
# n.chains=2
# alpha=c(0.05)
# initial.proposal.sigma = 0.1
# target.acc.rate =0.3
# #################### to be deleted


# fit the bayesian model 
IC.BAYESIAN_model_1 <- IC.BAYESIAN.fit.model(model.id,Data,
                                             captureformula,thetaDM,lambdaDM,
                                             captureOFFSET,thetaOFFSET,lambdaOFFSET,
                                             n.post.burnin=200, n.burn.in =200,n.thin=2,
                                             nstops=10,n.chains=2,alpha=c(0.05),
                                             priors,
                                             initial.proposal.sigma = 0.1,
                                             target.acc.rate =0.3)
IC.BAYESIAN.print.output(IC.BAYESIAN_model_1)
IC.BAYESIAN.p.value(IC.BAYESIAN_model_1)
###############################################################################
###############################################################################











# 
# 
# ###############################################################################
# #model identification: Category proportions are equal (lambda_m = lambda_f=0.5)
# model.id = paste("{ p(c*t), theta(t), lambda(0.5) }")
# 
# # give  the required design matrices for capture recapture probabilities,
# # theta(sampling(sexing) fractions ) and lambda(Category proportion)
# captureDM = create.DM(c(1,2,3,4))    
# thetaDM   = create.DM(c(1,2))       
# lambdaDM  = matrix(, ncol=0,nrow=1) 
# 
# #give the offset vectors(vectors of zeros should be given since no restriction)
# captureOFFSET = c(0,0,0,0) 
# thetaOFFSET   = c(0,0)
# lambdaOFFSET  = c(logit(0.5))
# 
# 
# # priors - prior distributions for beta parameters(as in MARK)
# priors <- NULL 
# 
# # priors for capture history is normal distribution in logit scale
# priors$p.mu <- c(-2.5, -4.5,-4.8,-3.8 ) #vector of length ncol(design$p) 
# priors$p.sigma <-c(0.2,0.3,0.2,0.2) #vector of length ncol(design$p)
# 
# ## since category proportions are fixed, no need to give the prior for lambda
# 
# # priors for sub-sample proportion(theta) is normal distribution in logit scale
# priors$theta.mu <- c(5, -2.4) #vector of length ncol(design$theta) 
# priors$theta.sigma <-c(0.5,0.3) #vector of length ncol(design$theta)
# 
# # Prior distribution for the population size(N)  in log normal scale
# priors$N.mu <- c(12.25) 
# priors$N.sigma <-c(0.4) 
# 
# 
# # fit the bayesian model 
# bayes_model_2 <- bayesian.fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
#                                     captureOFFSET,thetaOFFSET,lambdaOFFSET,
#                                     n.post.burnin=100000, n.burn.in =60000,n.thin=50,
#                                     nstops=20,n.chains=3,alpha=c(0.05),
#                                     priors,
#                                     initial.proposal.sigma = 0.1,
#                                     target.acc.rate =0.3)
# 
# bayesian.print.output(bayes_model_2)
# 
# 
# ###############################################################################
# ###############################################################################
# #model identification:  Capture Probabilities vary by category but not time. 
# model.id = paste("{ p(c), theta(t), lambda(c) }")
# 
# # give  the required design matrices for capture recapture probabilities,
# # theta(sampling(sexing) fractions ) and lambda(Category proportion)
# captureDM = create.DM(c(1,2,1,2)) 
# thetaDM   = create.DM(c(1,2)) 
# lambdaDM  = create.DM(c(1)) 
# 
# #give the offset vectors(vectors of zero's should be given since no restriction)
# captureOFFSET = c(0,0,0,0) 
# thetaOFFSET   = c(0,0)
# lambdaOFFSET  = c(0)
# 
# # priors - prior distributions for beta parameters(as in MARK)
# priors <- NULL 
# 
# # priors for capture history is normal distribution in logit scale
# priors$p.mu <- c(-3, -4) #vector of length ncol(design$p) 
# priors$p.sigma <-c (0.2,0.2) #vector of length ncol(design$p)
# 
# # priors for category proportion(lambda) is normal distribution in logit scale
# # use Dirichlet prior on the category proportion(lambda)
# priors$lambda.mu <- c(-0.7) #vector of length ncol(design$lambda) 
# priors$lambda.sigma <-c(0.3) #vector of length ncol(design$lambda)
# 
# # priors for sub-sample proportion(theta) is normal distribution in logit scale
# priors$theta.mu <- c(5, -2.4) #vector of length ncol(design$theta) 
# priors$theta.sigma <- c(0.5,0.3) #vector of length ncol(design$theta)
# 
# # Prior distribution for the population size(N)  in log normal scale
# priors$N.mu <- c(12.25) 
# priors$N.sigma <-c(0.4) 
# 
# 
# # fit the bayesian model 
# bayes_model_3 <- bayesian.fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
#                                     captureOFFSET,thetaOFFSET,lambdaOFFSET,
#                                     n.post.burnin=100000, n.burn.in =60000,n.thin=50,
#                                     nstops=20,n.chains=3,alpha=c(0.05),
#                                     priors,
#                                     initial.proposal.sigma = 0.1,
#                                     target.acc.rate =0.3)
# 
# bayesian.print.output(bayes_model_3)
# 
# ###############################################################################
# ##############################################################################
# 
# 
# 
# ###############################################################################
# ######## Bayesian Model comparison table.######################################
# 
# # Model comparison table
# # enter which method have to be used in DIC calculation 
# # enter  "Spiegelhalter" or "Gelman" for two different methods of DIC calculation
# #   method = "Spiegelhalter" to use the "pD" for DIC calculation
# #   method = "Gelman" to use the "pv" for DIC calculation
# 
# 
# #### Two different ways of DIC calculation ########
# #
# #  DEVIANCE=  D(theta)= -2 * log(likelihood) 
# #
# #   
# #        "Dbar" is the posterior mean of the deviance,
# #               (i.e. the average of D(theta) over the samples theta)
# #        "D(theta_bar)" is the value of the Deviance evaluated at the average of
# #                       the samples of theta
# #        "pD" is  the effective number of parameters
# #         pD= "posterior mean deviance - deviance of posterior means"
# #
# ### Method 1. ####  
# #       pD = Dbar - D(theta_bar)     # Spiegelhalter et al. (2002)
# #           
# ### Method 2. ####               
# #       pD = pv = 1/2  *  var(D(theta))   # Gelman et al (2004)
# #
# #  DIC = DIC = Dbar + pD   or  DIC = D(thata_bar) + 2 pD
# 
# ## 
# 
# bayesian.model.comparison.table(method="Spiegelhalter", 
#                                 models = list(bayes_model_1,bayes_model_2,
#                                               bayes_model_3)) 
# 
# ###############################################################################
###############################################################################
######################### END #################################################



