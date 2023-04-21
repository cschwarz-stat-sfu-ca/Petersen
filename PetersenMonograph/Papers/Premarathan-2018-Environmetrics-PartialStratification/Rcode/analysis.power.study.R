###############################################################################
############################################################################### 
##########  power of the study and verify the power using simulation ##########
##########             Mille Lacs Walleye dataset                    ##########
###############################################################################

source("load.R")  # load required functions and packages

options(width=350)

# generate some data where the probabilities vary between sexes/sample times and 
# fit a model where the # probabilities are equal then use the power of the 
# study to detect such differences Verify the power using simulation study

# give the  categories in the population
categories = c("M","F")

# give the different  capture probabilities for p1M, p1F, p2M and p2F

# delta is the difference between the capture probability of male and female 
# at each sample time
delta = 0.02

pm1 = 0.08  # capture probability of male at time 1
pm2 = 0.04   # capture probability of male at time 2

cap.prob = c( pm1, pm1+delta, pm2, pm2+delta)

# give the category proportions( total should add up to 1)
lambda= c(0.6, 0.4)

# give the sub-sample proportions for the time 1 and 2
theta = c(0.8, 0.5)

given.param = NULL
given.param$p = matrix(cap.prob, ncol=2, byrow = FALSE)
given.param$lambda = lambda
given.param$theta  = theta

# total numbers individuals capture for the study 
sample.size = 2000


#############  alternative (unrestricted)  model ##############################
#model identification : alternative model (unrestricted model)
unrest.model.id = paste("{ p(c*t), theta(t), lambda(c) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
unrest.captureDM = create.DM(c(1,2,3,4)) 
unrest.thetaDM   = create.DM(c(1,2)) 
unrest.lambdaDM  = create.DM(c(1))  

#give the offset vectors(vectors of zeros should be given since no restriction)
unrest.captureOFFSET = c(0,0,0,0) 
unrest.thetaOFFSET   = c(0,0)
unrest.lambdaOFFSET  = c(0)

############ Null (restricted) model ##########################################
#model identification : nullmodel
#                       (capture probabilities not vary by category but time )
rest.model.id = paste("{ p(t), theta(t), lambda(c) }")

# Design matrices restricted model
rest.captureDM = create.DM(c(1,1,2,2))  # capture probabilities are equal
rest.thetaDM   = create.DM(c(1,2))
rest.lambdaDM  = create.DM(c(1)) 

#offset vectors for restricted model
rest.captureOFFSET = c(0,0,0,0) 
rest.thetaOFFSET   = c(0,0)
rest.lambdaOFFSET  = c(0)

########### power analysis plots #################################

# "n.simulation"  is the number of simulations to verify the power of the study.
# "alpha" is the type I error probability
power.study.result = power.study.simulation(given.param,categories,
                                            unrest.model.id,unrest.captureDM,unrest.thetaDM,
                                            unrest.lambdaDM,unrest.captureOFFSET,
                                            unrest.thetaOFFSET,unrest.lambdaOFFSET,
                                            rest.model.id,rest.captureDM,rest.thetaDM,
                                            rest.lambdaDM,rest.captureOFFSET,
                                            rest.thetaOFFSET,rest.lambdaOFFSET,
                                            sample.size,n.simulation=1000, alpha=0.05)

power.study.result


############# verify the power using Devineau method and simulations ##########
## considered sample size  2000
## alpha = 0.05
power.data_1  = read.csv("power_analysis.csv", header = TRUE, as.is=TRUE)
power.analysis.plot_1 = ggplot(power.data_1, aes(DeltaP,PowerDevineau)) + 
                          geom_line(aes(y = PowerDevineau,colour = "Devineau"),size=1) + 
                          geom_point(aes(y = PowerSimulation,shape = "simulation"),
                                    size = 3) +
                          theme(legend.title=element_blank())+
                          xlab("Delta_p") + ylab("Power")+
                          ggtitle(paste("Verify Power using Devineau Method and simulations",
                                  "\n sample size = 2000, alpha = ", 0.05,
                                  "\n Model under H0 : { p(t), theta(t), lambda(c) }",
                                  "\n Model under Ha : { p(c*t), theta(t), lambda(c) }",
                                   sep=" "))+
                          theme(axis.text=element_text(size=10,face="bold"),
                                  axis.title=element_text(size=14,face="bold"),
                                  plot.title = element_text(size = 14,face="bold"))
power.analysis.plot_1

#### power analysis for different samples sizes ####
## considered sample sizes  n = 1000, 2000, 3000, 4000
## alpha = 0.05
power.data_2 = read.csv("power_analysis_for_different_samples.csv", header = TRUE, as.is=TRUE)
power.analysis.plot_2 = ggplot(power.data_2, aes(DeltaP)) + 
                            geom_line(aes(y = power1000,linetype="n=1000"),size=1) + 
                            geom_line(aes(y = power2000,linetype="n=2000"),size=1) + 
                            geom_line(aes(y = power3000,linetype="n=3000"),size=1) + 
                            geom_line(aes(y = power4000,linetype="n=4000"),size=1) + 
                                  labs(linetype='Sample Size') +
                            xlab("Delta_p") + ylab("Power")+
                            ggtitle(paste("Power Analysis for Different Sample Sizes at alpha = ",
                                  0.05, "\n Model under H0 : { p(t), theta(t), lambda(c) }",
                                  "\n Model under Ha : { p(c*t), theta(t), lambda(c) }",
                                  sep=" "))+
                            theme(axis.text=element_text(size=10,face="bold"),
                                  axis.title=element_text(size=14,face="bold"),
                                  plot.title = element_text(size = 14,face="bold"))
power.analysis.plot_2


###############################################################################
##power of the study under  H0 model : ( p(c*t), theta(t), lambda(.) ) and 
#                           Ha model : ( p(c*t), theta(t), lambda(c) )
###############################################################################

# generate some data where the probabilities vary between sexes and sample times 
# detect such differences in lambda(.) model
# (null) H0 model : ( p(c*t), theta(t), lambda(.) )
# (alternative) Ha model : ( p(c*t), theta(t), lambda(c) )
# Verify the power using simualation study
# lambda_M in generated data is 0.4
# consider Delda_lambda = 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 
#Then the corresponding values of lambda under H0 are 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7


# give the  categories in the population
categories = c("M","F")

# give the different  capture probabilities for p1M, p1F, p2M and p2F
cap.prob = c( 0.05, 0.08, 0.1, 0.12)

# give the category proportions( total should add up to 1)
lambda= c(0.4, 0.6)

# give the subsample proportions for the time 1 and 2
theta = c(0.8, 0.5)

given.param = NULL
given.param$p = matrix(cap.prob, ncol=2, byrow = FALSE)
given.param$lambda = lambda
given.param$theta  = theta

# total numbers individuals capture for the study 
sample.size = 2000



######  alternative (unrestricted) model #####################################
#model identification : alternative model
unrest.model.id = paste("{ p(c*t), theta(t), lambda(c) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
unrest.captureDM = create.DM(c(1,2,3,4)) 
unrest.thetaDM   = create.DM(c(1,2)) 
unrest.lambdaDM  = create.DM(c(1))  

#give the offset vectors(vectors of zeros should be given since no restriction)
unrest.captureOFFSET = c(0,0,0,0) 
unrest.thetaOFFSET   = c(0,0)
unrest.lambdaOFFSET  = c(0)

##### null (restricted) model #################################################
#model identification : null model
rest.model.id = paste("{ p(c*t), theta(t), lambda(0.4+Delta) }")

# Design matrices restricted model
rest.captureDM = create.DM(c(1,2,3,4))  
rest.thetaDM   = create.DM(c(1,2))
rest.lambdaDM  = matrix(, ncol=0,nrow=1)

#offset vectors for restricted model
rest.captureOFFSET = c(0,0,0,0) 
rest.thetaOFFSET   = c(0,0)
rest.lambdaOFFSET  = c(logit(0.45))


########### power analysis plots #################################
# "alpha" is the type I error probability
power.study.result = power.study.simulation(given.param,categories,
                                            unrest.model.id,unrest.captureDM,unrest.thetaDM,
                                            unrest.lambdaDM,unrest.captureOFFSET,
                                            unrest.thetaOFFSET,unrest.lambdaOFFSET,
                                            rest.model.id,rest.captureDM,rest.thetaDM,
                                            rest.lambdaDM,rest.captureOFFSET,
                                            rest.thetaOFFSET,rest.lambdaOFFSET,
                                            sample.size,n.simulation=1000, alpha=0.05) 
power.study.result


# verify the power using Devineau method and simulations ####
## considered sample size  2000,  alpha = 0.05
power.data_2  = read.csv("power_analysis_for_lambda.csv", header = TRUE, as.is=TRUE)
power.analysis.plot_3 = ggplot(power.data_2, aes(Delta_lambda)) + 
                        geom_line(aes(y = power_n_2000_Devineau,
                                colour = "Devineau"),size=1) + 
                        geom_point(aes(y = power_n_2000_simulation,
                                shape = "simulation"), size = 3) +
                        theme(legend.title=element_blank())+
                        xlab("Delta") + ylab("Power")+
                        ggtitle(paste("Verify Power using Devineau Method and simulations",
                              "\n sample size = 2000, alpha = ", 0.05,
                              "\n Model under H0 : { p(c*t), theta(t), lambda(0.4+Delta) } ",
                              "\n Model under Ha : { p(c*t), theta(t), lambda(c) }",
                              sep=" "))+
                        theme(axis.text=element_text(size=10,face="bold"),
                              axis.title=element_text(size=14,face="bold"),
                              plot.title = element_text(size = 14,face="bold"))
power.analysis.plot_3


###############################################################################
############### End of power of the study #####################################
###############################################################################