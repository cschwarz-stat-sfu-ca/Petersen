###############################################################################
############## Analysis of the Mille Lacs Walleye dataset #####################
###############################################################################
# We can create different models and  the results are stored in
# MLE_model_1, MLE_model_2,...etc.
# Run Each Model, one at a time

#setwd("D:\\SFU/Research/Rcode/Rcode")
#setwd("U:\\Lasantha/Research/Rcode")
#setwd("c:\\Lasantha/SFU/Research/Rcode")

source("load.R")  # load required functions and packages
set.seed(123)

options(width=350)

# following notation used for model id
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
# if capture probabilities are not fixed then the captureoffset is
# c(0,0,0,0)
###############################################################################

Data = get.data("Mille_Lacs_Walleye.csv") # get  Mille Lacs Walleye dataset
str(Data)

###############################################################################
###############################################################################
#model identification : unrestricted model
model.id = paste("{ p(c*t), theta(t), lambda(c) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM = create.DM(c(1,2,3,4))
thetaDM   = create.DM(c(1,2))
lambdaDM  = create.DM(c(1))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = c(0,0,0,0)
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

MLE.model_1 =  fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
                         thetaOFFSET,lambdaOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
print.output(MLE.model_1)


MLE.model_1.IC =  ic.fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
                         thetaOFFSET,lambdaOFFSET)
print.output(MLE.model_1.IC)

#parametric bootstrap goodness of fit
bootstrap.gof.model_1 = parametric.bootstrap(MLE.model_1,n.boots=1000)
bootstrap.gof.model_1

###############################################################################
###############################################################################
#model identification: Category proportions are equal (lambda_m = lambda_f=0.5)
model.id = paste("{ p(c*t), theta(t), lambda(0.5) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM = create.DM(c(1,2,3,4))
thetaDM   = create.DM(c(1,2))
lambdaDM  = matrix(0, ncol=0,nrow=1)

#give the offset vectors(vectors of zeros should be given since no restriction)
captureOFFSET = c(0,0,0,0)
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(logit(0.5))

MLE.model_2 =  fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
                         thetaOFFSET,lambdaOFFSET)
# print results (Model information,  MLEs, SE, residual plot, etc..)
print.output(MLE.model_2)

#parametric bootstrap goodness of fit
bootstrap.gof.model_2 = parametric.bootstrap(MLE.model_2,n.boots=1000)
bootstrap.gof.model_2

###############################################################################
###############################################################################
#model identification:  Capture Probabilities vary by category but not time.
model.id = paste("{ p(c), theta(t), lambda(c) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM = create.DM(c(1,2,1,2))
thetaDM   = create.DM(c(1,2))
lambdaDM  = create.DM(c(1))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = c(0,0,0,0)
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

MLE.model_3=  fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
                        thetaOFFSET,lambdaOFFSET)
# print results (Model information,  MLEs, SE,  residual plot, etc..)
print.output(MLE.model_3)

#parametric bootstrap goodness of fit
bootstrap.gof.model_3 = parametric.bootstrap(MLE.model_3,n.boots=1000)
bootstrap.gof.model_3

###############################################################################
###############################################################################
#model identification :All capture probabilities are equal
model.id = paste("{ p(.), theta(t), lambda(c) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM = create.DM(c(1,1,1,1))
thetaDM   = create.DM(c(1,2))
lambdaDM  = create.DM(c(1))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = c(0,0,0,0)
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

MLE.model_4 = fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
                         thetaOFFSET,lambdaOFFSET)
# print results (Model information,  MLEs, SE, residual plot, etc..)
print.output(MLE.model_4 )

#parametric bootstrap goodness of fit
bootstrap.gof.model_4 = parametric.bootstrap(MLE.model_4,n.boots=1000)
bootstrap.gof.model_4

###############################################################################
###############################################################################
#model identification : variation in capture probability by category and by time.
# Additive model with no interaction
model.id = paste("{ p(c+t), theta(t), lambda(c) }")

# give  the required design matrices
# Design matrix for capture recapture probabilities
captureDM = matrix(c(1,1,1,
                     1,0,1,
                     1,1,0,
                     1,0,0), nrow=4,byrow = TRUE)
thetaDM   = create.DM(c(1,2)) #Design matrix for theta(sampling(sexing) fractions)
lambdaDM  = create.DM(c(1))   #Design matrix for lambda(Category proportion)

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = c(0,0,0,0)
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

MLE.model_5 = fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
                         thetaOFFSET,lambdaOFFSET)
# print results (Model information,  MLEs, SE, residual plot, etc..)
print.output(MLE.model_5)

#parametric bootstrap goodness of fit
bootstrap.gof.model_5 = parametric.bootstrap(MLE.model_5,n.boots=1000)
bootstrap.gof.model_5

###############################################################################
###############################################################################
# Simple Petersen model
#model identification :All capture probabilities are equal
model.id = paste("{ p(t), theta(t), lambda(c) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM = create.DM(c(1,1,2,2))
thetaDM   = create.DM(c(1,2))
lambdaDM  = create.DM(c(1))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = c(0,0,0,0)
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

model_6 =  fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
                     thetaOFFSET,lambdaOFFSET)
# print results (Model information,  MLEs, SE, residual plot, etc..)
print.output(model_6 )

###############################################################################
######## Model comparison table.##############################################

# Model comparison table
model.comparison.table(models = list(MLE.model_1,MLE.model_2,MLE.model_3,
                                            MLE.model_4, MLE.model_5,model_6))

# "model.comparison.table" will produce a table which is ordered according to
# the AICc values for all the models and also maximum of 4 residual plots in
# a 2x2 grid for the models corresponding to minimum AICc values.
# If number of models considered is less than 4, it will produce all the
# residual plots for the considered models in 2x2 grid


############### End of Model comparison table.  ###############################
###############################################################################


###############################################################################
###############################################################################
###### Sex ratio is fixed at MLE (lambda_m = 0.3298394, lambda_f=0.6701606) ###
###############################################################################
###############################################################################
#model identification
model.id = paste("{ p(c*t), theta(t), lambda(MLE) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM = create.DM(c(1,2,3,4))
thetaDM   = create.DM(c(1,2))
lambdaDM  = matrix(, ncol=0,nrow=1)

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = c(0,0,0,0)
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(logit(MLE.model_1$est$lambda[1]))

MLE.model_1_lambda_MLE = fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,
                                    captureOFFSET,thetaOFFSET,lambdaOFFSET)
print.output(MLE.model_1_lambda_MLE )


###############################################################################
# Comparison of of the SE of the unrestricted model and the model with the sex
# ratio is fixed at MLE for the Mille Lacs Walleye dataset
###############################################################################
# parameters are - capture probabilities (p1M,p2M)
#                 - category proportions (lambda_1,lambda_2)
#                 - sub sample proportions (theta) and
#                 - population total and category totals

# MLEs for the unrestricted model
MLE_unrestricted = round(c(MLE.model_1$est$full, MLE.model_1$est$N_lambda),3)


# SE of the  estimates for the unrestricted model
se_model_1 = round(c(as.vector(MLE.model_1$se$p),MLE.model_1$se$lambda,
                      MLE.model_1$se$theta,MLE.model_1$se$N,
                      MLE.model_1$se$N_lambda),3)

# SE of the  estimates for the model when the sex ratio is fixed at the MLE
se_model_1_lambda_MLE = round(c(as.vector(MLE.model_1_lambda_MLE$se$p),
                                MLE.model_1_lambda_MLE$se$lambda,
                                MLE.model_1_lambda_MLE$se$theta,
                                MLE.model_1_lambda_MLE$se$N,
                                MLE.model_1_lambda_MLE$se$N_lambda),3)

model.cats = MLE.model_1$rawdata$category # categories in the data set

# name the parameters
parameter.names = c(paste("p1",model.cats,sep=""),paste("p2",model.cats,sep=""),
                      paste("lambda_",model.cats,sep=""),"theta_1", "theta_2",
                      "N",paste("N_",model.cats,sep=""))

Table_1 = data.frame(parameter.names, MLE_unrestricted, se_model_1,
                       se_model_1_lambda_MLE)
colnames(Table_1)=c("Parameter", "MLE", paste("SE of",MLE.model_1$model.id,
                                              sep=" "),
                    paste("SE of",MLE.model_1_lambda_MLE$model.id, sep=" "))

# Comparison of of the SE of the unrestricted model and the model with the sex
# ratio is fixed  at MLE
Table_1


###############################################################################

###############################################################################
###### Identify the violation of assumption through the gof plots  ############
###############################################################################
# generate some data where the probabilities vary between categories and/or
# sample times and fit a model where they are forced to equal.  Show the
# residuals and gof plots and see what features of these plots show extent of
# violation of assumptions.

# To generate data : need number of categories, capture probabilities,
# category proportions, sub-sample proportions and the sample size

# give the  categories in the population
categories = c("M","F")

# give the different  capture probabilities for p1M, p1F, p2M and p2F
cap.prob = c( 0.07, 0.05, 0.08, 0.1)

# give the category proportions( total should add up to 1)
lambda= c(0.6, 0.4)

# give the subsample proportions for the time 1 and 2
theta = c(0.6, 0.5)

# total numbers individuals capture for the study
sample.size = 4000


# generate data
gen.Data = generate.data(categories,cap.prob,lambda,theta,sample.size)

######### fit the required model to the generated data ################

#model identification :All capture probabilities are equal
model.id = paste("( p(.), theta(t), lambda(c) )")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM = create.DM(c(1,1,1,1)) # with equal capture probabilities
thetaDM   = create.DM(c(1,2))
lambdaDM  = create.DM(c(1))

#give the offset vectors(vectors of zeros should be given since no restriction)
captureOFFSET = c(0,0,0,0)
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

model.gen.data_1 =  fit.model(model.id,gen.Data,captureDM,thetaDM,lambdaDM,
                                        captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results (Model information,  MLEs, SE,  residual plot, etc..)
print.output(model.gen.data_1)

###############################################################################
###############################################################################



###############################################################################
###############################################################################
# #model identification : variation in capture probability by category and by time.
# # Additive model with no interaction
# model.id = paste("{ p(c+t), theta(t), lambda(0.32) }")
#
# # give  the required design matrices
# # Design matrix for capture recapture probabilities
# captureDM = matrix(c(1,1,1,
#                      1,0,1,
#                      1,1,0,
#                      1,0,0), nrow=4,byrow = TRUE)
# thetaDM   = create.DM(c(1,2)) #Design matrix for theta(sampling(sexing) fractions)
# lambdaDM  = matrix(, ncol=0,nrow=1)    #Design matrix for lambda(Category proportion)
#
# #give the offset vectors(vectors of zero's should be given since no restriction)
# captureOFFSET = c(0,0,0,0)
# thetaOFFSET   = c(0,0)
# lambdaOFFSET  = c(logit(0.32))
#
# test_ic_model_7 = fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
#                         thetaOFFSET,lambdaOFFSET)
# # print results (Model information,  MLEs, SE, residual plot, etc..)
# print.output(test_ic_model_7)
#
#
# ###############################################################################
# #model identification : variation in capture probability by category and by time.
# # Additive model with no interaction
# model.id = paste("{ p(c*t), theta(t), lambda(0.32) }")
#
# # give  the required design matrices
# # Design matrix for capture recapture probabilities
# captureDM = create.DM(c(1,2,3,4))
# thetaDM   = create.DM(c(1,2)) #Design matrix for theta(sampling(sexing) fractions)
# lambdaDM  = matrix(, ncol=0,nrow=1)    #Design matrix for lambda(Category proportion)
#
# #give the offset vectors(vectors of zero's should be given since no restriction)
# captureOFFSET = c(0,0,0,0)
# thetaOFFSET   = c(0,0)
# lambdaOFFSET  = c(logit(0.32))
#
# test_ic_model_8 = fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
#                             thetaOFFSET,lambdaOFFSET)
# # print results (Model information,  MLEs, SE, residual plot, etc..)
# print.output(test_ic_model_8)
#

###############################################################################
#model identification : variation in capture probability by category and by time.
# Additive model with no interaction
model.id = paste("{ p(c*t), theta(t), lambda(0.32) }")

# give  the required design matrices
# Design matrix for capture recapture probabilities
captureDM = matrix(, ncol=0,nrow=1)
thetaDM   = create.DM(c(1,2)) #Design matrix for theta(sampling(sexing) fractions)
lambdaDM  = matrix(, ncol=0,nrow=1)    #Design matrix for lambda(Category proportion)

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = c(0.07580455,0.01157313,0.007892997,0.020999333)
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(logit(0.3298384))

test_ic_model_9 = fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
                            thetaOFFSET,lambdaOFFSET)
# print results (Model information,  MLEs, SE, residual plot, etc..)
print.output(test_ic_model_9)























