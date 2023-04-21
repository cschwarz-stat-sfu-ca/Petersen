############################################################################### 
#####                           Analysis of                              ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################

# We can create different models and  the results are stored in
# MLE_PSkSCR_model_1 , MLE_PSkSCR_model_2, ,...etc.
# Run Each Model, one at a time

#############  path below should points to the Functions directory  ###########
# setwd("U:\\Lasantha/Research/Rcode")
# setwd("C:\\Lasantha/SFU/Research/Rcode")
# setwd("D:\\SFU/Research/Rcode")

source("load.R")  # load required functions and packages
set.seed(123)

###############################################################################
# get  data to analysis
data <- PSkSCR.get.data("PSkSCR_Test_data_1.csv")
str(data)

# temp = get.data("Mille_Lacs_Walleye.csv") # get  Mille Lacs Walleye dataset 
# str(temp)
# newdata <- NULL
# newdata$history <- Data$History
# newdata$counts <- Data$counts
# newdata$category <- Data$category
# data <- newdata
# str(data)
# ###############################################################################


###############################################################################
###############################################################################

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
# if capture probabilities are not fixed then the capture offset is
# c(0,0,0,0)


###############################################################################
###############################################################################

#model identification : unrestricted model
model.id <- paste("{ p(c*t), theta(t), lambda(c), nu(t) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ),  lambda(Category proportion)
# and loss on capture probabilities 
captureDM <- create.DM(c(1,2,3,4,5,6)) 
thetaDM   <- create.DM(c(1,2,3))
lambdaDM  <- create.DM(c(1)) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(0)
p_lossOFFSET  <- c(0,0,0)

MLE_PSkSCR_model_1 <- PSkSCR.fit.model(model.id,data,
                                       captureDM,thetaDM,lambdaDM,p_lossDM,
                                       captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
PSkSCR.print.output(MLE_PSkSCR_model_1)

# # ##################  this is just for checking codes
# # 
# new.data <- NULL
# length.new.data <- length(MLE_PSkSCR_model_1$obs.exp.counts$History)-1
# new.data$history <- as.character( MLE_PSkSCR_model_1$obs.exp.counts$History[1: length.new.data])
# new.data$counts <- MLE_PSkSCR_model_1$obs.exp.counts$Expected.counts[1:length.new.data]
# new.data$category<- data$category
# 
# expected.data.MLE_PSkSCR_model_1 <- PSkSCR.fit.model(model.id,data=new.data,
#                                        captureDM,thetaDM,lambdaDM,
#                                        captureOFFSET,thetaOFFSET,lambdaOFFSET,
#                                        dataprepare="YES")
# ## need to use dataprepare="YES" because in the data preparation step
# ## it round off the  expected counts to integers.
# ## this new.data is already prepared since  each history is
# ## duplicated maximum two times(if there is negative value for the counts)
# 
# ### data preparation is only needed for raw data ###
# 
# PSkSCR.print.output(expected.data.MLE_PSkSCR_model_1)
# # 
# # ################## end of  checking codes

#parametric bootstrap goodness of fit
PSkSCR.bootstrap.gof.model_1 <- PSkSCR.parametric.bootstrap(MLE_PSkSCR_model_1,n.boots=20)
PSkSCR.bootstrap.gof.model_1



###############################################################################
###############################################################################


#model identification : Capture probability vary by category, but not by time
model.id <- paste("{ p(c), theta(t), lambda(c), nu(t) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) , lambda(Category proportion)
# and loss on capture probabilities 
captureDM <- create.DM(c(1,2,1,2,1,2)) 
thetaDM   <- create.DM(c(1,2,3))
lambdaDM  <- create.DM(c(1)) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(0)
p_lossOFFSET  <- c(0,0,0)

MLE_PSkSCR_model_2 <- PSkSCR.fit.model(model.id,data,
                                       captureDM,thetaDM,lambdaDM,p_lossDM,
                                       captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
PSkSCR.print.output(MLE_PSkSCR_model_2) 

# #parametric bootstrap goodness of fit
PSkSCR.bootstrap.gof.model_2 <- PSkSCR.parametric.bootstrap(MLE_PSkSCR_model_2,n.boots=20)
PSkSCR.bootstrap.gof.model_2

###############################################################################
###############################################################################


#model identification : Category proportions are equal (0.5)
model.id <- paste("{ p(c*t), theta(t), lambda(0.5), nu(t) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) , lambda(Category proportion)
# and loss on capture probabilities 
captureDM <- create.DM(c(1,2,3,4,5,6)) 
thetaDM   <- create.DM(c(1,2,3))
lambdaDM  = matrix(, ncol=0,nrow=1) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(logit(0.5))
p_lossOFFSET  <- c(0,0,0)

MLE_PSkSCR_model_3 <- PSkSCR.fit.model(model.id,data,
                                       captureDM,thetaDM,lambdaDM,p_lossDM,
                                       captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
PSkSCR.print.output(MLE_PSkSCR_model_3) 

# #parametric bootstrap goodness of fit
PSkSCR.bootstrap.gof.model_3 <- PSkSCR.parametric.bootstrap(MLE_PSkSCR_model_3,n.boots=30)
PSkSCR.bootstrap.gof.model_3

###############################################################################
###############################################################################

#model identification : Capture probability vary by time, but not by category
model.id <- paste("{ p(t), theta(t), lambda(c), nu(t) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) , lambda(Category proportion)
# and loss on capture probabilities 
captureDM <- create.DM(c(1,1,2,2,3,3)) 
thetaDM   <- create.DM(c(1,2,3))
lambdaDM  <- create.DM(c(1)) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(0)
p_lossOFFSET  <- c(0,0,0)

MLE_PSkSCR_model_4 <- PSkSCR.fit.model(model.id,data,
                                       captureDM,thetaDM,lambdaDM,p_lossDM,
                                       captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
PSkSCR.print.output(MLE_PSkSCR_model_4) 

# #parametric bootstrap goodness of fit
PSkSCR.bootstrap.gof.model_4 <- PSkSCR.parametric.bootstrap(MLE_PSkSCR_model_4,n.boots=100)
PSkSCR.bootstrap.gof.model_4

###############################################################################
###############################################################################

#model identification : Capture probability not vary by time and  by category
model.id <- paste("{ p(.), theta(t), lambda(c), nu(t) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ), lambda(Category proportion)
# and loss on capture probabilities 
captureDM <- create.DM(c(1,1,1,1,1,1)) 
thetaDM   <- create.DM(c(1,2,3))
lambdaDM  <- create.DM(c(1))
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(0)
p_lossOFFSET  <- c(0,0,0)

MLE_PSkSCR_model_5 <- PSkSCR.fit.model(model.id,data,
                                       captureDM,thetaDM,lambdaDM,p_lossDM,
                                       captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
PSkSCR.print.output(MLE_PSkSCR_model_5) 

# #parametric bootstrap goodness of fit
PSkSCR.bootstrap.gof.model_5 <- PSkSCR.parametric.bootstrap(MLE_PSkSCR_model_5,n.boots=500)
PSkSCR.bootstrap.gof.model_5

###############################################################################
###############################################################################


#model identification : Category proportions are fixed at MLEs
model.id <- paste("{ p(t), theta(t), lambda(0.7), nu(t) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) , lambda(Category proportion)
# and loss on capture probabilities 
captureDM <- create.DM(c(1,1,2,2,3,3)) 
thetaDM   <- create.DM(c(1,2,3))
lambdaDM  = matrix(, ncol=0,nrow=1) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(logit(0.7))
p_lossOFFSET  <- c(0,0,0)

MLE_PSkSCR_model_6 <- PSkSCR.fit.model(model.id,data,
                                       captureDM,thetaDM,lambdaDM,p_lossDM,
                                       captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
PSkSCR.print.output(MLE_PSkSCR_model_6) 

# #parametric bootstrap goodness of fit
PSkSCR.bootstrap.gof.model_6 <- PSkSCR.parametric.bootstrap(MLE_PSkSCR_model_6,n.boots=10)
PSkSCR.bootstrap.gof.model_6
###############################################################################
###############################################################################

## Compare all the given models
PSkSCR.model.comparison.table(list(MLE_PSkSCR_model_1, MLE_PSkSCR_model_2,
                                   MLE_PSkSCR_model_3, MLE_PSkSCR_model_4,
                                   MLE_PSkSCR_model_5,MLE_PSkSCR_model_6))



###############################################################################
###############################################################################


###############################################################################
###############################################################################
###### Sex ratio is fixed at MLE (lambda_m = 0.7315861, lambda_f=0.2684139) ###
###############################################################################
###############################################################################


#model identification : Category proportions are fixed at MLEs
model.id <- paste("{ p(c*t), theta(t), lambda(MLE), nu(t) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) , lambda(Category proportion)
# and loss on capture probabilities 
captureDM <- create.DM(c(1,2,3,4,5,6)) 
thetaDM   <- create.DM(c(1,2,3))
lambdaDM  = matrix(, ncol=0,nrow=1) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(logit(MLE_PSkSCR_model_1$est$lambda[1]))
p_lossOFFSET  <- c(0,0,0)

MLE_PSkSCR_model_7 <- PSkSCR.fit.model(model.id,data,
                                       captureDM,thetaDM,lambdaDM,p_lossDM,
                                       captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
PSkSCR.print.output(MLE_PSkSCR_model_7) 

###############################################################################
###############################################################################
