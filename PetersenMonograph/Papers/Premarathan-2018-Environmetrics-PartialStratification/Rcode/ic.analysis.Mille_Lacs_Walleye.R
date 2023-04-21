###############################################################################
###############################################################################
### Analysis of the Mille Lacs Walleye dataset with individual covariate ######
###############################################################################

#setwd("D:\\SFU/Research/Rcode/Rcode")
#setwd("U:\\Lasantha/Research/Rcode")
#setwd("c:\\Lasantha/SFU/Research/Rcode")

# We can create different models and  the results are stored in
# IC_MLE_model_1, IC_MLE_model_2,...etc.
# Run Each Model, one at a time

source("load.R")  # load required functions and packages
set.seed(123)

options(width=350)

#Data = ic.get.data("ic.test.data.csv") # get  Mille Lacs Walleye dataset 
Data = ic.get.data("ic.millelacs-2013.csv")
str(Data)
n.rows <- length(Data$History) # number of rows of data

###############################################################################
######## take a sample of 1000 rows from the millelacs-2013 data set ##########
######## this takes less time for testing the code                   ##########
# ## This is to run the code and check 
# df <- data.frame(Data[1:4])
# sample_df <-df[sample(nrow(df), 500),]
# 
# Data$History <- as.character(sample_df$History)
# Data$counts   <- sample_df$counts 
# Data$length   <- sample_df$length
# Data$lengthsq <- sample_df$lengthsq
# str(Data)
# n.rows <- length(Data$History) # number of rows of data
###############################################################################

###############################################################################
#################### Look at histogram of lengths #############################
hist(Data$length)

# convert the data into a dataframe except "categories"
data.df <- data.frame(Data[1:4])
head(data.df)

# histogram of the lenths considering all the capture histories
hist.length <- ggplot(data.df, aes(x=length)) +
               geom_histogram(fill = 'grey40', colour = 'white')+
               ggtitle("Histogram of Lengths")
hist.length

# histogram of the lengths for each of the capture histories
hist.length.hist <- ggplot(data.df, aes(x=length)) +
                    geom_histogram(fill = 'grey40', colour = 'white')+
                    geom_histogram(aes(y = ..density..))+
                    facet_wrap( ~ History, ncol=1)+
                    ggtitle("Histogram of Lengths of Each Capture History")
hist.length.hist

# histogram of the lengths for each of the capture histories except"UU"
# consider "density" because we cannot see the histogram when there are very
# fewer counts for some capture histories compared to others
data.df.new <- data.df[data.df$History!= "UU",]
hist.length.hist.new <- ggplot(data.df.new, aes(x=length)) +
                        geom_histogram(aes(y = ..density..),fill = 'grey40', 
                                       colour = 'white')+
                        facet_wrap( ~ History, ncol=1)+
                        ggtitle("Histogram of Lengths of selected Capture Histories")
hist.length.hist.new

################### End of Histograms of lengths ##############################


###############################################################################
######################## standardize the covariate ############################

## see the correlation between length & length**2. correlations very high
cor(Data$length, Data$lengthsq) 

true.length <- Data$length  # lengths before standardizing
# standardize the covariates. just subtract the mean
# if a different standardization is used, then change it in the 
# function "ic.plot" where converting the sequence of point to regular scale
Data$length <- Data$length - mean(Data$length)
Data$lengthsq <- Data$length**2
hist(Data$length)
cor(Data$length, Data$lengthsq)  # now the correlation is very low
str(Data)

# ###### Trim the data to keep fish only between -5 and 5 of the mean   #########
# select <- Data$length > -5  & Data$length < 5
# Data$History <- Data$History[select]
# Data$counts  <- Data$counts [select]
# Data$length  <- Data$length [select]
# Data$lengthsq<- Data$lengthsq[select]
# n.rows <- length(Data$History)
# str(Data)

####################  End of standardize the covariat #########################

###############################################################################
########################## Models Fitting ##################################### 
###############################################################################

### MODEL 1 ###

# give the capture formula
captureformula <- ~length*category*time+ lengthsq*category*time+ category*time

# design Matrices
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(Data$category) *2) )
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

#model identification :All capture probabilities are equal
model.id = paste("{ p(",captureformula,"), theta(t), lambda(c) }",sep="")[2]

#Rprof("IC_MLE_model_1.out", line.profiling = TRUE)
IC_MLE_model_1 <- ic.fit.model(model.id ,Data,true.length,
                               captureformula, thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results 
ic.print.output(IC_MLE_model_1)
#Rprof(NULL)

#summaryRprof("IC_MLE_model_1.out",lines = "show")


###############################################################################
###############################################################################

### MODEL 2 ###

# give the capture formula
captureformula <- ~length+ lengthsq+ category*time

# design Matrices
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(Data$category) *2) )
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

#model identification :All capture probabilities are equal
model.id = paste("{ p(",captureformula,"), theta(t), lambda(c) }",sep="")[2]


IC_MLE_model_2 <- ic.fit.model(model.id ,Data,true.length,
                               captureformula, thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results 
ic.print.output(IC_MLE_model_2)


###############################################################################
###############################################################################
### MODEL 3 ###

# give the capture formula
captureformula <- ~category*time

# design Matrices
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(Data$category) *2) )
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

#model identification :All capture probabilities are equal
model.id = paste("{ p(",captureformula,"), theta(t), lambda(c) }",sep="")[2]

IC_MLE_model_3 <- ic.fit.model(model.id ,Data,true.length,
                               captureformula, thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results 
ic.print.output(IC_MLE_model_3)

###############################################################################
###############################################################################
### MODEL 4 ###

# give the capture formula
captureformula <- ~ category + time

# design Matrices
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(Data$category) *2) )
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

#model identification :All capture probabilities are equal
model.id = paste("{ p(",captureformula,"), theta(t), lambda(c) }",sep="")[2]

IC_MLE_model_4 <- ic.fit.model(model.id ,Data,true.length,
                               captureformula, thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results 
ic.print.output(IC_MLE_model_4)

###############################################################################
###############################################################################
### MODEL 5 ###

# give the capture formula
captureformula <-~length+ I(length^2)+ category +time
captureformula <-~length+ lengthsq+ category +time

# design Matrices
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(Data$category) *2) )
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

#model identification :All capture probabilities are equal
model.id = paste("{ p(",captureformula,"), lambda(c) }",sep="")[2]

IC_MLE_model_5 <- ic.fit.model(model.id ,Data,true.length,
                               captureformula, thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results 
ic.print.output(IC_MLE_model_5)

###############################################################################
###############################################################################
### MODEL 6 ###

# give the capture formula
captureformula <- ~length*category*time + lengthsq

# design Matrices
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = matrix(, ncol=0,nrow=1) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(Data$category) *2) )
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(logit(0.32))

#model identification :All capture probabilities are equal
model.id = paste("{ p(",captureformula,"), theta(t), lambda(0.32) }",sep="")[2]

IC_MLE_model_6 <- ic.fit.model(model.id ,Data,true.length,
                               captureformula, thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results 
ic.print.output(IC_MLE_model_6)


###############################################################################
###############################################################################
### MODEL 7 ###

# give the capture formula
captureformula <- ~length*category*time+ lengthsq + category*time

# design Matrices
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(Data$category) *2) )
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

#model identification :All capture probabilities are equal
model.id = paste("{ p(",captureformula,"), theta(t), lambda(c) }",sep="")[2]

#Rprof("IC_MLE_model_1.out", line.profiling = TRUE)
IC_MLE_model_7 <- ic.fit.model(model.id ,Data,true.length,
                               captureformula, thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results 
ic.print.output(IC_MLE_model_7)
#Rprof(NULL)

#summaryRprof("IC_MLE_model_1.out",lines = "show")
###############################################################################
######## Model comparison table.###############################################
# "ic.model.comparison.table" will produce a table which is ordered according to 
# the AICc values for all the models 
ic.model.comparison.table(models = list(IC_MLE_model_1,IC_MLE_model_2,
                                        IC_MLE_model_3)) 

###############################################################################
###############################################################################



######### Delete later  ###############################

# compare the results in MODEL 3 with MLE method 


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

MLE.model =  fit.model(model.id,Data,captureDM,thetaDM,lambdaDM,captureOFFSET,
                         thetaOFFSET,lambdaOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
print.output(MLE.model) 

######### End Delete later #############################