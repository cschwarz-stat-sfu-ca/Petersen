###############################################################################
###############################################################################
### Generate data with individual covariate and given formula #################
##############  ( with two categories)   #####################################
############################################################################### 

########################################################
########################################################
source("load.R")  # load required functions and packages
set.seed(123)

N <- 10000 # population size 
theta <- c(0.3,0.1) # sub sample proportions
category <- c("M","F") # category

# category proportions. sum should be equal to 1 and the and the values 
# corresponds to the categories defined above
lambda <- c(0.3,0.7) 

# give the capture formula
#captureformula <- ~ category+ time  ## simple model:  model 1
#captureformula <- ~length + category+ time   ## model 2
captureformula <- ~length+ lengthsq + category+ time  # model 3
#captureformula <- ~length*category*time+ lengthsq*category*time+ category*time # model 4


# number of beta parameters related to capture formula
#  need add 1, because of the intercept
n.beta.param.cap <- length(colnames(attr(terms(captureformula), "factors"))) +1


### convert to logit scale  #####

beta.logit.theta <- logit(theta)
beta.logit.lambda <- lambda.to.cumulative.logit.lambda(lambda) 
# give the beta parameters according to the capture formula
#beta.logit.p <- c(-2, -0.1,-2 )  ## this is for the simple model:  model 1
#beta.logit.p <- c(-2, -0.01,-0.1,-0.5 )  ## this is for the simple model:  model 2
beta.logit.p <- c(-3,0.1,-0.004,-0.5,-0.8)  # this is for nice curve with length square : for model 3
#beta.logit.p <- c(-3.4 , -0.4 ,-0.6,  -1.9, -0.07, 0.5, 0.1, 2.4, 0.05, 0.01, -0.2, 0.01) #  for model 4




# logit estimates are put in the same order as used in the other functions
logit.est <- c( beta.logit.p,beta.logit.lambda,beta.logit.theta)


###################################

simulate_data<- NULL
cat_size <- round(N*lambda) # number of histories to be generated from each category

simulate_data$length <- c(rnorm(N, mean=20,sd=2.8)) # this is closer to Mille Lacs Data set

simulate_data$counts    <- rep(1,N) # counts of the corresponding capture histories

history <- rep("00", N) # initialize histories as "00"
simulate_data$History  <- history   # capture histories of the raw data
simulate_data$category <- category 
simulate_data$lengthsq <- simulate_data$length^2 # length square
str(simulate_data)

# design MatriX FOR CAPTURE PROBABILITIES
captureDM <- ic.create.DM(simulate_data,captureformula)

# design Matrices FOR THETA AND LAMBDA
thetaDM   = create.DM(c(1,2)) 
lambdaDM  = create.DM(c(1)) 


#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET = matrix(0, nrow=N, ncol= (length(simulate_data$category) *2) )
thetaOFFSET   = c(0,0)
lambdaOFFSET  = c(0)

#model identification :All capture probabilities are equal
model.id = paste("{ p(",captureformula,"), theta(t), lambda(c) }",sep="")[2]


########################################## 
#estimates in regular form
 regular.est <- ic.unpack.parm(logit.est,simulate_data,
                            captureDM,thetaDM,lambdaDM,
                            captureOFFSET,thetaOFFSET,lambdaOFFSET)

# capture_Prob has four columns as M_time1, F_time_1, M_time2 and F_time2
capture_Prob <- regular.est$p
head(regular.est$p)
capture_Prob[99990:100000,]

# this is just to look at the plots( capture probabilities) with generated
#lengths without standardizing.
ic.generated.data.plots (logit.est,simulate_data,
                          captureformula,thetaDM,lambdaDM,
                          thetaOFFSET,lambdaOFFSET)
############################################


###############################################################################
######## Generate data at time 1 ##############################################

########## generate data for category "M" ##########

history_M <- history[1:cat_size[1]]

# generate a vector indicating the capture status(0/1 from Bernoulli process)
cap_M_1 <- rbinom(cat_size[1],1,capture_Prob[1:cat_size[1] , 1]) # 1-capture  , 0-not capture

# generate a vector indicating the category status(0/1 from Bernoulli process)
sex_M_1 <- rbinom(cat_size[1],1,theta[1]) # 1-stratify  , 0-not stratify

# 0-not captured, "U"- captured but not stratified, "category[1]"- stratified (eg:"M")
substr(history_M[cap_M_1 & cap_M_1 ],1,1) <- "U" 
substr(history_M[cap_M_1 & sex_M_1 ],1,1) <- category[1] 


########## generate data for category "F" ##########

history_F <- history[(cat_size[1]+1):N]

# generate a vector indicating the capture status(0/1 from Bernoulli process)
cap_F_1 <- rbinom(cat_size[2],1,capture_Prob[(cat_size[1]+1):N , 2]) # 1-capture  , 0-not capture

# generate a vector indicating the category status(0/1 from Bernoulli process)
sex_F_1 <- rbinom(cat_size[2],1,theta[1]) # 1-stratify  , 0-not stratify

# 0-not captured, "U"- captured but not stratified, "cat"- stratified (eg:"M")
substr(history_F[cap_F_1 & cap_F_1 ],1,1) <- "U" 
substr(history_F[cap_F_1 & sex_F_1 ],1,1) <- category[2] 



#####################################################################
########## Generated data at time 2 #################################

########## generate data for category "M" ##########

# generate a vector indicating the capture status(0/1 from Bernoulli process)
cap_M_2 <- rbinom(cat_size[1],1,capture_Prob[1:cat_size[1] , 3]) # 1-capture  , 0-not capture

# generate a vector indicating the category status(0/1 from Bernoulli process)
sex_M_2 <- rbinom(cat_size[1],1,theta[2]) # 1-stratify  , 0-not stratify

# 0-not captured, "U"- captured but not stratified, "category[1]"- stratified (eg:"M")

substr(history_M[cap_M_2 & cap_M_2 ],2,2) <- "U"  
substr(history_M[cap_M_1 & sex_M_1 & cap_M_2 ],2,2) <- category[1] 
substr(history_M[!cap_M_1 & cap_M_2 & sex_M_2 ],2,2) <- category[1] 


########## generate data for category "F" ##########

# generate a vector indicating the capture status(0/1 from Bernoulli process)
cap_F_2 <- rbinom(cat_size[2],1,capture_Prob[(cat_size[1]+1):N , 4]) # 1-capture  , 0-not capture

# generate a vector indicating the category status(0/1 from Bernoulli process)
sex_F_2 <- rbinom(cat_size[2],1,theta[2]) # 1-stratify  , 0-not stratify

# 0-not captured, "U"- captured but not stratified, "category[1]"- stratified (eg:"M")

substr(history_F[cap_F_2 & cap_F_2 ],2,2) <- "U" 
substr(history_F[cap_F_1 & sex_F_1 & cap_F_2 ],2,2) <- category[2] 
substr(history_F[!cap_F_1 & cap_F_2 & sex_F_2 ],2,2) <- category[2] 


###############################################################################
##############  generates data as a list (without the history "00) ############

#combine generates histories
generated_history <- c(history_M, history_F)

# summary of generated histories
table(history_M)
table(history_F)
table(generated_history)


sim.data <- NULL

sim.data$History  <- generated_history[(generated_history!="00")]
sim.data$length   <- simulate_data$length[(generated_history!="00")]
sim.data$lengthsq <- sim.data$length^2
sim.data$counts   <- simulate_data$counts[(generated_history!="00")]
sim.data$category <- category

str(sim.data) 
sim.table <- table(sim.data$History)
print("simulated histories ")
print( sim.table)
paste("total simulated histories = ", sum(sim.table))

################################################# 
# this is just to check the simple model. ie considering no length  
###           captureformula <- ~ category* time   #############
## Expected counts for the histories
n0F <- N*lambda[2]*(1-capture_Prob[1,2])*capture_Prob[1,4]*theta[2]
n0M <- N*lambda[1]*(1-capture_Prob[1,1])*capture_Prob[1,3]*theta[2]
n0U <- N*( lambda[1]*(1-capture_Prob[1,1])*capture_Prob[1,3]* (1-theta[2]) +  
             lambda[2]*(1-capture_Prob[1,2])*capture_Prob[1,4]* (1-theta[2])  ) 

nF0 <- N*lambda[2]*capture_Prob[1,2]*theta[1]*(1-capture_Prob[1,4])
nFF <- N*lambda[2]*capture_Prob[1,2]*theta[1]*capture_Prob[1,4]
nM0 <- N*lambda[1]*capture_Prob[1,1]*theta[1]*(1-capture_Prob[1,3])
nMM <- N*lambda[1]*capture_Prob[1,1]*theta[1]*capture_Prob[1,3]

nU0 <- N*( lambda[1]*capture_Prob[1,1]*(1-theta[1])*(1-capture_Prob[1,3]) +  
             lambda[2]*capture_Prob[1,2]*(1-theta[1])*(1-capture_Prob[1,4]) ) 
nUU <- N*( lambda[1]*capture_Prob[1,1]*(1-theta[1])*capture_Prob[1,3] +  
             lambda[2]*capture_Prob[1,2]*(1-theta[1])*capture_Prob[1,4] ) 



expectedcounts <- round(c( n0F,n0M ,n0U, nF0, nFF,nM0,nMM,nU0 , nUU), 2)
names(expectedcounts) <-c("0F","0M" ,"0U", "F0", "FF","M0","MM","U0" , "UU"  )
print("Expected histories")
print(expectedcounts)
paste(" Expected total histories = ", sum(expectedcounts ))

############
print("simulated histories ")
print( sim.table)
paste("total simulated histories = ", sum(sim.table))


###############################################################################
#################################################################################


## test with generated data ###############################

n.rows <- length(sim.data$History) # number of rows of data
captureOFFSET = matrix(0, nrow=n.rows, ncol= (length(sim.data$category) *2) )


true.length <- sim.data$length  # lengths before standardizing
# standardize the length
sim.data$length   <- sim.data$length - mean(sim.data$length)
sim.data$lengthsq <- sim.data$length**2

cor(sim.data$length, sim.data$lengthsq)


IC_MLE_model_Test <- ic.fit.model(model.id ,sim.data,true.length,
                               captureformula, thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results 
ic.print.output(IC_MLE_model_Test)

###########################################################


# Data<- NULL
# Data$History <- sim.data$History 
# Data$counts <- sim.data$counts 
# Data$length <- sim.data$length 
# Data$lengthsq <- sim.data$lengthsq 
# Data$category <- sim.data$category 
# 
# str(Data)
# 
# str(sim.data)
