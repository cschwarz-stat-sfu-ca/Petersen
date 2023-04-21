###############################################################################
###### Identify the violation of assumption through the gof plots  ############
######## Generate data with violation of assumptions ( Heterogeneity) #########
###############################################################################
# generate some data where the assumptions of the model are violated and see 
# if the residual plot and  the p-value of the gof test assessment capture the 
# problem ( violation of assumptions : failure of non-death assumption, entering 
# animals to the sampled population between the first and second period, etc..)
# 
# phi = probability that animal alive at the time of the first sampling 
#       occasion is still alive and present in the population at the time of 
#       the second sampling occasion
# lambda_B = fraction of the total net births that enter the system between 
#            time 1 and time 2
# 
# p1 - capture probability at time 1
# p2 - capture probability at time 2

########################################################
source("load.R")  # load required functions and packages
options(width=350)
set.seed(123)

########################################################
#### consider 2 male groups and two female groups ######
########################################################
       N <- 100000 # population size
   theta <- c(0.8,0.6) # sub sample proportions
     phi <- .9
lambda_B <- 1.e-10 # since this is converted to log form, have to give
                   # a very small number instead of zero

##################################
########## Male group 1 ##########
##################################
      cat1 <- "M" # category
M_lambda_1 <- 0.3 # category proportion for a particular group of males

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.2
p2 <- 0.3

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0.03 # standard deviation of p1 at time 1
      s_p2 <- 0.06 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0.7 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
           D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input1 <-list(N=N,mu=mu,VCV=VCV,cat=cat1,theta=theta,lambda=M_lambda_1)


##################################
########## Male group 2 ##########
##################################
      cat1 <- "M" # category
M_lambda_2 <- 0.2 # category proportion for a particular group of males

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.3
p2 <- 0.4

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0.05 # standard deviation of p1 at time 1
      s_p2 <- 0.03 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0.8 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
           D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input2 <- list(N=N,mu=mu,VCV=VCV,cat=cat1,theta=theta,lambda=M_lambda_2)


#########################
#### Female  group 1 ####
#########################
cat2 <- "F" # category
F_lambda_1 <- 0.4 # category proportion for a particular group of females

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.2
p2 <- 0.4

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0.01 # standard deviation of p1 at time 1
      s_p2 <- 0.05 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0.6 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
           D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input3 <- list(N=N,mu=mu,VCV=VCV,cat=cat2,theta=theta,lambda=F_lambda_1)


#########################
#### Female  group 2 ####
#########################
cat2 <- "F" # category
F_lambda_2 <- 0.1 # category proportion for a particular group of females

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.3
p2 <- 0.3

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
s_p1 <- 0.02 # standard deviation of p1 at time 1
s_p2 <- 0.06 # standard deviation of p2 at time 2
s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0.5 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
           D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input4 <- list(N=N,mu=mu,VCV=VCV,cat=cat2,theta=theta,lambda=F_lambda_2)


###################################################################
#### check whether the proportion for each category is correct ####
###################################################################
tot_prop <- M_lambda_1 + M_lambda_2 + F_lambda_1 + F_lambda_2
if((tot_prop != 1)){stop("Error: proportion for each category is wrong. 
                         Total does not equal to 1")}

###############################################################################
##### fit the model for the generated data(with heterogeneity) with the ####### 
##### following design matrices and offset vectors ############################

#model identification : unrestricted model
model.id = paste("{ p(c*t), theta(t), lambda(c) }")

captureDM <- create.DM(c(1,2,3,4)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)


###############################################################################
############### fitting the model to data with heterogeneity ##################
# create list of list as input data 
# there are 4 list as input1,input2,input3,input4 since 2 male groups and two female groups 
input =list(input1,input2,input3,input4)

#################################################################
## analysis and residual plot for one data set with capture heterogeneity
gen.data.set <- ldply(input, generate.data.heterogeneity)

data <- NULL
data$History <- as.vector(gen.data.set$History)
data$counts <- gen.data.set$counts
data$category <- c(cat1,cat2)

model_res <-  fit.model(model.id,data,captureDM,thetaDM,lambdaDM,
                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
print.output(model_res) 

##################################################################

# "n.sim" is the  number of simulations of heterogeneity data sets

# simulate "n.sim" number of data set with heterogeneity and fit the model
# to each data set and save the estimates in a data frame and then calculate
# the mean vales of the estimates

#histogram of estimates of population size and the mean values of the all estimates
# the function "simulation.results" is in the file "generate.data.heterogeneity.R"
heterogeneity.results <- simulation.results(n.sim=1000,input,model.id,
                                            captureDM,thetaDM,lambdaDM,
                                            captureOFFSET,thetaOFFSET,lambdaOFFSET)
heterogeneity.results

###############################################################################
############################ End scenario 1 ###################################
###############################################################################

###############################################################################
#################  Generate data with different scenarios   ###################
###### With larger standard deviation for capture probabilities and larger  ###
###### correlation and category proportions are different considerably      ###
###### values are same as above situation except for the standard deviation ###
###############################################################################

########################################################
#### consider 2 male groups and two female groups ######
########################################################
N <- 100000 # population size
theta <- c(0.8,0.6) # sub sample proportions
phi <- .9
lambda_B <- 1.e-10 # since this is converted to log form, have to give
# a very small number instead of zero

##################################
########## Male group 1 ##########
##################################
cat1 <- "M" # category
M_lambda_1 <- 0.3# category proportion for a particular group of males

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.2
p2 <- 0.3

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ### give values between 20 -30 percent
s_p1 <- 0.05 # standard deviation of p1 at time 1
s_p2 <- 0.1 # standard deviation of p2 at time 2
s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
R_p1_p2 <- 0.7 # correlation between p1 and p2
R_p1_phi <- 0 # correlation between p1 and phi
R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
R_p2_phi <- 0 # correlation between p2 and phi
R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input1 <-list(N=N,mu=mu,VCV=VCV,cat=cat1,theta=theta,lambda=M_lambda_1)


##################################
########## Male group 2 ##########
##################################
cat1 <- "M" # category
M_lambda_2 <- 0.2 # category proportion for a particular group of males

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.3
p2 <- 0.4

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
s_p1 <- 0.1 # standard deviation of p1 at time 1
s_p2 <- 0.1 # standard deviation of p2 at time 2
s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
R_p1_p2 <- 0.8 # correlation between p1 and p2
R_p1_phi <- 0 # correlation between p1 and phi
R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
R_p2_phi <- 0 # correlation between p2 and phi
R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input2 <- list(N=N,mu=mu,VCV=VCV,cat=cat1,theta=theta,lambda=M_lambda_2)


#########################
#### Female  group 1 ####
#########################
cat2 <- "F" # category
F_lambda_1 <- 0.4 # category proportion for a particular group of females

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.2
p2 <- 0.4

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
s_p1 <- 0.06 # standard deviation of p1 at time 1
s_p2 <- 0.1 # standard deviation of p2 at time 2
s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
R_p1_p2 <- 0.6 # correlation between p1 and p2
R_p1_phi <- 0 # correlation between p1 and phi
R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
R_p2_phi <- 0 # correlation between p2 and phi
R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input3 <- list(N=N,mu=mu,VCV=VCV,cat=cat2,theta=theta,lambda=F_lambda_1)


#########################
#### Female  group 2 ####
#########################
cat2 <- "F" # category
F_lambda_2 <- 0.1 # category proportion for a particular group of females

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.3
p2 <- 0.3

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
s_p1 <- 0.1 # standard deviation of p1 at time 1
s_p2 <- 0.1 # standard deviation of p2 at time 2
s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
R_p1_p2 <- 0.5 # correlation between p1 and p2
R_p1_phi <- 0 # correlation between p1 and phi
R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
R_p2_phi <- 0 # correlation between p2 and phi
R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input4 <- list(N=N,mu=mu,VCV=VCV,cat=cat2,theta=theta,lambda=F_lambda_2)


###################################################################
#### check whether the proportion for each category is correct ####
###################################################################
tot_prop <- M_lambda_1 + M_lambda_2 + F_lambda_1 + F_lambda_2
if((tot_prop != 1)){stop("Error: proportion for each category is wrong. 
                         Total does not equal to 1")}

###############################################################################
##### fit the model for the generated data(with heterogeneity) with the ####### 
##### following design matrices and offset vectors ############################

#model identification : unrestricted model
model.id = paste("{ p(c*t), theta(t), lambda(c) }")

captureDM <- create.DM(c(1,2,3,4)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)


###############################################################################
############### fitting the model to data with heterogeneity ##################
# create list of list as input data 
# there are 4 list as input1,input2,input3,input4 since 2 male groups and two female groups 
input =list(input1,input2,input3,input4)

#################################################################
## analysis and residual plot for one data set with capture heterogeneity
gen.data.set <- ldply(input, generate.data.heterogeneity)

data <- NULL
data$History <- as.vector(gen.data.set$History)
data$counts <- gen.data.set$counts
data$category <- c(cat1,cat2)

model_res <-  fit.model(model.id,data,captureDM,thetaDM,lambdaDM,
                        captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
print.output(model_res) 

##################################################################

# "n.sim" is the  number of simulations of heterogeneity data sets

# simulate "n.sim" number of data set with heterogeneity and fit the model
# to each data set and save the estimates in a data frame and then calculate
# the mean vales of the estimates

#histogram of estimates of population size and the mean values of the all estimates
# the function "simulation.results" is in the file "generate.data.heterogeneity.R"
heterogeneity.results <- simulation.results(n.sim=1000,input,model.id,
                                            captureDM,thetaDM,lambdaDM,
                                            captureOFFSET,thetaOFFSET,lambdaOFFSET)
heterogeneity.results

#####Fitting a different model ################################################
##### fit the model for the generated data(with heterogeneity) with the ####### 
##### following design matrices and offset vectors ############################

#model identification : unrestricted model
model.id = paste("{ p(t), theta(t), lambda(c) }")

captureDM <- create.DM(c(1,1,2,2)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)

heterogeneity.results <- simulation.results(n.sim=1000,input,model.id,
                                            captureDM,thetaDM,lambdaDM,
                                            captureOFFSET,thetaOFFSET,lambdaOFFSET)
heterogeneity.results

###############################################################################
########### END of Generate data with different scenarios   ###################
###### With larger standard deviation for capture probabilities and larger  ###
###### correlation and category proportions are different considerably      ###
###############################################################################


###############################################################################
###############################################################################
########## check the code using simple Lincoln Petersen estimates #############
###############################################################################
       N <- 100000 # population size
   theta <- c(1,1) # sub sample proportions
     phi <- .99999999 # Since this is converted to logit, give a value
                      # closer to 1 instead of value 1
lambda_B <- 1.e-10 # since this is converted to log form, have to give
                   # a very small number instead of zero
    cat1 <- "M" # category
lambda_M <- 0.4 # category proportion for a particular group of males

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.2
p2 <- 0.3

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0 # standard deviation of p1 at time 1
      s_p2 <- 0 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input_M <-list(N=N,mu=mu,VCV=VCV,cat=cat1,theta=theta,lambda=lambda_M)

###########################################
    cat2 <- "F" # category
lambda_F <- 0.6 # category proportion for a particular group of males

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.2
p2 <- 0.3

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0 # standard deviation of p1 at time 1
      s_p2 <- 0 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input_F <-list(N=N,mu=mu,VCV=VCV,cat=cat2,theta=theta,lambda=lambda_F)


#### fit the model for the generated data(no heterogeneity  in the data) #################
#model identification : unrestricted model
model.id = paste("{ p(t), theta(t), lambda(c) }")

captureDM <- create.DM(c(1,1,2,2)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)

#########################################
# create list of list as input data
input =list(input_M,input_F)

n.sim =500 #the  number of simulations of heterogeneity data sets
LP.sim.data <- rdply(n.sim,fit.model.to.heterogeneity.data(input,model.id,
                                                           captureDM,thetaDM,lambdaDM,
                                                           captureOFFSET,thetaOFFSET,lambdaOFFSET))
names(LP.sim.data) <- c("sample","N","N_M","N_F", "p_1M","p_1F", "p_2M", 
                               "p_2F", "lambda_M", "lambda_F", "theta_1", "theta_2")


# mean of the estimates of the population size is displyed inside the graph
x.pos = min(LP.sim.data$N) + 
                       0.75*(max(LP.sim.data$N) - min(LP.sim.data$N))

#histogram for the estimates of the  population size
histogram_N_hat_LP <- ggplot(LP.sim.data, aes(x=N,..density.. ))+   
                        geom_histogram(aes(y=..density.. ))+
                        xlab("Estimates of N") + ylab("Density")+
                        ggtitle(paste("Histogram of population estimates for non-heterogeneity Data", 
                          "\n model: ", model.id,
                          "\n number of simulation = ", n.sim, ", initial pop size = N =",N, sep=" "))+
                        theme(axis.text=element_text(size=10,face="bold"),
                            axis.title=element_text(size=14,face="bold"), 
                            plot.title = element_text(size = 14,face="bold"))+
                        annotate("text", label = paste("mean = ",
                                 round(mean(LP.sim.data$N)),sep=""), 
                        x = x.pos, hjust = 0, y = Inf, vjust = 2, color = "darkred")

N_M <- N * lambda_M
N_F <- N * lambda_F
p_1M <- p1
p_1F <- p1
p_2M <- p2
p_2F <- p2

initial.value <- prettyNum(round(c(N,N_M,N_F,p_1M,p_1F,p_2M,p_2F, lambda_M,
                                             lambda_F, theta),3),, big.mark = ",")
names(initial.value) <- c("N","N_M","N_F", "p_1M","p_1F", "p_2M", 
                        "p_2F", "lambda_M", "lambda_F", "theta_1", "theta_2")
sim.means.LP <- prettyNum(round(colMeans(LP.sim.data[2:ncol(LP.sim.data)]),3), 
                                                                   big.mark = ",")
df.initial_sim.means <- as.data.frame(cbind(initial.value,sim.means.LP))
names(df.initial_sim.means) <- c( "initial.param.values" ,"mean_estimate")

#histogram for the estimates of the  population size
histogram_N_hat_LP
# initial parameter values and mean of the estimates using simulations
df.initial_sim.means 


###############################################################################
######## end of check the code using simple Lincoln Petersen estimates ########
###############################################################################

###############################################################################
###############################################################################
###### generate data similar to Mille Lacs Walleye data set and analysis ######
###############################################################################
###############################################################################
       N <- 205500 # population size
   theta <- c(0.9939,0.08291) # sub sample proportions
     phi <- 0.999999
lambda_B <- 1.e-10 # since this is converted to log form, have to give
                   # a very small number instead of zero

##################################
########## Male group   ##########
##################################
      cat1 <- "M" # category
M_lambda_1 <- 0.33 # category proportion for a particular group of males

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.0758
p2 <- 0.00789

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0 # standard deviation of p1 at time 1
      s_p2 <- 0 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
            D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input_MLD_M <-list(N=N,mu=mu,VCV=VCV,cat=cat1,theta=theta,lambda=M_lambda_1)

#########################
#### Female  group   ####
#########################
cat2 <- "F" # category
F_lambda_1 <- 0.67 # category proportion for a particular group of females

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.01157
p2 <- 0.02099

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0 # standard deviation of p1 at time 1
      s_p2 <- 0 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input_MLD_F <- list(N=N,mu=mu,VCV=VCV,cat=cat2,theta=theta,lambda=F_lambda_1)

###############################################################################
##### fit the model for the generated data(with heterogeneity) with the ####### 
##### following design matrices and offset vectors ############################

#model identification : unrestricted model
model.id = paste("{ p(c*t), theta(t), lambda(c) }")

captureDM <- create.DM(c(1,2,3,4)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)


###############################################################################
############### fitting the model to data with heterogeneity ##################
# create list of list as input data 
# there are 4 list as input1,input2,input3,input4 since 2 male groups and two female groups 
input =list(input_MLD_M,input_MLD_F)


###################################################################
#### check whether the proportion for each category is correct ####
###################################################################
tot_prop <- M_lambda_1 + F_lambda_1 
if((tot_prop != 1)){stop("Error: proportion for each category is wrong. 
                         Total does not equal to 1")}

#################################################################
## analysis and residual plot for one data set with capture heterogeneity
gen.data.set <- ldply(input, generate.data.heterogeneity)

data <- NULL
data$History <- as.vector(gen.data.set$History)
data$counts <- gen.data.set$counts
data$category <- c(cat1,cat2)

model_res <-  fit.model(model.id,data,captureDM,thetaDM,lambdaDM,
                        captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
print.output(model_res) 

##################################################################

# "n.sim" is the  number of simulations of heterogeneity data sets

# simulate "n.sim" number of data set with heterogeneity and fit the model
# to each data set and save the estimates in a data frame and then calculate
# the mean vales of the estimates

#histogram of estimates of population size and the mean values of the all estimates
# the function "simulation.results" is in the file "generate.data.heterogeneity.R"
heterogeneity.results <- simulation.results(n.sim=500,input,model.id,
                                            captureDM,thetaDM,lambdaDM,
                                            captureOFFSET,thetaOFFSET,lambdaOFFSET)
heterogeneity.results


#####Fitting a different model (simple petersen) ##############################
##### fit the model for the generated data(with heterogeneity) with the ####### 
##### following design matrices and offset vectors ############################

#model identification : unrestricted model
model.id = paste("{ p(t), theta(t), lambda(c) }")

captureDM <- create.DM(c(1,1,2,2)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)

heterogeneity.results <- simulation.results(n.sim=500,input,model.id,
                                            captureDM,thetaDM,lambdaDM,
                                            captureOFFSET,thetaOFFSET,lambdaOFFSET)
heterogeneity.results



###############################################################################
##End of generate data similar to Mille Lacs Walleye data set and analysis ####
###############################################################################
###############################################################################

###############################################################################
################# Generate data with different scenarios  #####################
###### With smaller standard deviation for capture probabilities and larger ###
###### correlation and category proportions are different considerably      ###
###############################################################################

########################################################
#### consider 2 male groups and two female groups ######
########################################################
       N <- 100000 # population size
   theta <- c(0.8,0.6) # sub sample proportions
     phi <- .9
lambda_B <- 1.e-10 # since this is converted to log form, have to give
# a very small number instead of zero

##################################
########## Male group 1 ##########
##################################
      cat1 <- "M" # category
M_lambda_1 <- 0.15 # category proportion for a particular group of males

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.1
p2 <- 0.05

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0.001 # standard deviation of p1 at time 1
      s_p2 <- 0.001 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0.8 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input1 <-list(N=N,mu=mu,VCV=VCV,cat=cat1,theta=theta,lambda=M_lambda_1)


##################################
########## Male group 2 ##########
##################################
      cat1 <- "M" # category
M_lambda_2 <- 0.15 # category proportion for a particular group of males

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.3
p2 <- 0.15

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0.001 # standard deviation of p1 at time 1
      s_p2 <- 0.001 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0.8 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input2 <- list(N=N,mu=mu,VCV=VCV,cat=cat1,theta=theta,lambda=M_lambda_2)


#########################
#### Female  group 1 ####
#########################
      cat2 <- "F" # category
F_lambda_1 <- 0.35 # category proportion for a particular group of females

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.05
p2 <- 0.15

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0.001 # standard deviation of p1 at time 1
      s_p2 <- 0.001 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0.8 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input3 <- list(N=N,mu=mu,VCV=VCV,cat=cat2,theta=theta,lambda=F_lambda_1)


#########################
#### Female  group 2 ####
#########################
      cat2 <- "F" # category
F_lambda_2 <- 0.35 # category proportion for a particular group of females

## give the mean, standard deviation and correlation between each pair of   ###
## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###

### mean for p1, p2, phi, lambda_B ###
p1 <- 0.15
p2 <- 0.45

mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 

### standard deviation of each variable ###
      s_p1 <- 0.001 # standard deviation of p1 at time 1
      s_p2 <- 0.001 # standard deviation of p2 at time 2
     s_phi <- 0 # standard deviation of phi
s_lambda_B <- 0 # standard deviation of lambda_B 

### correlation between each pair of variables ####
       R_p1_p2 <- 0.8 # correlation between p1 and p2
      R_p1_phi <- 0 # correlation between p1 and phi
 R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
      R_p2_phi <- 0 # correlation between p2 and phi
 R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
R_phi_lambda_B <- 0 # correlation between phi and lambda_B

## correlation matrix ##
R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
               R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
               R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
               R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)

## vector with standard deviation of each variable ##
standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
D <- diag(standard_dev) # create a diagonal matrix 

## convert correlation matrix to variance-covariance matrix ##
VCV <- D%*%R%*%D 

input4 <- list(N=N,mu=mu,VCV=VCV,cat=cat2,theta=theta,lambda=F_lambda_2)


###################################################################
#### check whether the proportion for each category is correct ####
###################################################################
tot_prop <- M_lambda_1 + M_lambda_2 + F_lambda_1 + F_lambda_2
tot_prop
#if((tot_prop != 1)) {stop("Error: proportion for each category is wrong. 
#                         Total does not equal to 1")}

###############################################################################
##### fit the model for the generated data(with heterogeneity) with the ####### 
##### following design matrices and offset vectors ############################

#model identification : unrestricted model
model.id = paste("{ p(c*t), theta(t), lambda(c) }")

captureDM <- create.DM(c(1,2,3,4)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)


###############################################################################
############### fitting the model to data with heterogeneity ##################
# create list of list as input data 
# there are 4 list as input1,input2,input3,input4 since 2 male groups and two female groups 
input =list(input1,input2,input3,input4)

#################################################################
## analysis and residual plot for one data set with capture heterogeneity
gen.data.set <- ldply(input, generate.data.heterogeneity)

data <- NULL
data$History <- as.vector(gen.data.set$History)
data$counts <- gen.data.set$counts
data$category <- c(cat1,cat2)

model_res <-  fit.model(model.id,data,captureDM,thetaDM,lambdaDM,
                        captureOFFSET,thetaOFFSET,lambdaOFFSET)
# print results (Model information, MLEs, SE,residual plot, etc..)
print.output(model_res) 

##################################################################

# "n.sim" is the  number of simulations of heterogeneity data sets

# simulate "n.sim" number of data set with heterogeneity and fit the model
# to each data set and save the estimates in a data frame and then calculate
# the mean vales of the estimates

#histogram of estimates of population size and the mean values of the all estimates
# the function "simulation.results" is in the file "generate.data.heterogeneity.R"
heterogeneity.results <- simulation.results(n.sim=1000,input,model.id,
                                            captureDM,thetaDM,lambdaDM,
                                            captureOFFSET,thetaOFFSET,lambdaOFFSET)
heterogeneity.results

###################### Fitting a different model ##############################
##### fit the model for the generated data(with heterogeneity) with the ####### 
##### following design matrices and offset vectors ############################

#model identification : unrestricted model
model.id = paste("{ p(t), theta(t), lambda(c) }")

captureDM <- create.DM(c(1,1,2,2)) 
thetaDM   <- create.DM(c(1,2)) 
lambdaDM  <- create.DM(c(1)) 

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0) 
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- c(0)

heterogeneity.results <- simulation.results(n.sim=1000,input,model.id,
                                            captureDM,thetaDM,lambdaDM,
                                            captureOFFSET,thetaOFFSET,lambdaOFFSET)
heterogeneity.results


###############################################################################
############### END of Generate data with different scenarios   ###############
###### With smaller standard deviation for capture probabilities and larger ###
###### correlation and category proportions are different considerably      ###
###############################################################################


