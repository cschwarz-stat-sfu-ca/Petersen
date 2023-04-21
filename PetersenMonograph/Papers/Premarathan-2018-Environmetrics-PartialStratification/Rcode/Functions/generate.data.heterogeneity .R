###############################################################################
###############################################################################
######## Generate data with violation of assumptions ( Heterogeneity) #########
###############################################################################
# generate some data where the assumptions of the model are violated and see 
# if the residual plot and  the p-value of the gof test assessment capture the 
# problem ( violation of assumptions : failuer of non-death assumption, entering 
# animals to the sampled population between the first and second period, etc.. )
# 
# phi = probability that animal alive at the time of the first sampling occasion 
#       is still alive and present in the population at the time of the second 
#       sampling occasion
# lambda_B = fraction of the total net births that enter the system between 
#            time 1 and time 2
# 
# p1 - capture probability at time 1
# p2 - capture probability at time 2

# Function:generate.data.heterogeneity
# Input : a list with N-population size, vector of means for p1,p2, phi and lambda_b,
#         variance  covariance matrix, category, sub sample proportions(theta), 
#         category(group) proportion (lambda) 
# output : generated data for observable capture histories with heterogeneity


generate.data.heterogeneity = function(input){
           N <- input$N 
          mu <- input$mu
         VCV <- input$VCV
         cat <- input$cat
       theta <- input$theta
      lambda <- input$lambda
      
  ###### generate multivariate normal r.v. for p1,p2,phi and lambda_B #########
  
  # formula representing the transformation.
  logitprob <-  llply(c(paste("~(log(x",1:(nrow(VCV)-1),")","-log(1-x",1:(nrow(VCV)-1),"))",sep=""),
                        paste("~log(x",nrow(VCV),")",sep="")),as.formula)
  
  #  covariance matrix in logit form
  logit_VCV = deltamethod(logitprob,mu, VCV, ses=FALSE) 
  
  # means of the multivariate normal distribution in logit/log form
  logit_mu <- c(logit(mu[1:(nrow(VCV)-1)]), log(mu[nrow(VCV)])) 
  
  # generate multivariate normal r.v. in logit form
  simdata_logit <- mvrnorm(N*lambda, mu=logit_mu, Sigma=logit_VCV)
  
  # convert to the regular form (p1,p2,phi and lambda_B)
  simData <- cbind(expit(simdata_logit[ ,1:(nrow(VCV)-1)]), 
                                  exp(simdata_logit[ ,nrow(VCV)]))

  
  ######## Generate data at time 1 ################################
  
  size <- round(N*lambda) # number of histories to be generated
  
  history <- rep("00", size) # initialize histories as "00"
  
  # generate a vector indicating the capture status(0/1 from bernoulli process)
  cap_1 <- rbinom(size,1,simData[ , 1]) # 1-capture  , 0-not capture
  
  # generate a vector indicating the category status(0/1 from bernoulli process)
  sex_1 <- rbinom(size,1,theta[1]) # 1-stratify  , 0-not stratify
  
  # 0-not captured, "U"- captured but not stratified, "cat"- stratified (eg:"M")
  substr(history[cap_1 & cap_1 ],1,1) <- "U" 
  substr(history[cap_1 & sex_1 ],1,1) <- cat 
  
  
  ########### Generated data at time 2 ############################
  
  # generate a vector indicating the capture status(0/1 from bernoulli process)
  cap_2 <- rbinom(size,1,simData[ , 2]) # 1-capture  , 0-not capture
  
  # generate a vector indicating the category status(0/1 from bernoulli process)
  sex_2 <- rbinom(size,1,theta[2]) # 1-stratify  , 0-not stratify
  
  # generate a vector indicating whether the individual is alive or 
  # dead (0/1 from bernoulli process)
  alive_2 <- rbinom(size,1,simData[ , 3]) # 1-alive  , 0-dead
  
  substr(history[cap_2 & alive_2 ],2,2) <- "U"
  substr(history[cap_1 & sex_1 & cap_2 & alive_2],2,2) <- cat
  substr(history[!cap_1 & cap_2 & alive_2 & sex_2 ],2,2) <- cat
  
  
  ####### generate data for new births at time 2 #######
  
  # generate lambda_B for every individual. Then sum these and generate poisson 
  # random variables for the new recruitments. here the possible histories
  # are "00", "0U", and "0'cat'" (eg: "00", "0U" and "0M")
  
  sum_poisson_lambda <- sum(simData[ , 4]) # sum of means
  
  n_new_births <- rpois(1,sum_poisson_lambda)
  
  if(n_new_births > 0){ # generate possible histories for new individuals
    new_history <- rep("00", n_new_births) # initialize histories as "00" for new births
  
    # generate a vector indicating the capture status for 
    # new births ( 0/1 from bernoulli process)
    cap_2_new <- rbinom(new_history,1,mu[2]) # 1-capture  , 0-not capture
  
    # generate a vector indicating the category status for 
    # new births(0/1 from bernoulli process)
    sex_2_new <- rbinom(new_history,1,theta[2]) # 1-stratify  , 0-not stratify
  
    # 0-not captured, "U"- captured but not stratified, "cat"- stratified (eg : "M")
    substr(new_history[cap_2_new & cap_2_new ],2,2) <- "U" 
    substr(new_history[cap_2_new & sex_2_new ],2,2) <- cat 
    
    gen.history <- c( history, new_history) # generated histories
    
  } else {gen.history <- history }
  
  #####################################################
  
  # generated capture histories
  summary <- table(gen.history)
  data.summary <- as.data.frame(summary)
  names(data.summary) <- c("History", "counts")
  
  #remove "00" history. This gives counts for observable capture histories
  if(data.summary[1,1]=="00"){
    obs.summary <- data.summary[2:nrow(data.summary),]
  } else {obs.summary <- data.summary}
  

  return(obs.summary)
  
} # end of "generate.data.heterogeneity"


###############################################################################
###############################################################################
## fit model to generated data with heterogeneity

# Function : fit.model.to.heterogeneity.data
# Inputs : input as a list of list 
#           (inner list contains N, mu, VCV, cat, theta, lambda)
#          and modelid, design matrices and offset vectors for fitting model
# Output : estimated for parameters
#          ( population size, population size for each category,
#            capture probabilities, category proportions, sub-sample proportions)

fit.model.to.heterogeneity.data <- function(input, model.id,captureDM,thetaDM,lambdaDM,
                                            captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  gen.data <- ldply(input , generate.data.heterogeneity)
  
  data <- NULL
  data$History <- as.vector(gen.data$History)
  data$counts <- gen.data$counts
  data$category <- c(cat1,cat2)
  
  res <-  fit.model(model.id,data,captureDM,thetaDM,lambdaDM,
                                         captureOFFSET,thetaOFFSET,lambdaOFFSET)
         N <- res$est$N
       N_M <- res$est$N_lambda[1]
       N_F <- res$est$N_lambda[2]
      p_1M <- res$est$p[1,1]
      p_1F <- res$est$p[2,1]
      p_2M <- res$est$p[1,2]
      p_2F <- res$est$p[2,2]
  lambda_M <- res$est$lambda[1]
  lambda_F <- res$est$lambda[2]
   theta_1 <- res$est$theta[1]
   theta_2 <- res$est$theta[2]
  
  return(c(N, N_M, N_F, p_1M, p_1F, p_2M, p_2F, lambda_M, lambda_F, theta_1, theta_2))
} # end of  fit.model.to.heterogeneity.data 


##############################################################################
##############################################################################
# simulate "n.sim" number of data set with capture heterogeneity and fit the 
# model to each data set and save the estimates in a data frame and then 
# calculate the mean vales of the estimates

# Function : simulation.results
# Inputs : number of data sets to be simulated, input as a list of list 
#           (inner list contains N, mu, VCV, cat, theta, lambda)
#          and modelid, design matrices and offset vectors for fitting model
# OUtput : histogram of estimates of population size and the 
#          mean values of the all estimates

simulation.results <- function(n.sim=1000,input,model.id,captureDM,thetaDM,lambdaDM,
                                 captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  estimates.sim.data <- rdply(n.sim,fit.model.to.heterogeneity.data(input, 
                                             model.id,captureDM,thetaDM,lambdaDM,
                                             captureOFFSET,thetaOFFSET,lambdaOFFSET))
  names(estimates.sim.data) <- c("sample","N","N_M","N_F", "p_1M","p_1F", "p_2M", 
                                 "p_2F", "lambda_M", "lambda_F", "theta_1", "theta_2")
  
  # mean of the estimates of the population size is displyed inside the graph
  x.pos = min(estimates.sim.data$N) + 
                   0.75*(max(estimates.sim.data$N) - min(estimates.sim.data$N))
  
  #histogram for population estimates
  histogram_N_hat <- ggplot(estimates.sim.data, aes(x=N,..density.. ))+   
                      geom_histogram(aes(y=..density.. ))+
                      xlab("Estimates of N") + ylab("Density")+
                      ggtitle(paste("Histogram of population estimates for heterogeneity Data", 
                          "\n model: ", model.id,
                          "\n number of simulation = ", n.sim, ", initial pop size = N =",N, sep=" "))+
                      theme(axis.text=element_text(size=10,face="bold"),
                            axis.title=element_text(size=14,face="bold"), 
                            plot.title = element_text(size = 14,face="bold"))+
                      annotate("text", label = paste("mean = ",
                            round(mean(estimates.sim.data$N)),sep=""), 
                            x = x.pos, hjust = 0, y = Inf, vjust = 2, color = "darkred")
 
  # mean values of the all estimates
  sim.means <- as.data.frame(prettyNum(round(colMeans(
                             estimates.sim.data[2:ncol(estimates.sim.data)]),3), big.mark = ","))
  names(sim.means) <- c( "mean_estimate")
  
  # returs histogram of estimates of population size and the mean values of the all estimates
  results <- NULL
  results$histogram <- histogram_N_hat
  results$sim.means <- sim.means
  
  return(results)
  
}



###############################################################################
############ Test the function "generate.data.heterogeneity" ##################
###############################################################################
# test.generate.data.heterogeneity <-function(){
#   N <- 1000 # population size
#   theta <- c(0.8,0.6) # sub sample proportions
#   phi <- 0.9
#   lambda_B <- 1.e-10 # since this is converted to log form  have  to give 
#                        a very small number instead of zero,
#   
#   # consider male group
#   cat <- "M" # category
#   M_lambda_1 <- 0.3 # category proportion for a particular group of males
#   
#   ## give the mean, standard deviation and correlation between each pair of   ###
#   ## variables (p1, p2, phi, lambda_B) for generating data with heterogeneity ###
#   
#   ### mean for p1, p2, phi, lambda_B ###
#   p1 <- 0.2
#   p2 <- 0.3
#   
#   mu <- c(p1, p2, phi, lambda_B) #  mean values as a vector 
#   
#   ### standard deviation of each variable ###
#         s_p1 <- 0.03 # standard deviation of p1 at time 1
#         s_p2 <- 0.06 # standard deviation of p2 at time 2
#        s_phi <- 0 # standard deviation of phi
#   s_lambda_B <- 0 # standard deviation of lambda_B 
#   
#   ### correlation between each pair of variables ####
#          R_p1_p2 <- 0.7 # correlation between p1 and p2
#         R_p1_phi <- 0 # correlation between p1 and phi
#    R_p1_lambda_B <- 0 # correlation between p1 and lambda_B
#         R_p2_phi <- 0 # correlation between p2 and phi
#    R_p2_lambda_B <- 0 # correlation between p2 and lambda_B
#   R_phi_lambda_B <- 0 # correlation between phi and lambda_B
#   
#   ## correlation matrix ##
#   R <- matrix(c( 1,            R_p1_p2,      R_p1_phi,      R_p1_lambda_B,
#                  R_p1_p2,      1,            R_p2_phi,      R_p2_lambda_B,
#                  R_p1_phi,     R_p2_phi,     1,             R_phi_lambda_B,
#                  R_p1_lambda_B,R_p2_lambda_B,R_phi_lambda_B, 1), ncol = 4)
#   
#   ## vector with standard deviation of each variable ##
#   standard_dev <- c(s_p1, s_p2, s_phi, s_lambda_B)
#   D <- diag(standard_dev) # create a diagonal matrix 
#   
#   ## convert correlation matrix to variance-covariance matrix ##
#   VCV <- D%*%R%*%D 
#   
#   input <-list(N=N,mu=mu,VCV=VCV,cat=cat,theta=theta,lambda=M_lambda_1)
#   
#   
#   return(generate.data.heterogeneity(input))
#     
# } # end of test.generate.data.heterogeneity
# 
# test.generate.data.heterogeneity()
# 
# 
####################### end of testing the function ###########################
###############################################################################



