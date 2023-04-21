############################################################################### 
########### initial estimates for the data with individual covariates ######### 
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
###############################################################################
 
# The following function creates initial estimates for the optimization routine.

# Function "ic.initial.estimates"
# Input : Data
# Output : initial values for capture probabilities, category proportions, 
#          sub-sample proportions and population size

ic.initial.estimates = function(Data){ 
  cats = Data$category
  
  # n1 is the  total number of individuals caught in  sample 1
  # This is calculated counting the histories without a 0 in the first
  #  digit of the history (e.g. U0,UU, M0,MM, ...)
  n1 = sum(Data$counts[ ! (substr(Data$History,1,1)=="0" )])
  
  # m2 is the number of marked individuals caught in sample 2 ( e.g. UU, MM, FF,...)
  m2 = n1 - sum(Data$counts[ (substr(Data$History,2,2)=="0" )]) 
  
  # n2 is the  total number of individuals caught in  sample 2 
  # ( m2 + count the histories with a 0 in the first digit of the history)
  n2 = m2 + sum(Data$counts[ (substr(Data$History,1,1)=="0")])
  
  # n1.star is the number of individuals in sub-sample  in sample 1 
  # ( subtract total count of the  histories with a U in the first digit of
  #   the history from n1)
  n1.star = n1 - sum(Data$counts[substr(Data$History,1,1)=="U" ])
  
  
  # n2.star is the number of individuals in sub-sample  in sample 2 
  # ( subtract total count of the  histories with  0U from histories with
  #   a 0 in the first digit of the history )
  n2.star =  sum(Data$counts[ (substr(Data$History,1,1)=="0")]) -
                                    sum(Data$counts[Data$History=="0U"]) 
  
  
  #########################
  # count the number of individuals in each category in sub-sample 1 (ie.  n1.star)
  cat.count.time1 <- rep(0,length(cats))
  for(i in 1:(length(cats)-1)){
    cat.count.time1[i] <- sum(Data$counts[substr(Data$History,1,1)==cats[i] ])
  }
  cat.count.time1[length(cats)] <- n1.star- sum(cat.count.time1[1:(length(cats)-1)])
  
  # count the number of individuals in each category in sub-sample 2 (ie.  n2.star)
  cat.count.time2 <- rep(0,length(cats))
  for(i in 1:(length(cats)-1)){
    cat.count.time2[i] <- sum(Data$counts[(substr(Data$History,1,1)=="0") & (substr(Data$History,2,2)==cats[i])])
  }
  cat.count.time2[length(cats)] <- n2.star- sum(cat.count.time2[1:(length(cats)-1)])
  
  # count the number of individuals in each category in both sub-sample times
  cat.count <- cat.count.time1 + cat.count.time2
  
  
  ############ Initial values ###################
  # initial value for the population size; this is the simple Lincoln-Petersen estimate
  N.init = n1*n2/m2 

  theta_1 = n1.star/n1 # initial value for sub-sample proportion in sample 1
  theta_2 = n2.star/(n2-m2) # initial value for sub-sample proportion in sample 2
  
  # initial values for capture probabilities
  sample.1.cap.prob = rep(n1/N.init , length(cats))  # at time 1
  sample.2.cap.prob = rep(n2/N.init , length(cats))  # at time 2
  
  # initial values for category proportions:   
  #### use equal category proportions as initial values # this was the previous method
  lambda.init  = rep(1/length(cats),length(cats)) ## use this as the initial value
  #lambda.init  <- cat.count/(n1.star+n2.star)
  #lambda.init  <- c(0.3,0.7)
    
  # all initial capture probabilities as a vector
  cap.prob.init = c(sample.1.cap.prob,sample.2.cap.prob) 
  theta.init = c(theta_1 ,theta_2) # initial sub-sample proportions
  
  initial.estimates = NULL
  initial.estimates$p = matrix(cap.prob.init,ncol=2, byrow = FALSE)
  colnames(initial.estimates$p) = c("time1", "time2")
  rownames(initial.estimates$p) = c(paste(cats))
  
  initial.estimates$lambda = lambda.init
  names( initial.estimates$lambda)= c(paste("lambda_",cats))
  
  initial.estimates$theta = theta.init
  names(initial.estimates$theta) =  c("theta_1", "theta_2")
  
  
  initial.estimates$full = c(cap.prob.init,lambda.init,theta.init)
  
  
  return(initial.estimates) 
  
} # end of function "ic.initial.estimates"


###############################################################################
##########  Test Function #####################################################
###############################################################################
# ### test the function initial.estimates
#  
#  test.ic.initial.estimates = function(){
#   Data = NULL
#   Data$History   = c("U0", "U0", "M0" ,"MM", "F0" ,"F0" ,"FF", "0M",
#                      "0F", "0U", "0U", "0U","UU", "UU" )
#   Data$counts    = c(20, 21,  16 ,  4,  3, 20,  7, 12, 25, 43,  2, 18,  5,  4)
#   Data$category  = c("M","F") 
#   ##  here U0 = 41, UU= 9, M0 = 16, MM = 4, F0 = 23, FF= 7, 0M = 12, 
#   ##       0F= 25, 0U = 63
#   
#   return(ic.initial.estimates(Data))
#  }
#  test.ic.initial.estimates()

###############################################################################
########## End  Test Function #################################################
############################################################################### 