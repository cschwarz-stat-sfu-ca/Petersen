###############################################################################
#####      initial estimates for the optimization routine for            ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ###### 
############################################################################### 


# The following function creates initial estimates for the optimization routine.

# Function "PSkSCR.initial.estimates"
# Input : data
# Output : initial values for capture probabilities, category proportions, 
#          sub-sample proportions and population size

PSkSCR.initial.estimates <- function(data, indicator){
  
  cats = data$category
  
  ##number of sampling occasions
  st <-  nchar(data$history[1])
  
  # total number of animals captures at each sampling time
  n_t <- rep(0,st)
  
  # total number of animals with known stratification at each sampling time
  n_ss <- rep(0,st)
  
  for( i in 1:st){
    n_t[i] <- sum( abs( data$counts[ ! (substr(data$history,i,i)=="0" )]) )
    n_ss[i] <- n_t[i] - sum( abs( data$counts[substr(data$history,i,i)=="U"] ) )
  }
  
  ######################
  # initial sub-sample proportions ( this is a rough estimate) 
  theta.init = c(n_ss/n_t)
  names(theta.init) <-  c(paste("theta_",c(1:st), sep=""))
  ######################
  
    
  # n_t[1] is the  total number of individuals caught in  sample 1
  # n_t[2] is the  total number of individuals caught in  sample 2 
   
  # m2 is the number of marked individuals caught in sample 2 ( e.g. UU, MM, FF,...)
  m2 <- n_t[1] - ( sum(abs(data$counts[(substr(data$history,2,2)=="0" )])) -
                sum(abs(data$counts[(substr(data$history,1,2)=="00")])) )
  
  ######################
  # initial value for the population size;  This is calculated using the firt 
  # two sampling occasion and then use the  simple Lincoln Petersen estimate
  N.init <- n_t[1] *n_t[2] /m2 
  names(N.init) <- c("N")
  ######################
  
  ######################
  # initial values for capture probabilities; consider equal capture 
  # probability for all the categories in all capture occasions 
  # all initial capture probabilities as a vector
  cap.prob.init <- rep(n_t[1]/ N.init, st* length(data$category) )
  names(cap.prob.init) <- c(paste(rep(cats, st),"_Time_",
                                  rep(c(1:st),each=length(cats)), sep=""))
  ######################
  
  ######################
  # initial values for category proportions:   
  #                 use equal category proportions as inital values 
  lambda.init  <- rep(1/length(cats),length(cats))
  names( lambda.init) <- c(paste("lambda_",cats,sep=""))
  ######################
  
  ######################
  # initial values for loss on capture probabilities in each sampling time:  
  
  ### indicator variable "h" is to identify the capture history
  ###  1 - captured
  ###  0 - Not not captured
  h <- indicator$h
  
  
  ### indicator variable "dead" is to identify whether the 
  ### captured animal is alive or dead
  ###  1 - animal is alive
  ###  0 - animal is dead
  dead <- indicator$dead
  
  total.captured <- rep(0,st) # total captured at each sampling occasion
  total.dead <- rep(0,st) # total dead at each sampling occasion
  
  
  for( i in 1:st){
    total.captured[i] <- sum( (abs(data$counts))  * h[,i])
    temp <- sum( (abs(data$counts))  * h[,i] * dead[,i]) # total not dead
    total.dead[i] <- (total.captured[i] - temp)
  }
  
  # estimate for Proportion of dead animals captured at each sampling occasion 
  p_loss.init <- total.dead/total.captured
  names(p_loss.init) <-  c(paste("loss_Time_",c(1:st), sep=""))
  
  ######################
 
  
  initial.estimates <- NULL
  initial.estimates$p <- matrix(cap.prob.init,ncol=st, byrow = FALSE)
  colnames(initial.estimates$p) <- c(paste("Time_",c(1:st), sep=""))
  rownames(initial.estimates$p) <- c(paste(cats))
  
  initial.estimates$lambda <- lambda.init
  
  initial.estimates$theta = theta.init
  
  initial.estimates$p_loss <- p_loss.init
    
  initial.estimates$N <- N.init
  
   
  initial.estimates$full <- c(cap.prob.init,lambda.init,theta.init,p_loss.init,N.init)
  
  
  return(initial.estimates) 
  
} # end of function "PSkSCR.initial.estimates"


###############################################################################
##########  Test Function #####################################################
###############################################################################
# # test the function initial.estimates
#  
#  test.PSkSCR.initial.estimates = function(){
#   data = NULL
#   data$history   = c("U0", "U0", "M0" ,"MM", "F0" ,"F0" ,"FF", "0M",
#                      "0F", "0U", "0U", "0U","UU", "UU" )
#   data$counts    = c(20, 21,  16 ,  4,  -3, 20,  7, 12, 25, 43,  -2, 18,  5,  -4)
#   data$category  = c("M","F") 
#   ##  here U0 = 41, UU= 9, M0 = 16, MM = 4, F0 = 23, FF= 7, 0M = 12, 
#   ##       0F= 25, 0U = 63
# 
#   indicator <- PSkSCR.create.indicator(data)
#   
#   return(PSkSCR.initial.estimates(data,indicator))
#  }
#  test.PSkSCR.initial.estimates()
#
###############################################################################
########## End  Test Function #################################################
############################################################################### 
