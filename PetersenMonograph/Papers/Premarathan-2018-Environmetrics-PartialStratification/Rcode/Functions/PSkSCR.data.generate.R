
###############################################################################
########## Generate data for Partial Stratification in  #######################
##########    k-Sample Capture-Recapture Experiments    #######################
##########             with known dead removals         ####################### 
###############################################################################

## function "PSkSCR.data.generate" generates data for each category
# if the count is a gegative value, that means number of losses on capture

# st= number of sampling occasions
# N = population size
# theta =sub sample proportions
# category = give the category  ie. "M"
# lambda = category proportion for the given category
# p = capture probabilities for each sampling time for given category
# p_loss =  probability that animal is dead at the sampling occasion  j wherej = 1, ... , t
#           this is used to get the know removals in each sampling time


PSkSCR.data.generate <- function(st, N, theta, category, lambda, p, p_loss){
  
  # number of histories to be generated
  size <- N*lambda
  #initialize an individual history
  ind.history <- do.call(paste, c(as.list(rep(0,st)), sep=""))
  # initialize histories eg.  as "000" for 3 capture occasions..
  history <- rep(ind.history, size)
  
  # initialize the vector "update_cap". this keeps whether the animal is 
  # captured or not during  the experiment, 1-captured 0- not
  update_cap <- rep(0,size)
  # initialize the vector "update_sex". this keeps whether the animal is 
  # categoried(sub-sampled) or not during  the experiment, 1- sub-sampled 0-not
  update_sex <- rep(0,size)
  # initialize the vector "update_available". this keeps whether the animal is 
  # available to capture or not available during  the experiment, 
  # animal may be alive or dead. 1-available 0- not available
  update_available <- rep(1,size)
  
  for(i in 1:st){
    
    # generate a vector indicating the capture status(0/1 from bernoulli process)
    cap <- rbinom(size,1,p[i]) # 1-capture  , 0-not capture
    
    # generate a vector indicating the category status(0/1 from bernoulli process)
    sex <- rbinom(size,1,theta[i]) # 1-stratify  , 0-not stratify
    
    # generate a vector indicating whether the individual is alive or 
    # dead (0/1 from bernoulli process)
    alive <- rbinom(size,1,(1-p_loss[i])) # 1-alive  , 0-dead
    
    # saub-sampled in the current sampling occasion
    sex_currect <- as.numeric((sex+cap) > 1)
    
    
    # 0-not captured, "U"- captured but not stratified, "cat"- stratified (eg:"M")
    # assign "U" for all captured animals
    substr(history[update_available & cap ],i,i) <- "U" 
    # if subsampled earlier, then category is known
    substr(history[update_available & update_sex & cap ],i,i) <- category
    # sub-samples the animals those are NOT captured earlier
    substr(history[update_available & (!update_cap) & cap & sex ],i,i) <- category 
    
    
    ## even the animal is dead, if it has not been captured at the current time, 
    ## it is possible to capture at the next time.
    temp1 <- as.numeric(( (!cap) | alive))
    # then the availables animals for next time are as follows
    update_available <- as.numeric((update_available + temp1) > 1)
    
    update_sex <- as.numeric((!update_cap) & sex_currect) + update_sex
    update_cap <- as.numeric((update_cap + cap) > 0)
        
  }
  
  # variable "cap.ind"  shows as follows
  ##  cap.ind = 0 means not captured
  ##  cap.ind = -1 means captured a dead animal ( can be any sampling time t)
  ##  cap.ind = 1 means captured a alive animal
  
  cap.ind <- rep(1,size)  # initialize all with 1
  cap.ind[which(update_available==0)]<- -1 # indicate the dead animals by -1
  cap.ind[ which(history==ind.history)] <-  0  # assign 0's for not captured animals
  
  # create a two-way table with histories and cap.ind
  summary <- table(history,cap.ind) 
  
  data.summary <- as.data.frame(summary)
  names(data.summary) <- c("history", "cap.ind", "frequency")
  # remove the histories "000"
  obs.summary <- data.summary[data.summary$history!=ind.history,] 
  
  new.cap.ind <- rep(0, nrow(obs.summary))
  
  new.cap.ind[ obs.summary$cap.ind==-1] <- -1
  new.cap.ind[ obs.summary$cap.ind!=-1] <- 1
  
  obs.summary$counts <- new.cap.ind * obs.summary$frequency
  
  # get the observed histories 
  obs.summary <- obs.summary[obs.summary$counts != 0, ]
    
  # only the observed histories and their counts
  obs.data.summary <- NULL
  obs.data.summary$history <- as.character(obs.summary$history)
  obs.data.summary$counts <- obs.summary$counts
  
  df.obs.data.summary <- as.data.frame(obs.data.summary)
      
  return(df.obs.data.summary)
  
} # end of function "PSkSCR.data.generate"

###############################################################################
###########  test the function  PSkSCR.data.generate   ########################
# 
#  test.PSkSCR.data.generate <- function(){
#   N <-100 # population size
#   st <- 3 # number of sampling occasions
#   theta <- c(0.4,0.3, 0.3) # sub sample proportions
#   p_loss <- c(0, 0.2, 0.1)
#   category <- "M" # category
#   lambda <- 0.7 # category proportion for males
#   p <- c( 0.2, 0.3, 0.4) # capture probabilities for each sampling time #
#   
#   out <- PSkSCR.data.generate(st, N, theta, category, lambda, p, p_loss)
#   
#   return(out)
#  }
# 
#  test.PSkSCR.data.generate()

###################   End of test function    #################################
###############################################################################