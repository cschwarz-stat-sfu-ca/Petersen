###############################################################################
###############################################################################
# Function neg.log.likelihood 
#  Input : parameters are in logit or log form, raw Data, Desigm Matrices and 
#          OFFSET vectors
#          (estimates are in the order:logit.capture.rates(p),
#            logit.cumulative.category.proportion(lambda), 
#            logit.sampling.fraction(theta), log.population.size(N))
#           (Data is a list with  capture history, counts, category)
#  Output: negative log likelihood velue

neg.log.likelihood = function(logit.est, Data,
                              captureDM,thetaDM,lambdaDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET){
  
  parm   = unpack.parm(logit.est,Data,
                       captureDM,thetaDM,lambdaDM,
                       captureOFFSET,thetaOFFSET,lambdaOFFSET)
   
  # N hat
  N = parm$N
  
  hist = Data$History       # capture history from raw data
  ct   = Data$counts        # counts from raw data
  cats = Data$category      # categories
    
  History = c(hist,"00")      # add "00" capture history to end of the
                              # vector of capture history from raw data
      
  counts  = c(ct,(N-sum(ct))) # add total animals did not catch(relating to 
                              # capture history "00") to the end of 
                              # the vector of counts from the raw data

  # call the function 'log.prob.history'. These probabilities and the counts
  # are in the exact order related to each capture history. 
  log.prob = log.prob.history(parm,History,cats)  
    
  
  # log-ikelihood has 2 parts
  part1 = lgamma(N+1) - sum(lgamma(counts+1)) # factorial part
  
  part2 = sum(counts * log.prob)
  
  # log-likelihhod
  LLH = part1 + part2  
  
  # negative log-likelihood
  NLLH = - LLH  # -(log-likelihood),convert to negative value for optimization purpose
  
  return(NLLH)
  
} # end of "neg.log.likelihood"

###############################################################################  
###############################################################################
##### CJS functions to compute the probability of each history passed to it. ##

log.prob.history <- function(parm, history, cats){
# This computes the log(probability) of each history based on the parms list
# of parameters and the vector of category labels.
#
# We assume that parms has been unpacked and has elements (at a minimum)
#  p - length(cat) x 2 matrix
#  lambda - proportions of categories in population (must add to 1)
#  theta  - probability that an unmarked fish captured is “categorized”
#
# cats is a vector of category codes. The value of “U” CANNOT be used as a valid
#   code as this indicates an unknown category. Upper and lower case categories
#   are ok, but the user must use caution.
#
# history - a vector of histories. Each history is a character string of length 2
#  e.g. M0 indicates a male captured at time 1 and never seen again
# 
# This routine will compute the probability of history (00) if you pass this to this
# routine even though this history is unobservable.
#
# Returns a vector of the log(prob) of each history. If a history has probability of 0
# a 0 is returned rather than a -Inf or NA. [This will allow graceful handling
# of cases when certain parameters are fixed at 1 or 0

 
# We construct a (length(cat)+2) x (length(cat)+2) matrix representing all of the
# possible capture histories log(probabilites) for (cat, "0" ,”U”) combinations
  prob <- matrix(0, nrow=length(cats)+2, ncol=length(cats)+2)

# Categories are captured at both occasions (e.g. MM)
  prob[matrix(rep(1:length(cats),2),ncol=2)] <- parm$lambda*parm$p[,1]*parm$theta[1]*parm$p[,2]

# Categories captured at first occasion and then never seen (e.g M0)
  prob[matrix(c(1:length(cats),
                rep(length(cats)+1,
                    length(cats))),ncol=2)] <- parm$lambda*parm$p[,1]*parm$theta[1]*(1-parm$p[,2])

# Categories not captured at first occassion, but captured at time 2 (e.g. 0M)
  prob[matrix(c(rep(length(cats)+1,
                    length(cats)),1:length(cats)),
              ncol=2)] <- parm$lambda*(1-parm$p[,1])*parm$p[,2]*parm$theta[2]
  
# History (00)
  prob[length(cats)+1,length(cats)+1] <- sum(parm$lambda*(1-parm$p[,1])*(1-parm$p[,2]))
  
# History (U0)
  prob[length(cats)+2,
       length(cats)+1] <- sum(parm$lambda*parm$p[,1]*(1-parm$theta[1])*(1-parm$p[,2]))
  
# History (0U)
  prob[length(cats)+1,
       length(cats)+2] <- sum(parm$lambda*(1-parm$p[,1])*parm$p[,2]*(1-parm$theta[2]))
  
# History (UU)
  prob[length(cats)+2,
       length(cats)+2] <- sum(parm$lambda*parm$p[,1]*(1-parm$theta[1])*parm$p[,2])

  if(sum(prob)<.9999999){stop("Error in computing cell probabilities")}  

# Now index each history into the matrix of probabilities and return the log(prob) 
  index <- matrix(c(match(substr(history,1,1),c(cats,"0","U")),
                    match(substr(history,2,2),c(cats,"0","U"))),ncol=2)
  #print(index)
  log_p_hist <- log(prob[index]+(prob[index]==0))
  #print(log_p_hist)

  return(log_p_hist)
}  # end of log.prob.history


###############################################################################
###############################################################################
# #Following code is to check the above functions
# 
# # data and patameter estimates for checking functions
# testing.data <- function(){
#   
#   data <- NULL
#   data$History  = c("U0","UU","M0","MM","F0","FF","0M","0F","0U")
#   data$counts    = c(41,9,16,4,23,7,12,25,63)
#   data$category = c("M","F")
#   str(data)
#   #print(data)
#   
#   # Set up some values of the parameters for testing
#   parm = NULL
#   parm$p      = matrix(seq(.1,.4,.1),nrow=2)
#   parm$lambda = c(0.4, 0.6)
#   parm$theta  = c(0.3, 0.6)
#   parm$N      = 500
#   str(parm)
#   #print(parm)
#   
#   return(list(parm,data))
# }
# 
# 
# # test  both 'neg.log.like' and 'log.prob.history' functions
# # 'log.prob.history'should produce sum of the probabilities equal to 1
# 
# check.both.functions  = function(){
#   # Check the log.prob.history function
#   parm = testing.data()[[1]]
#   data = testing.data()[[2]]
#   
#   res1  = neg.lg.ld(parm,data)
#   print(paste("Negative log Likelihhod: " , res1 , sep = ""))
#   
#   res2 = log.prob.history(parm,c(data$History,"00"),data$category)
#   print(paste("sum of the probabilities: " , sum(exp(res2)), sep = "")) # to see add to 1
# }
# 
# check.both.functions()
###############################################################################

###############################################################################
# ##### CJS function to check probability #####################################
# 
# check.log.prob.history <- function(){
# # Check the log.prob.history function
# 
#    data <- NULL
#    data$History <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U")
#    data$Count   <- c(41,9,16,4,23,7,12,25,63)
#    data$Category <- c("M","F")
#    str(data)
#    print(data)
# 
#    # Set up some values of the parameters for testing
#    parm <- NULL
#    parm$p <- matrix(seq(.1,.4,.1),nrow=2)
#    parm$lambda <- c(0.4, 0.6)
#    parm$theta <- c(0.3, 0.6)
#    print(parm)
# 
#    # Test out on all possible capture histories (including the 00 history)
#    res <- log.prob.history(parm, c(data$History,"00"), data$Category)
#    sum(exp(res))  # see if add to 1
# }
# 
# check.log.prob.history()
# 
# 
###############################################################################
###############################################################################
