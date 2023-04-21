#####################################################################################################
###### Identify the violation of assumption through the gof plots  ##################################
#####################################################################################################
# generate some data where the probabilities vary between categories and/or sample times and fit a model 
# where they are forced to equal.  Show the residuals and gof plots and see what features of these plots
# show extent of violation of assumptions.

generate.data = function(categories,cap.prob,lambda,theta,sample.size){
  
  given.param = NULL
  given.param$p = matrix(cap.prob, ncol=2, byrow = FALSE)
  given.param$lambda = lambda
  given.param$theta  = theta
  
  cap.histories = c(paste(categories,"0",sep=""),paste(categories,categories,sep=""),paste("0",categories,sep=""))
  
  gen.Data = NULL  # generate data
  gen.Data$History = c("U0", "UU", cap.histories,"0U") # all possible observable capture histories("00" is unobservable)
  gen.Data$category = categories
  
  hist = c(gen.Data$History,"00")# add "00" capture history to end of the vector of capture history of expected.data
  
  # the function "log.prob.history" is in the file "neg.log.likelihood" which was used to find the log probabilities
  gen.cats = gen.Data$category
  expected.log.prob = log.prob.history(given.param,hist,gen.cats)  # call the function 'log.prob.history' in neg.log.likelihhod 
  # These log probabilities and the counts are in the exact order related to each capture history. 
  
  expected.prob = exp(expected.log.prob) - (expected.log.prob==0) # expected probabilities
  
  
  obs.hist.exp.prob  =  expected.prob[1:(length(expected.prob)-1)] # expected probabilities for observable capture historis( without history"00")
  p_00 = expected.prob[length(expected.prob)]  # probability of the capture history "00"
  cond.obs.hist.exp.prob = obs.hist.exp.prob/p_00  # conditional probabilities for  observable capture histories
  
  gen.Data$counts = rmultinom(1, size = sample.size, prob = cond.obs.hist.exp.prob) # generate data 
  
  return(gen.Data)
  
}



##################### testing the function ######################################

# test.generate.data = function(){
# 
#   # give the  categories in the population
#   categories = c("M","F")
# 
#   # give the different  capture probabilities for p1M, p1F, p2M and p2F
#   cap.prob = c( 0.08, 0.06, 0.07, 0.06)
# 
#   # give the categoty proportions( total should add up to 1)
#   lambda= c(0.6, 0.4)
# 
#   # give the subsample proportions for the time 1 and 2
#   theta = c(0.8, 0.5)
# 
#   # total numbers individuals capture for the study 
#   sample.size = 4000
# 
#   res = generate.data(categories,cap.prob,lambda,theta,sample.size)
#   
#   return(res)
# }
# 
# test.generate.data()
#
###################### end testing ##############################################
