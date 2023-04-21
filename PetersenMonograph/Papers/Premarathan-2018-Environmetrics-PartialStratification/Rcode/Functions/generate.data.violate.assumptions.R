#####################################################################################################
###### Identify the violation of assumption through model assesment #################################
#####################################################################################################
# generate some data where the assumptions of the model are violated and see if the residual plot and 
# the p-value of the gof test assesment capture the problem
# ( violation of assumptions : falier of non-death assumption, entering animals to the sampled population
# between the first and second period, etc.. )

#  phi = probability that animal alive at the time of the first sampling occation is still alive and 
#        present in the population at the time of the second sampling occation
#
#  B = number of animals entering the sampled population between the first and second perionds

# use expected value method.

generate.data.violate.assumptions = function(categories,cap.prob,lambda,theta,N, phi, B){
  
  given.param = NULL
  given.param$p = matrix(cap.prob, ncol=2, byrow = FALSE)
  given.param$lambda = lambda
  given.param$theta  = theta
  given.param$N      = N
  
  cap.histories = c(paste(categories,"0",sep=""),paste(categories,categories,sep=""),paste("0",categories,sep=""))
  cap.hist.1 = paste(categories,"0",sep="")  # Categories captured at first occasion and then never seen (e.g M0)
  cap.hist.2 = paste(categories,categories,sep="") # Categories are captured at both occasions (e.g. MM)
  cap.hist.3 = paste("0",categories,sep="") # Categories not captured at first occassion, but captured at time 2 (e.g. 0M)
  
  gen.Data = NULL  # generate data
  gen.Data$History = c("U0",cap.hist.1, cap.hist.2,"UU", cap.hist.3,"0U") # all possible observable capture histories("00" is unobservable)
  gen.Data$category = categories
  
  hist = c(gen.Data$History,"00")# add "00" capture history to end of the vector of capture history of expected.data
  
  # the function "log.prob.history" is in the file "neg.log.likelihood" which was used to find the log probabilities
  gen.cats = gen.Data$category
  expected.log.prob = log.prob.history(given.param,hist,gen.cats)  # call the function 'log.prob.history' in neg.log.likelihhod 
  # These log probabilities and the counts are in the exact order related to each capture history. 
  
  expected.prob = exp(expected.log.prob) - (expected.log.prob==0) # expected probabilities
  

  expected.counts.all = N * expected.prob # expected counts for all the capture histories including history "00"
  gen.Data$counts = expected.counts.all[1:length(expected.counts.all)-1] # expected couts for all the capture histories without history "00"   
  
  ncats = length(categories)
  expected.count.1 = N * expected.prob[1:(ncats+1)]  # histories with captured in the first occasion an newer seen in the second time
  expected.count.2 = N*phi*expected.prob[(ncats+2):(2*ncats+2)]
  expected.count.3 = (N*phi + B)*expected.prob[(2*ncats+3):(3*ncats+3)]
    
  gen.Data$counts = c(expected.count.1,expected.count.2,expected.count.3)
  
  return(gen.Data)
  
} 


####################################################################################################  
#################### Testing #######################################################################  
#   
# test = function(){
#   # give the  categories in the population
#   categories = c("M","F")
#   
#   # give the different  capture probabilities for p1M, p1F, p2M and p2F
#   cap.prob = c( 0.07, 0.05, 0.08, 0.1)
#   
#   # give the categoty proportions( total should add up to 1)
#   lambda= c(0.6, 0.4)
#   
#   # give the subsample proportions for the time 1 and 2
#   theta = c(0.6, 0.5)
#   
#   # give the population size
#   N = 200000
#   
#   phi = 0.8
#   B = 3500
# 
# 
#   #model identification :All capture probabilities are equal
#   model.id = paste("( p(c*t), theta(t), lambda(c) )")
#   
#   # give  the required design matrices
#   captureDM = create.DM(c(1,2,3,4)) # Design matrix for capture recapture probabilities (with equal probabilities)
#   thetaDM   = create.DM(c(1,2)) # Design matrix for theta(sampling(sexing) fractions)
#   lambdaDM  = create.DM(c(1)) # Design matrix for lambda(Category proportion)
#   
#   #give the offset vectors(vectors of zeros should be given since no restriction)
#   captureOFFSET = c(0,0,0,0) 
#   thetaOFFSET   = c(0,0)
#   lambdaOFFSET  = c(0)
#   
#   generate.data = generate.data.violate.assumptions(categories,cap.prob,lambda,theta,N, phi, B)
#   
#   model.gen.data_2 =  fit.model(model.id,generate.data,captureDM,thetaDM,lambdaDM,captureOFFSET,thetaOFFSET,lambdaOFFSET)
#   return(print.output(model.gen.data_2))
# }
# 
# test()
# 
# ####################################################################################################  
####################################################################################################  








1-pchisq(119.3,90)