###############################################################################
###############################################################################
####### negative log likelihood for the data with individual covariates #######
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
###############################################################################   

# Function ic.neg.log.likelihood 
#  Input : parameters are in logit form, raw Data, indicator matrices,
#          Design Matrices and 
#          OFFSET vectors
#          (estimates are in the order: beta parameters acording to the capture formula,
#           lambda and theta)
#           (Data is a list with  capture history, individual covariate,
#             counts, category)
#  Output: negative log likelihood value


#logit.est=ic.initial.pack.parm.est  # delete this line later
ic.neg.log.likelihood <- function(logit.est, Data,indicator,
                                 captureDM,thetaDM,lambdaDM,
                                 captureOFFSET,thetaOFFSET,lambdaOFFSET){

  parm   <- ic.unpack.parm(logit.est,Data,
                          captureDM,thetaDM,lambdaDM,
                          captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  
  history <- Data$History    # capture history from raw data
  #counts  <- Data$counts     # counts from raw data 
                                   #no need this because we have data for individuals
  cats <- Data$category      # categories
  
  # call the function 'ic.log.prob.history'. These probabilities are 
  # in the exact order according to the capture histories. 
  # if phi.flag is "YES" then calculate phi hat (estimate for capture probability) 
  #   for each history and use it to calculate the N hat (estimate for population size)
  ic.log.prob = ic.log.prob.history(parm,history,cats,indicator,phi.flag="YES")  
  
  
  ## log-likelihood has 2 parts
  ##    1. conditional probability part for the captured individuals
  ##    2. probabilities for the not captured individuals
  #
  ## since probabilities for the not captured individuals cannot be calculated,
  ## no need to consider this
  
  part1 = sum(ic.log.prob$log_prop_hist) - sum(ic.log.prob$log_p_capture_i)
  
  ## in part1, the likelihood contribution is calculated as follows  
  ## the conditional likelihood
  ## sum( log(p_i/ (1- p00_i)) )= sum( log(p_i) ) - sum( log(1- p00_i) )
  
  
  ### log-likelihood ###
  LLH = part1
  
  # negative log-likelihood
  NLLH = - LLH  # -(log-likelihood),convert to negative value for optimization purpose
  #print(NLLH)
  return(NLLH)
  
} # end of "ic.neg.log.likelihood"

###############################################################################  
###############################################################################
############ compute the probability of each history passed to it. ############

ic.log.prob.history <- function(parm, history, cats, indicator, phi.flag="NO"){
  # This computes the log(probability) of each history based on the 'parm' list

        h <- indicator$h   # indicator variable to identify the capture history
                           ####  1 - captured,   0 - Not not captured
        s <- indicator$s   # stratified (sub-sampled) or not at the current sampling time.
                           ####  1 - stratified,  0 - Not stratified
       cs <- indicator$cs # captured in the sample or not at the current sampling time 
                          # from the animals not marked earlier.
                          ####  1 - captured in the sample ( i.e. "U", "M", "F" )
                          ####  0 - Not captured  in the sample ( "0")
    s.ind <- indicator$s.ind  # whether the captured animal is in particular category or "U"
                              #### column 1 indicates the captures animal is in
                              #### cats[1] or "U" ( ie. "M" or "U") by 1 and otherwise 0 
                              #### column 2 indicates the captures animal is in 
                              #### cats[2] or "U" ( ie. "F" or "U") by 1 and otherwise 0
                              ####
                              ####  and so on.....
  
  #create an array  for capture probability for each individual acording to
  #number of categories
  p <- array(parm$p,dim=c(length(history),nrow=length(cats), ncol=2))
  
  lambda<- parm$lambda
  theta <- parm$theta
  
  
  prob_hist<- rep(0, nrow(h)) # prob for each captured history
  prob.hist.00 <- rep(0, nrow(h)) # probability of not capturing the each individual
  for( j in 1:nrow(h) ){
    for(i in 1:length(lambda)){
      temp <- lambda[i]*prod( (p[j,i,]^(h[j,])) * 
                                (theta^(s[j,])) * 
                                ((1-p[j,i,])^(1-h[j,])) * 
                                ((1-theta)^((1-s[j,])*cs[j,])) )* s.ind[j,i]
      prob_hist[j] <- temp + prob_hist[j] 
    }
    
    # probabilities of  capture history "00" for each captured individual
    # i.e. probability that the each individual is not catchable
    prob.hist.00[j] <- sum(apply((1-p[j,,]),1,prod) *lambda)
  }
  

  # total probability of capture an individual(at least capture once)
  p_capture_i <- 1- prob.hist.00
  
  log_p_capture_i <- log(p_capture_i + (p_capture_i==0))
  
  # log probability of capture each individual
  log_prop_hist <- log(prob_hist + (prob_hist==0))
    
  if(phi.flag == "YES"){
      phi <- p_capture_i
  }
    
  out <- NULL
  out$log_prop_hist <- log_prop_hist
  out$log_p_capture_i <- log_p_capture_i
  
  if(phi.flag == "YES"){
    out$phi <- phi
  }
  
  return(out)
}  # end of ic.log.prob.history

#############################################
#############################################  New code - February 2016

###############################################################################



###############################################################################
###############################################################################
############### test the function  "test.ic.neg.log.likelihood" ###############
# 
# # data and patameter estimates for checking functions
test.ic.neg.log.likelihood  <- function(){
   
   Data <- NULL
   #including the "00" history
   Data$History <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U", "MM", "FF", "0M")
   Data$counts <- rep(1,12)
 ###  Data$length <- rnorm (12,18,4 )
   Data$length <- c(15.49, 11.25, 21.35, 18.61, 13.44, 23.01, 19.70, 16.81,
                    21.58, 21.51, 21.28, 20.75)
   Data$length <- Data$length - mean(Data$length)
   Data$category <- c("M","F")
   str(Data)
   print(Data)
   
   theta <- c(0.9,0.6) # sub sample proportions
   # category proportions. sum should be equal to 1 and the and the values 
   # corresponds to the categories defined above
   lambda <- c(0.7,0.3) 
   
   # give the capture formula
   captureformula <- ~length+ I(length^2) + category+ time
   captureformula <- ~length+ category+ time
   
   ### convert to logit scale  #####
   beta.logit.theta <- logit(theta)
   beta.logit.lambda <- lambda.to.cumulative.logit.lambda(lambda) 
   # give the beta parameters acording to the capture formula
   beta.logit.p <- c(0,2,-0.06,-0.3,-0.5) # for length length**2 category*time
   beta.logit.p <- c(0,2,-0.3,-0.5)  # for length + category*time
   
   # logit estimates are put in the same order as used in the other functions
   logit.est <- c( beta.logit.p,beta.logit.lambda,beta.logit.theta)
   
 
   # design MatriX FOR CAPTURE PROBABILITIES
   captureDM <- ic.create.DM(Data,captureformula)
 
   # design Matrices FOR THETA AND LAMBDA
   thetaDM   <- create.DM(c(1,2)) 
   lambdaDM  <- create.DM(c(1)) 
     
   #give the offset vectors(vectors of zero's should be given since no restriction)
   captureOFFSET <- matrix(0, length(Data$History), ncol= (length( Data$category) *2) )
   thetaOFFSET   <- c(0,0)
   lambdaOFFSET  <- c(0)
   
   indicator <- ic.create.indicator(Data)

   
   parm   = ic.unpack.parm(logit.est,Data,
                           captureDM,thetaDM,lambdaDM,
                           captureOFFSET,thetaOFFSET,lambdaOFFSET)
   print(parm)

   res1 <- ic.neg.log.likelihood(logit.est, Data,indicator,
                                captureDM,thetaDM,lambdaDM,
                                captureOFFSET,thetaOFFSET,lambdaOFFSET)
 
   print(paste("Negative log Likelihhod: " , res1  , sep = ""))
   
   # compare the hessian functions from optimHess and ic.inv.fisher.info
   hess <- ic.inv.fisher.info(logit.est,Data, indicator,
                              captureDM, thetaDM,lambdaDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET) 
   hess2 <-  optimHess(par=logit.est, fn=ic.neg.log.likelihood, 
               Data=Data, indicator=indicator, 
               captureDM =captureDM, thetaDM=thetaDM,lambdaDM=lambdaDM,
               captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
               lambdaOFFSET=lambdaOFFSET)
   browser()
   hess2
   solve(hess) # need to compare the appropriate matrices
   

   res2 = ic.log.prob.history(parm,Data$History,Data$category,indicator, phi.flag="YES")
   print(" probabilities for each of the histories")
   print(Data$History)
   print(exp(res2$log_prop_hist)) 
   
   # this is just to look at the plots for generated data with given parameters
   ic.plots(logit.est,Data, Data$length, 
            captureformula,thetaDM,lambdaDM,
            thetaOFFSET,lambdaOFFSET)
   
 }
# 
# test.ic.neg.log.likelihood()
# 

###############################################################################
################ check the function  "ic.log.prob.history" ####################

# check.ic.log.prob.history <- function(){
# 
#    Data <- NULL
#    # including the "00" history
#    Data$History <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U", "MM", "FF", "0M")
#    Data$category <- c("M","F")
#    str(Data)
#    print(Data)
# 
#    # Set up some values of the parameters for testing
#    parm <- NULL
#    parm$p <- matrix(c(0.1,0.2,0.3,0.4,
#                       0.2,0.3,0.1,0.3,
#                       0.1,0.4,0.3,0.2,
#                       0.3,0.3,0.2,0.2,
#                       0.4,0.2,0.3,0.1,
#                       0.2,0.1,0.1,0.2,
#                       0.2,0.2,0.3,0.1,
#                       0.4,0.2,0.3,0.1,
#                       0.2,0.3,0.3,0.2,
#                       0.3,0.3,0.2,0.2,
#                       0.3,0.1,0.1,0.2,
#                       0.1,0.1,0.1,0.2),ncol=4, byrow=TRUE)
#    parm$lambda <- c(0.4, 0.6)
#    parm$theta <- c(0.8, 0.6)
#    print(parm)
#    indicator <- ic.create.indicator(Data)
#    # Test out on all possible capture histories (not including the 00 history)
#    res <- ic.log.prob.history(parm, Data$History, Data$category,indicator,phi.flag="YES")
#    print(res)
#    print(" probabilities for each of the histories")
#    print(Data$History)
#    print(exp(res$log_prop_hist)) 
#    print(exp(res$log_p_capture_i))
# }
# 
# check.ic.log.prob.history()

###############################################################################
###############################################################################

# # ### following are just to ckeck the  probabilities 
# 
# pm1 = parm$p[,1]
# pf1 = parm$p[,3]
# pm2 = parm$p[,2]
# pf2 = parm$p[,4]
# 
# lambda_1 = 0.4
# lambda_2 = 0.6 
# 
# theta_1 = 0.8
# theta_2 = 0.6 
# 
# p00 = lambda_1*(1-pm1)*(1-pm2) + lambda_2* (1-pf1)*(1-pf2)
# 
# #### conditional probabilities
# 
# pMM= lambda_1*pm1*theta_1*pm2
# pFF= lambda_2*pf1*theta_1*pf2
# 
# pM0 = lambda_1*pm1*theta_1*(1-pm2)
# pF0 = lambda_2*pf1*theta_1*(1-pf2)
# 
# p0M = lambda_1*(1-pm1)*pm2*theta_2
# p0F = lambda_2*(1-pf1)*pf2*theta_2
# 
# pUU = (lambda_1*pm1*(1-theta_1)*pm2 +lambda_2*pf1*(1-theta_1)*pf2)
# p0U = (lambda_1*(1-pm1)*(1-theta_2)*pm2 +lambda_2*(1-pf1)*(1-theta_2)*pf2)
# pU0 = (lambda_1*pm1*(1-theta_1)*(1-pm2) +lambda_2*pf1*(1-theta_1)*(1-pf2))
# 
# 
# total = pMM +  pFF + pM0+ pF0+ p0M+ p0F+ pUU+ p0U+ pU0+ p00
# 
# conditional_tot = (pMM +  pFF + pM0+ pF0+ p0M+ p0F+ pUU+ p0U+ pU0)/(1-p00)
# 


