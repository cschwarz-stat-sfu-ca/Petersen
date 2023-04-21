###############################################################################
#####               calculate negative log likelihood for                ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################
# FunctionPSkSCR.neg.log.likelihood 
#  Input : parameters are in logit or log form, raw Data, Design Matrices and 
#          OFFSET vectors
#          (estimates are in the order:logit.capture.rates(p),
#            logit.cumulative.category.proportion(lambda), 
#            logit.sampling.fraction(theta), log.population.size(N))
#           (Data is a list with  capture history, counts, category)
#  Output: negative log likelihood value

PSkSCR.neg.log.likelihood <- function(logit.est, data,indicator,
                                      captureDM,thetaDM,lambdaDM,p_lossDM,
                                      captureOFFSET,thetaOFFSET,
                                      lambdaOFFSET,p_lossOFFSET){
  
  parm <- PSkSCR.unpack.parm(logit.est,data,
                             captureDM,thetaDM,lambdaDM,p_lossDM,
                             captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
  
  N <- parm$N
    
  hist <- data$history # capture history from raw data
  abs.ct   <- abs(data$counts)   # absolute values of counts from raw data
  ct.00 <- (N-sum(abs.ct)) # count for capture history "00"
  # counts of observable capture histories and capture history "00"
  ct <- c(abs.ct,ct.00)
  
  
  # call the function 'PSkSCR.log.prob.history'. These probabilities and the counts
  # are in the exact order related to each observable capture history.
  # probabilities of observable capture histories and capture history of unseen animals
  # (i.e. for two sampling occasions "00", for 3 sampling occasion "000" , ..)
  log.prob.hist <- PSkSCR.log.prob.history(parm,indicator) 
  
  #######################################
  # log-ikelihood has 2 parts
  part1 <- lgamma(N+1) - sum(lgamma(ct+1)) # factorial part
  
  part2 <- sum(ct * log.prob.hist)
  
  # log-likelihood
  LLH <-  part1 +  part2  
  
  # negative log-likelihood
  NLLH <- - LLH  # -(log-likelihood),convert to negative value for optimization purpose
  
  return(NLLH)
  
} # end of "PSkSCR.neg.log.likelihood "



###############################################################################
####### functions to compute the probability of each history passed to it. ####
############################################################################### 

# histories are passed through indicator variables (matrices)

PSkSCR.log.prob.history <- function(parm,indicator){
  
        h <- indicator$h
        z <- indicator$z
  cat.ind <- indicator$cat.ind
        s <- indicator$s
       cs <- indicator$cs
    s.ind <- indicator$s.ind
     dead <- indicator$dead
  
       p <- parm$p
  lambda <- parm$lambda
   theta <- parm$theta
  p_loss <- parm$p_loss
  
  
  
  prob_hist<- rep(0, nrow(h))
  for( j in 1:nrow(h) ){
    for(i in 1:length(lambda)){
      temp <- lambda[i]*prod( (p[i,]^(h[j,]*z[j,])) * 
                                ((1-p_loss)^(h[j,]*dead[j,])) * 
                                (theta^(s[j,]))* 
                                ((1-p[i,])^((1-h[j,])*z[j,])) * 
                                (p_loss^(h[j,]*(1-dead[j,]))) *
                                ((1-theta)^((1-s[j,])*cs[j,])) ) * s.ind[j,i]
      prob_hist[j] <- temp + prob_hist[j] 
    }
  }
  
    
  # probability of observed capture histories
  prob.hist.observed <- prob_hist
    
  # probability for capture history of unseen animals
  # (i.e. for two sampling occasions "00", for 3 sampling occasion "000" , ..)
  prob.hist.00 <- sum(apply((1-p),1,prod) *lambda)
  
  # probabilities of observable capture histories and capture history "00"
  prob.hist <- c(prob.hist.observed,prob.hist.00)
  
  if(sum(prob.hist)>1){
    print(" Error in probability of histories calculations")
     stop
  }
  
  # log- probabilities of observable capture histories and capture history "00"
  log_prob.hist <- log(prob.hist + (prob.hist==0))
  
  return(log_prob.hist)
  
  
} # end of function "PSnSCR.log.prob.history"


###############################################################################
###############################################################################
############ test the  function PSkSCR.log.prob.history #######################
# 
# test.PSkSCR.log.prob.history <- function(){
# 
#    data <- NULL
#    data$history <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U", "M0","FF","0F")
#    data$counts   <- c(41,9,16,4,23,7,12,25,63,-2,-1,-3)
#    data$category <- c("M","F")
#    str(data)
#    print(data)
#    
#    indicator <-PSkSCR.create.indicator(data)
#    
#    # Set up some values of the parameters for testing
#    parm <- NULL
#    parm$p <- matrix(seq(.1,.4,.1),nrow=2)
#    parm$lambda <- c(0.4, 0.6)
#    parm$theta <- c(0.3, 0.6)
#    parm$N <- 500
#    print(parm)
# 
#    # Test out on all observed capture histories 
#    res <-PSkSCR.log.prob.history(parm, indicator)
#    print(res)
# }
# 
# test.PSkSCR.log.prob.history()

###############################################################################
############ test the  function PSkSCR.neg.log.likelihood #####################

# test.PSkSCR.neg.log.likelihood <- function(){
#   
#     data <- NULL
#     data$history <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U", "M0","FF","0F")
#     data$counts   <- c(41,9,16,4,23,7,12,25,63,-2,-1,-3)
#     data$category <- c("M","F")
#     str(data)
#     print(data)
#       
#     # Set up some values of the parameters for testing
#     parm <- NULL
#     parm$p <- matrix(seq(.1,.4,.1),nrow=2)
#     parm$lambda <- c(0.4, 0.6)
#     parm$theta <- c(0.3, 0.6)
#     parm$N <- 500
#     print(parm)
#          
#      # give  the required design matrices for capture recapture probabilities,
#     # theta(sampling(sexing) fractions ) and lambda(Category proportion)
#     captureDM = create.DM(c(1,2,3,4)) 
#     thetaDM   = create.DM(c(1,2)) 
#     lambdaDM  = create.DM(c(1)) 
#          
#     #give the offset vectors(vectors of zero's should be given since no restriction)
#     captureOFFSET = c(0,0,0,0) 
#     thetaOFFSET   = c(0,0)
#     lambdaOFFSET  = c(0)
#         
#     initial.est = PSkSCR.initial.estimates(data)
#     pack.init.est = PSkSCR.pack.parm(initial.est$full,data,captureDM,thetaDM,lambdaDM)
#     
#     indicator <-PSkSCR.create.indicator(data)
#     
#     Neg.Log.like <- PSkSCR.neg.log.likelihood(pack.init.est, data,indicator,
#                                               captureDM,thetaDM,lambdaDM,
#                                               captureOFFSET,thetaOFFSET,lambdaOFFSET)
#     
#      names(Neg.Log.like) <- c("Neg.Log.Like")
#     print(Neg.Log.like)
# 
# } # End of test.PSkSCR.neg.log.likelihood
# 
# test.PSkSCR.neg.log.likelihood()


############ End of  test the  function PSkSCR.neg.log.likelihood #############
###############################################################################

  