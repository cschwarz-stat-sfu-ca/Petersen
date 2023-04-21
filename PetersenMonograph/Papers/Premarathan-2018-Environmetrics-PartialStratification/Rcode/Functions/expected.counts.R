###############################################################################
############################################################################### 
# "expected.counts" generate the expected counts give the parameter set.
# Input : capture probabilities, category proportions,sub-sample proportions 
#         and categotires in the data set
# output : expected counts for the all possible observable capture histories
#           (except the capture history "00"- animals not captured)
#        : variance for each expected count
expected.counts = function(result){
 
   param  = result$param
   cats  = result$category 
  
   expected.Data = NULL
   # Create capture histories for Categories captured at first occasion and then
   # never seen (e.g. M0)and  Categories are captured at both occasions (e.g. MM) and 
   # Categories not captured at first occasion, but captured at time 2 (e.g. 0M)
   cap.hist = c(paste(cats,"0",sep=""),paste(cats,cats,sep=""),paste("0",cats,sep=""))
   
   # all possible observable capture histories("00" is unobservable)
   expected.Data$History = c("U0", "UU", cap.hist,"0U") 
   expected.Data$category = cats
   
   # add "00" capture history to end of the vector of capture history of expected.data
   hist = c(expected.Data$History,"00")
   
   # the function "log.prob.history" is in the file "neg.log.likelihood" which was
   # used to find thelog probabilities
   exp.cats = expected.Data$category
   expected.log.prob = log.prob.history(param,hist,exp.cats)#call the function 'log.prob.history'. 
   # These log probabilities and the counts are in the exact 
   # order related to each capture history. 
   
   # expected probabilities
   expected.prob = exp(expected.log.prob) - (expected.log.prob==0) 
   # have to use "-(expected.log.prob==0)" because of  "prob[index]==0" used 
   # inside the neg.log.likelihood
   ## sum(expected.prob) # this is just to check. Sum of the probabilities should equal to 1
   
   N =   param$N
   
   # expected counts for all the capture histories including history "00"
   expected.counts.all = N * expected.prob 
   # variance for all the  capture histories including history "00"
   expected.variance.all = N*expected.prob*(1-expected.prob ) 
   # variance for  the capture histories without history "00"
   expected.Data$variance =  expected.variance.all[1:length(expected.variance.all)-1]
   #print.expected.counts =cat(paste("Expected count for the capture History",hist[1:length(hist)],
   #                          sep = " ", ":",expected.counts.all[1:length(expected.counts.all)],'\n'))
   
   # expected couts for all the capture 
   expected.count = expected.counts.all[1:length(expected.counts.all)-1] 
   ## histories without history "00"
   ## data frame for expected histories and counts
   ## df.expected.counts = data.frame(expected.Data$History , expected.count) 
   ## names(df.expected.counts) = c("expected.History","expected.count")
   
   expected.Data$counts = expected.count  
   
   return(expected.Data)
  
}


###############################################################################
########### Test function #####################################################
###############################################################################
# 
# test.expected.counts= function(){
#    data <- NULL
#    data$History <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U")
#    data$count   <- c(41,9,16,4,23,7,12,25,63)
#    data$category <- c("M","F")
#    str(data)
#    print(data)
# 
#    # Set up some values of the parameters for testing
#    param <- NULL
#    param$p <- matrix(c( 0.1848103,0.1597991 ,0.1829310,0.2142439),nrow=2)
#    param$lambda <- c(0.3666383,(1-0.3666383))
#    param$theta <- c(0.5, 0.37)
#    param$N = 591
# 
#    print(param)
# 
#    test=NULL
#    test$param = param
#    test$category = data$category 
#    return(expected.counts(test))
# }
# 
# test.expected.counts()

###############################################################################
########### End Test function #################################################
###############################################################################





