###############################################################################
#####         calculate Expected counts for each capture history for     ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################

###############################################################################
###############################################################################
# "PSkSCR.expected.counts" generate the expected counts for the observed 
# histories and also the the possible expected values for the possible observable 
# but not observed histories. The all possible observable but not observed 
# histories are group together and call it "OTHER"
#
# Input : capture probabilities, category proportions,sub-sample proportions, 
#         data and indicators
# output : expected counts for the all observed capture histories, and
#          combined expected value for the all observable but not observed 
#          capture histories.
#           (eg."000"- animals not captured)
#        : variance for each expected count

PSkSCR.expected.counts = function(parm, data,indicator){
  
  history <- data$history
  counts  <- data$counts
  
  N <- parm$N # MLE for Population size
 
  # expected log probabilities
  expected.log.prob  <- PSkSCR.log.prob.history(parm,indicator)
  
  # expected probabilities
  expected.prob = exp(expected.log.prob) - (expected.log.prob==0) 
  

  # probability for capture history of unseen animals
  # (i.e. for two sampling occasions "00", for 3 sampling occasion "000" , ..)
  prob.hist.00 <- expected.prob[length(expected.prob)]
  
  # total probability of all other observable but not observed histories
  prob.hist.other <- 1 - ( sum(expected.prob))
  
  # probability for all capture possible observable histories
  prob.hist.all<- c( expected.prob[1:(length(expected.prob)-1)], prob.hist.other )
  
  get.sign<- c(sign(counts),1)
  
  expected.data = NULL
  # probability for all possible observable  histories
  # (i.e. for example  "000" is unobservable)
  expected.data$history <- c(history,"OTHER") 
  expected.data$counts <- N *prob.hist.all*get.sign
  expected.data$category <-  data$category
  
  # variance for  the capture histories without history "00"
  expected.data$variance <-  N*prob.hist.all*(1-prob.hist.all)
  expected.data$prob.hist.all <- prob.hist.all
  
  
  return(expected.data)
  
} # End of "PSkSCR.expected.counts"

###############################################################################
###############################################################################





