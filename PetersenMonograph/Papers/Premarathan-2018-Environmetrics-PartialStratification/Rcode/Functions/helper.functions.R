# Numerical version of "logit" and "expit" functions

# p scale to logit
logit = function(prob){
  result =log(prob/(1-prob))
  return(result)
}

#antilogit(i.e. from logit to p scale)
expit = function(logit.prob){ 
  # if logit value is too small or too big then the regular value is -Inf or Inf
  # Therefore to address that issue, keep the logit value between -10 and 10
  #logit.prob =  max(-10, min(10, logit.prob)) # this dosent work for vector
  
  for(i in 1:length(logit.prob)){
    logit.prob[i] = max(-10, min(10, logit.prob[i]))
  }
  result = 1/(1+exp(-logit.prob))
  return(result)
}

##############################################################################
# To pack lambda parameters

# regular probabilities to cumulative logit probabilities
# Input : lambda vector of length(prob)
# Output: Vector of logit cumulatives of lambda. lenght is length(prob)-1
lambda.to.cumulative.logit.lambda = function(prob){
  result = regular.to.cumulative.logit.prob(prob)
  return(result)
}

# The following transformation preserves the "sum to one" constraint. 
# lambda.to.cumulative.logit.lambda ; 
# Input: A vector of categary proportions  of length = length(prob)
# Ouput: A vector of length = length(prob)-1

# Takes a vector of length(prob) Probabilities and returns a transformed vector with length(prob)-1 terms.
# Where the terms are caculated as follows;
#
# P_i  = logit( P_i/(P_i + ... P_k) ) ,    i = 1..k-1
#
# Assumptions: The vector prob sums to 1

regular.to.cumulative.logit.prob = function(prob){
  return(logit(cum.prob(prob)))
}

# The following function transform  the vector of size length(prob)  to form of a vactor of size length(prob)-1 vector
# for converting the lambda vector to cumulative form

cum.prob = function(prob){
  cum.prob = rev(cumsum(rev(prob)))
  cp  = prob/cum.prob 
  cp[is.na(cp)]= 1
  return(cp[-length(cp)])
}   

##############################################################################
# To unpack lambda parameters

# Inverts logit  probabilities back to regular probabiliries
cumulative.logit.lambda.to.lambda = function(logit.prob){
  result = inv.logit.lambda.to.lambda(logit.prob)
  #result[result < .00001] = 0  # force all small entries to be zero
  return(result)
}

# The following transformation preserves the "sum to one" constraint.
# inv.logit.lambda.to.lambda; 
# Input: A vector of logit.prob probabilities with length(logit.prob) elements;
# Ouput: A vector of probabilities "p" of length(logit.prob)+1; this vector sums to 1.
#
# The input vector elements are in the following form
# p_i = logit( P_i/(P_i + ... P_k) ) where, i = 1,...,(k-1) 

inv.logit.lambda.to.lambda = function(logit.prob){
  prob = expit(logit.prob)
  p = inv.cum.prob(prob)
  return(p)
}

inv.cum.prob = function(prob){
  temp = cumprod(1-prob)
  temp = c(1,temp[-length(temp)])*prob
  return(c(temp,1-sum(temp)))
}


# using  Least Squares method
#Input : a Matrix X(Design Matrix) and Vector y 
#Output: Least Squares Solution
LS = function(x,y){
  b = solve(t(x) %*% x) %*% t(x) %*% y
  return(b)
}


# using Generalized inverse method
#Input : a Matrix X(Design Matrix) and Vector y 
#Output: generalized linear model Solution
GI = function(x,y){
  b = ginv(t(x) %*% x) %*% t(x) %*% y
  return(b)
}


