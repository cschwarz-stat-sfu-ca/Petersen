############################################################################### 
#################  Create a planning tool for optimization ####################
###############################################################################

# Functions in this file
#   "planning.tool.optimal.allocation" # with linear "equality" cost function
#   "print.optimal.allocation.equality.constraint" # print results from above function
#   "constraint.function"
#   "SE.N.cal" - this is used for calculating the standard 
#                error of N for optimal allocation 
#   "extract.parameters"
#   "starting.param.optimal.allocation"
#   "planning.tool.optimal.allocation.enequality.cost" # with linear "inequality" 
#                                                        cost function
#
#   "contour.plot.SE.N" - this is used for calculating the standard error of N 
#                         when creating contour plots
#
###############################################################################
###############################################################################

# with linear "equality" cost function

# FUNCTION: "planning.tool.optimal.allocation"
#   Inputs : vector of guesstimates(N, cat_prop,r1,r2), vector of categories,
#            cost vector ( c0, c1,c2, c1.star, c2.star), total cost, 
#            initial values for the parameters (i.e. initial counts for 
#            the capture Histories U0, UU, 0U, C0, CC, 0C), 
#            lower bounds, and the upper bounds for the for optimization routine
#
#   Output : optimized counts for the capture Histories U0, UU, 0U, C0, CC, 0C
#            where "C" represents all the categories,
#            Optimised values for n1. n1.star, n2, and n2.star and
#            standard error of N at the optimized values.


planning.tool.optimal.allocation <- function(N, cat_prop,r1,r2,
                                            category,
                                            cost.vec, total.cost, start.parm.vec,
                                            lower.b, upper.b,
                                            model.id, 
                                            captureDM,thetaDM, lambdaDM, 
                                            captureOFFSET, thetaOFFSET, lambdaOFFSET){
  
  res <- solnp(start.parm.vec, fun = SE.N.cal, 
               eqfun = constraint.function, eqB =c(total.cost),
               LB = lower.b, UB = upper.b ,
               N=N, cat_prop=cat_prop,
               r1=r1,r2=r2, category=category,cost.vec=cost.vec,
               model.id=model.id, 
               captureDM=captureDM,thetaDM=thetaDM, lambdaDM=lambdaDM, 
               captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,lambdaOFFSET=lambdaOFFSET)
  
  finalparam <- res$pars # optimized values for n_U0, n_UU, n_0U, n_C0, n_CC, n_0C 
  names(finalparam) <- c("n_U0", "n_UU", "n_0U", "n_C0", "n_CC", "n_0C")
  
  finalSE <- res$values # standard error of N at optimized values above
  
  # optimized values for n1, n1.star, n2, and n2.star
       n1 <- finalparam[1] + finalparam[2] + finalparam[4] + finalparam[5]
  n1.star <- finalparam[4] + finalparam[5]
       n2 <- finalparam[2] + finalparam[5] + finalparam[3] + finalparam[6]
  n2.star <- finalparam[6]
  
  
  optimal.allocation <- c(n1, n1.star, n2, n2.star)
  names(optimal.allocation) <- c("n1", "n1.star", "n2", "n2.star")
  
  final.cost <- sum(c(1, n1, n1.star, n2, n2.star)*cost.vec)
  
  ######################
  out <- NULL
  
  out$model.id <- model.id
  out$total.cost <- total.cost # initial avilable total cost
  out$start.parm.vec <- start.parm.vec  ### starting parameter vector
  out$cost.vec <- cost.vec # cost for sampling and sub-sampling 
                          #  at each time and fixed cost
  out$finalparam <- finalparam
  out$finalSE <- finalSE
  out$optimal.allocation <- optimal.allocation
  out$final.cost <- final.cost
  
  return(out)

} #end of the function "planning.tool.optimal.allocation"

###############################################################################
##### Print the Optimal allocation results under linear equality Constraint ###

# Function  : print.optimal.allocation.equality.constraint
#  Input : output from the function "planning.tool.optimal.allocation"
#  Output : print the relults
  
print.optimal.allocation.equality.constraint <- function(result){
  cat("\n")
  cat( "############################################################")
  cat("\n")
  cat("Optimal Allocation Results Under Linear Equality Constraint \n")
  cat(" (Equality total cost Constraint) \n")
  cat("\n")
  cat("Model ID: ",result$model.id)
  cat("\n")
  cat("\n")
  cat("cost for sampling and sub-sampling at each time and fixed cost: \n")
  print(result$cost.vec)
  cat("\n")
  cat("Maximum total cost available: ",result$total.cost)
  cat("\n")
  cat("\n")
  cat("initial values of n_U0, n_UU, n_0U, n_C0, n_CC, n_0C for optimization routing: \n")
  print(result$start.parm.vec)
  cat("\n")
  cat( "--------------------")
  cat("\n")
  cat("\n")
  cat("Optimized values for n_U0, n_UU, n_0U, n_C0, n_CC, n_0C: \n")
  print(round(result$finalparam,2))
  cat("\n")
  cat("\n")
  cat("Vector of Standard error of N during optimization: 
      Last one is the value at the optimal (SE_N): \n")
  print(round(result$finalSE),2)
  cat("\n")
  cat("\n")
  cat("Optimized values for n1, n1.star, n2, and n2.star: \n")
  print(round(result$optimal.allocation,2))
  cat("\n")
  cat("\n")
  cat("Total cost at the optimal allocation:" , round(result$final.cost,2))
  cat("\n")
  cat("\n")
  cat( "############################################################")
  cat("\n")
} #end of the function "print.optimal.allocation.equality.constraint"



###############################################################################
###############################################################################
# Function : constraint.function
#       constraint function returning the vector of evaluated 
#       equality constraints or  inequality constraints.
#   Input : param.count, N, cat_prop,r1,r2, category,cost.vec
#   Output : cost

constraint.function <- function(param.count,
                                N, cat_prop,r1,r2, 
                                category,cost.vec,
                                model.id, 
                                captureDM,thetaDM, lambdaDM, 
                                captureOFFSET, thetaOFFSET, lambdaOFFSET){
  
  n.count <- param.count
  
       n1 <- n.count[1] + n.count[2] + n.count[4] + n.count[5]
  n1.star <- n.count[4] + n.count[5]
       n2 <- n.count[2] + n.count[5] + n.count[3] + n.count[6]
  n2.star <- n.count[6]
  
  cost <- sum(c(1, n1, n1.star, n2, n2.star)*cost.vec)
  return(cost)
  
} # end of the function "constraint.function"

###############################################################################
###############################################################################

# FUNCTION SE.N 
#   Inputs : counts for the capture Histories U0, UU, 0U, C0, CC, 0C and 
#            guesstimates (N, cat_prop,r1,r2), vector of categories,
#            vector of costs values
#   Outputs : standard error of N 

SE.N.cal = function(param.count, 
                    N, cat_prop,r1,r2, 
                    category,cost.vec,
                    model.id, 
                    captureDM,thetaDM, lambdaDM, 
                    captureOFFSET, thetaOFFSET, lambdaOFFSET){
  
  res = NULL
  res$param = extract.parameters(param.count, N, cat_prop,r1,r2)
  res$category = category
  
  # get the expected histories and counts acording to parameters
  expected.Data = expected.counts(res) # expected histories and expected counts
  
  # fit MLE model to expected data
  expected.model.Result =  fit.model(model.id,expected.Data,
                                     captureDM,thetaDM,lambdaDM,
                                     captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  # get the standard error of population size N
  SE_N = sqrt(expected.model.Result$vcv$N)

  #print(expected.model.Result$est$N)
  #print(expected.model.Result$est$lambda)
  #print(expected.model.Result$est$p)
  
  #### this is just to see at the each iteration ####
  #print("parameters")
  #print(param.count)
  #print("SE")
  #print(SE_N)
  ###################################################
  return(SE_N)
  
} # end of the function "SE.N.cal"



###############################################################################
###############################################################################

# FUNCTION exatract.parameters 
#   Inputs : counts for the capture Histories U0, UU, 0U, C0, CC, 0C and 
#            guesstimates (N, cat_prop,r1,r2)
#   Outputs : vector of parameters ( capture probabilities (p1m,p1f,p2m,p2f), 
#             category proportions(lambda), sub-sample proportions (theta), N)

extract.parameters = function (param.count, N, cat_prop,r1,r2){
  
  n.count <- param.count
  
  n1 <-  n.count[1]+n.count[2] + n.count[4] + n.count[5]
  n1.star <- n.count[4] +n.count[5]
  n2 <- n.count[2] + n.count[5] + n.count[3] + n.count[6]
  n2.star <-  n.count[6]
  
  
  R1 = r1*cat_prop[1]/cat_prop[2]  # ratio of male/female in the first occasion
  R2 = r2*cat_prop[1]/cat_prop[2]  # ratio of male/female in the second occasion
  
  p.1m = ( n1 * (R1/(R1+1))) / (N * cat_prop[1])
  p.1f = ( n1 * (1/(R1+1))) / (N * cat_prop[2])
  p.2m = ( n2 * (R2/(R2+1))) / (N * cat_prop[1])
  p.2f = ( n2 * (1/(R2+1))) / (N * cat_prop[2])       
  
  
  theta.1 = n1.star/n1
  theta.2 = n2.star/(n2.star + n.count[3] )
  
  parameters = NULL
  
  parameters$p      = matrix(c(p.1m,p.1f,p.2m,p.2f), ncol=2,byrow=FALSE)
  parameters$lambda = cat_prop
  parameters$theta  = c(theta.1,theta.2)
  parameters$N      = N
  parameters$full   = c(p.1m,p.1f,p.2m,p.2f,parameters$lambda, parameters$theta,N)
  
  return(parameters)
  
} # end of the function "exatract.parameters" 

###############################################################################
## Function : starting.param.optimal.allocation
##     Input : Data set (Mille Lacs Walleye dataset)
##     output: starting parameter vector is in the following oder
##             ( n_U0, n_UU, n_0U, n_C0, n_CC, n_0C )

starting.param.optimal.allocation <- function(Data){
  
  
  n_U0 <- sum(Data$counts[(substr(Data$History,1,1)=="U") & (substr(Data$History,2,2)=="0")])
  
  n_UU <- sum(Data$counts[(substr(Data$History,1,1)=="U") & (substr(Data$History,2,2)=="U")])
  
  #total number of individuals caught in  sample 1
  total.time.1<- sum(Data$counts[ ! (substr(Data$History,1,1)=="0" )]) 
  
  # m2 is the number of marked individuals caught in sample 2 ( e.g. UU, MM, FF,...)
  m2 <- total.time.1 - sum(Data$counts[ (substr(Data$History,2,2)=="0" )]) 
  
  n_CC <- m2- n_UU 
  n_C0 <- total.time.1 - (n_U0+n_UU) - n_CC
  
  n_0U <- sum(Data$counts[(substr(Data$History,1,1)=="0") & (substr(Data$History,2,2)=="U")])
  n_0C <- sum(Data$counts[ (substr(Data$History,1,1)=="0")]) - sum(Data$counts[Data$History=="0U"])
  
  
  start.parm.vec <- c(n_U0, n_UU, n_0U, n_C0, n_CC, n_0C )
  
  return(start.parm.vec)
}# end of the function "starting.param.optimal.allocation" 

###############################################################################
###############################################################################

# with linear "inequality" cost function

# Function : planning.tool.optimal.allocation.enequality.cost

planning.tool.optimal.allocation.inequality.cost <- function(N, cat_prop,r1,r2,
                                                            category,
                                                            cost.vec, 
                                                            inequalityLB,inequalityUB, 
                                                            start.parm.vec,
                                                            lower.b, upper.b,
                                                            model.id, 
                                                            captureDM,thetaDM, lambdaDM, 
                                                            captureOFFSET, thetaOFFSET, lambdaOFFSET){
  
  res <- solnp(start.parm.vec, fun = SE.N.cal, 
               ineqfun = constraint.function,
               ineqLB =c(inequalityLB) , ineqUB =c(inequalityUB),
               LB = lower.b, UB = upper.b ,
               N=N, cat_prop=cat_prop,
               r1=r1,r2=r2, category=category,cost.vec=cost.vec,
               model.id=model.id, 
               captureDM=captureDM,thetaDM=thetaDM, lambdaDM=lambdaDM, 
               captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,lambdaOFFSET=lambdaOFFSET)
  
  finalparam <- res$pars # optimized values for n_U0, n_UU, n_0U, n_C0, n_CC, n_0C 
  names(finalparam) <- c("n_U0", "n_UU", "n_0U", "n_C0", "n_CC", "n_0C")
  
  finalSE <- res$values # standard error of N at optimized values above
  
  # optimized values for n1, n1.star, n2, and n2.star
       n1 <- finalparam[1] + finalparam[2] + finalparam[4] + finalparam[5]
  n1.star <- finalparam[4] + finalparam[5]
       n2 <- finalparam[2] + finalparam[5] + finalparam[3] + finalparam[6]
  n2.star <- finalparam[6]
  
  
  optimal.allocation <- c(n1, n1.star, n2, n2.star)
  names(optimal.allocation) <- c("n1", "n1.star", "n2", "n2.star")
  
  final.cost <- sum(c(1, n1, n1.star, n2, n2.star)*cost.vec)
  
  ######################
  out <- NULL
  
  out$model.id <- model.id
  out$inequalityLB <- inequalityLB # lower bound for available total cost
  out$inequalityUB <- inequalityUB # upper bound for available total cost
  out$start.parm.vec <- start.parm.vec  ### starting parameter vector
  out$cost.vec <- cost.vec # cost for sampling and sub-sampling 
                           #  at each time and fixed cost
  out$finalparam <- finalparam
  out$finalSE <- finalSE
  out$optimal.allocation <- optimal.allocation
  out$final.cost <- final.cost
  
  return(out)
  
} # end of the function "planning.tool.optimal.allocation.inequality.cost"  
  

###############################################################################
##### Print the Optimal allocation results under linear inequality Constraint ###

# Function  : print.optimal.alloc.inequality.constraint
#  Input : output from the function "planning.tool.optimal.allocation.inequality.cost"
#  Output : print the relults

print.optimal.alloc.inequality.constraint <- function(result){
  cat("\n")
  cat( "############################################################")
  cat("\n")
  cat("Optimal Allocation Results Under Linear Inequality Constraint \n")
  cat(" (inequality total cost Constraint) \n")
  cat("\n")
  cat("Model ID: ",result$model.id)
  cat("\n")
  cat("\n")
  cat("cost for sampling and sub-sampling at each time and fixed cost: \n")
  print(result$cost.vec)
  cat("\n")
  cat("Lower bound for available total cost: ",result$inequalityLB)
  cat("\n")
  cat("Upper bound for available total cost: ",result$inequalityUB)
  cat("\n")
  cat("\n")
  cat("initial values of n_U0, n_UU, n_0U, n_C0, n_CC, n_0C for optimization routing: \n")
  print(result$start.parm.vec)
  cat("\n")
  cat( "--------------------")
  cat("\n")
  cat("\n")
  cat("Optimized values for n_U0, n_UU, n_0U, n_C0, n_CC, n_0C: \n")
  print(round(result$finalparam,2))
  cat("\n")
  cat("\n")
  cat("Vector of Standard error of N during optimization: 
      Last one is the value at the optimal (SE_N): \n")
  print(round(result$finalSE),2)
  cat("\n")
  cat("\n")
  cat("Optimized values for n1, n1.star, n2, and n2.star: \n")
  print(round(result$optimal.allocation,2))
  cat("\n")
  cat("\n")
  cat("Total cost at the optimal allocation:" , round(result$final.cost,2))
  cat("\n")
  cat("\n")
  cat( "############################################################")
  cat("\n")
} #end of the function "print.optimal.alloc.inequality.constraint"

  
###############################################################################
###############################################################################
###############################################################################


######################## CONTOUR PLOT #########################################
###                                                                         ###
## Contour plot for standard error of N for given n1, n1.sta, n2, and n2.star##            
###  n1, n2 are fixed at the optimised values.                              ###
###  or n1.star and n2.star are fixed at the optimised values.              ###
###                                                                         ###
###############################################################################



###############################################################################
contour.plot.SE.N <- function (n.counts, N, cat_prop,r1,r2,
                               model.id, 
                               captureDM,thetaDM, lambdaDM, 
                               captureOFFSET, thetaOFFSET, lambdaOFFSET){
  
       n1 <- n.counts[1]
  n1.star <- n.counts[2]
       n2 <- n.counts[3]
  n2.star <- n.counts[4]
  
  R1 <- r1*cat_prop[1]/cat_prop[2]  # ratio of male/female in the first occasion
  R2 <- r2*cat_prop[1]/cat_prop[2]  # ratio of male/female in the second occasion
  
  p.1m <- ( n1 * (R1/(R1+1))) / (N * cat_prop[1])
  p.1f <- ( n1 * (1/(R1+1))) / (N * cat_prop[2])
  p.2m <- ( n2 * (R2/(R2+1))) / (N * cat_prop[1])
  p.2f <- ( n2 * (1/(R2+1))) / (N * cat_prop[2])       
  
  
  theta.1 <- n1.star/n1
  
  n_UU <- N* ( cat_prop[1]*p.1m*(1-theta.1)*p.2m + cat_prop[2]*p.1f*(1-theta.1)*p.2f )
  n_CC <- N* ( cat_prop[1]*p.1m*theta.1*p.2m + cat_prop[2]*p.1f*theta.1*p.2f )
  
  theta.2 <- n2.star/(n2- (n_UU +n_CC ) )
  
  parameters <- NULL
  
  parameters$p      <- matrix(c(p.1m,p.1f,p.2m,p.2f), ncol=2,byrow=FALSE)
  parameters$lambda <- cat_prop
  parameters$theta  <- c(theta.1,theta.2)
  parameters$N      <- N
  parameters$full   <- c(p.1m,p.1f,p.2m,p.2f,parameters$lambda, parameters$theta,N)
  
  ### following is due to the way I used in the function "expected.counts" ####
  result <- NULL
  result$param <- parameters
  result$category <- category
  #############################################################################
  
  expected.Data <- expected.counts(result) # expected histories and expected counts

  # fit MLE model 
  expected.model.Result =  fit.model(model.id,expected.Data,
                                     captureDM,thetaDM,lambdaDM,
                                     captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  # standard error of N
  SE_N = sqrt(expected.model.Result$vcv$N)
  #print(SE_N)
  return(SE_N)
  
} # end of the function "contour.plot.SE.N" 

###############################################################################
###############################################################################

