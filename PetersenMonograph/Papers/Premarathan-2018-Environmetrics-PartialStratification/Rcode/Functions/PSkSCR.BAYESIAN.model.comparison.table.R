###############################################################################
#####                 Bayesian Model Comparision Table                   ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################

#Function : PSkSCR.BAYESIAN.model.comparison.table
# Input: method and output of List of Models
#        Default method is "Spiegelhalter"
# Output: table(sorted by DIC value) with modelID, number of parameters, N, SE of N, DIC, 
#          DIC weights and Delta DIC(difference of DIC from the minimum)

PSkSCR.BAYESIAN.model.comparison.table = function(method="Spiegelhalter",models=list()){
  
  # make a table
  table = ldply(models, function(x){
    model.id <- x$model.id # Extract the name 
    np <- x$np
    #N <- round(x$est$N)
    #N.se <- round(x$se$N)
    Dbar <- round(x$DIC$Dbar,3)
    pD <- round(x$DIC$pD,3)
    pv <- round(x$DIC$pv,3)
    if(method=="Spiegelhalter"){DIC <- round(x$DIC$DIC$pD,2)}
    if(method=="Gelman"){DIC <- round(x$DIC$DIC$pv,2)}
    
   # result= data.frame(model.id=model.id,np=np, N=N, N.se=N.se,DIC=DIC )
    result= data.frame(model.id=model.id,np=np,Dbar=Dbar,pD=pD,pv=pv,DIC=DIC )
    result
  })
  
  
  #  sort by DIC
  table = table[order(table$DIC),]
  
  table$DeltaDIC = table$DIC - min(table$DIC)  #difference of AICc from the minimum
  
  table$DeltaDIC = round(table$DeltaDIC,2)
  
  table$DIC.weights = round(exp(-table$DeltaDIC/2)/sum(exp(-table$DeltaDIC/2)),2)
  nrow(table)
  
  
  return(table)
} # end of PSkSCR.BAYESIAN.model.comparison.table

###############################################################################
###############################################################################