############################################################################### 
####### model.comparison.table  for the data with individual covariates  ######
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
############################################################################### 

# model.comparison.table
# Input: List (List of Models)
# Output: table(sorted by AICc value) with modelID, number of parameters, 
#         log likelihood, N, SE of N, AICc, AICc weights and Delta AICc
#         (difference of AICc from the minimum)

###############################################################################
###############################################################################
 
ic.model.comparison.table = function(models=list()){
  
  model_comparison <- NULL
  
  # make a table
  table = ldply(models, function(x){
              model.id = x$model.id # Extract the name 
                    np = x$np
                 log.L = round(-x$NLL,2)
                     N = round(x$est$N)
                  N.se = round(x$se$N)
                  AICc = round(x$AICc,2)
                                
    result= data.frame(model.id=model.id,np=np, log.L=log.L, N=N,
                       N.se=N.se,AICc=AICc )
    result
  })
  
 
  #  sort by AICc 
  table = table[order(table$AICc),]
  
  # difference  of AICc from the minimum
  table$DeltaAICc = table$AICc - min(table$AICc)  
  table$DeltaAICc = round(table$DeltaAICc,2)  # round off the values
  
  # AICc weights 
  table$AICc.weights = round(exp(-table$DeltaAICc/2)/sum(exp(-table$DeltaAICc/2)),2)
    
  return(table)
} # end of ic.model.comparison.table

###############################################################################
###############################################################################
## AICc Weights
# Akaike weights are can be used in model averaging. They represent the 
# relative likelihood of a model. To calculate them, for each model first 
# calculate the relative likelihood of the model, which is
# just exp( -0.5 * DELTA(AIC) score for that model). The Akaike weight for 
# a model is this value divided by the sum of these values across all models.


