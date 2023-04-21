# model.comparision.table
# Input: List (List of Models)
# Output: table(sorted by AICc value) with modelID, number of parameters, 
#         log likelihood, N, SE of N, AICc, AICc weights and Dealta AICc
#         (diference of AICc from the minimum)
#         and maximum of 4 residual plots in 2x2 grid for the models with
#         minimum AICc values(top 4 of the sorted "table")

###############################################################################
###############################################################################

model.comparison.table = function(models=list()){
  
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
  
  #diference of AICc from the minimum
  table$DeltaAICc = table$AICc - min(table$AICc)  
  table$DeltaAICc = round(table$DeltaAICc,2)  # round off the values
  
  # AICc weights 
  table$AICc.weights = round(exp(-table$DeltaAICc/2)/sum(exp(-table$DeltaAICc/2)),2)
    
  model_comparison$table <- table
  
  ######################################################
  ### create maximum of 4 residual plots in 2x2 grid for the models with 
  ### minimum AICc values(top 4 of the sorted "table")
  
  n.models = length(models)  # total  number of models consider
  
  # consider only the top 4 models(according to the minimum AICc value) if 
  # number of models greater than 4. otherwise consider the all the given
  # models(less than 4 models)
  k = min(4, n.models) 
                        
  # extract names of all the models 
  model.names = ldply(models, function(x){return(x$model.id)})  
  
  # get the names of the top k models in the sorted table
  table.top.models = as.vector(table$model.id[1:k]) 
  
  # compare the names of list of top k models in the table with all the model 
  # names and return the index for the model
  model.index = outer(table.top.models, model.names[,1], '==') %*%  seq(1:n.models) 
                                                           
  # extract the residual plots for the top k models in the table
  extract.res.plots = llply(model.index, function(x,models){models[[x]]$res.plot},
                            models=models)
  
  
  # create the residual plots for the top k models in a grid
  residual.plots <- do.call(arrangeGrob, c(extract.res.plots,nrow = 2,
                                   main=paste("Residual Plots with minimum ",
                                  "AICc values (upto maximum 4 plots)", sep="")))
  
  model_comparison$plots <- residual.plots
  
  return(model_comparison)
} # end of  model.comparison.table

###############################################################################
###############################################################################
## AICc Weights
# Akaike weights are can be used in model averaging. They represent the 
# relative likelihood of a model. To calculate them, for each model first 
# calculate the relative likelihood of the model, which is
# just exp( -0.5 * DELTA(AIC) score for that model). The Akaike weight for 
# a model is this value divided by the sum of these values across all models.


