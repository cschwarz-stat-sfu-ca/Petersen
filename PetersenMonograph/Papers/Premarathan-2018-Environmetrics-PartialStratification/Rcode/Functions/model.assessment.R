###############################################################################
#########  Model assessment for the MLE model  ################################
###############################################################################
# Function residual.plot
# input : list of values from the fitted model
#output : residual plot and the capture histories, observed and expected counts

residual.plot = function(model) {

  model.id= model$model.id 
  Data =model$rawdata
  #Find expected counts for model1:
  results= NULL
  results$param = model$est
  results$category = model$rawdata$category
  exp.data = expected.counts(results)  # exp.data is a list with  expected capture 
                                       # histories, corresponding counts and categories

  # Observed.counts
  #  This produces the observed counts as same order of the expected capture history.
  #  If some histories are not present in the 
  #  observed data, then their counts will be zero

  #sum of the counts for each of history in Data$History
  obs.counts = outer(exp.data$History, Data$History, '==') %*% Data$count 
  residual = obs.counts - exp.data$counts # calculated the residuals

  # calculated the standardized residuals
  standardized.residual = residual/sqrt(exp.data$variance) 

  # Output.df gives the Histories, expected and observed counts and residuals
  Output.df =data.frame(History = exp.data$History,Observed.Counts =obs.counts,
                        Expected.counts =exp.data$counts, 
                        Residual = residual, 
                        Standardized.Residuals = standardized.residual)

  #plot the standardized residuals
  res.plot = ggplot( Output.df,aes(x = Expected.counts, y =Standardized.Residuals,
                                   label=History)) +
                  geom_point(size = 3,colour="#990000") +
                  labs(x = "Counts",y = "Standardized Residuals",
                       title = paste("Residual Plot: model", model.id,sep=" " )) +
                  theme(axis.text=element_text(size=10,face="bold"),
                        axis.title=element_text(size=14),
                        plot.title = element_text(size = 12,face="bold")) +
                  geom_hline(yintercept=c(-1.96,0,1.96),
                             colour=c("#660000","black","#660000"))+ 
                  geom_text(aes(label=History),hjust=0, vjust=-.3)

  model.out = NULL
  model.out$model.id = model.id
  model.out$obs.exp.counts = Output.df[1:3]
  model.out$res.plot = res.plot
             
  return(model.out)
} # end of residual.plot


###############################################################################
###############################################################################



