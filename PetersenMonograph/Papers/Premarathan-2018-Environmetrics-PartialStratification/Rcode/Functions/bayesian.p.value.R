
###############################################################################
######### Posterior Predictive plots and Bayesian p-value using ###############
########## two methods (likelihood and Freeman-Tukey statistic) ###############
###############################################################################

# Function: baysian.p.value
# inputs : model information after fitting the bayesian model
# output : Posterior Predictive plots and Bayesian p-value

bayesian.p.value <- function(bayes.model.info){
  
#   # thin.beta.matrix store all the (thinned) posterior samples in all the chains 
#   thin.beta.matrix <- matrix( ncol=bayes.model.info$np, nrow =0)
#   
#   for( i in 1:bayes.model.info$n.chains ){
#     thin.beta.matrix <- rbind(thin.beta.matrix, 
#                               bayes.model.info$MCMC$thin.beta.array[ , , i])
#   }
#   

#  thinned.beta.matrix store all the (thinned) posterior samples in all the chains 
  thinned.beta.matrix <- bayes.model.info$thinned.beta.matrix
  
  
  # model id and design matrices and offset vectors
  model.id <- bayes.model.info$model.id
  
  captureDM <- bayes.model.info$captureDM
  thetaDM <- bayes.model.info$thetaDM
  lambdaDM <- bayes.model.info$lambdaDM
  
  captureOFFSET <- bayes.model.info$captureOFFSET
  thetaOFFSET <- bayes.model.info$thetaOFFSET
  lambdaOFFSET <- bayes.model.info$lambdaOFFSET
  
  
  # get unique capture histories from all the histories in the data set
  uniq.history = unique(bayes.model.info$Data$History)  
  
  #sum of the counts for each of history in Data$History
  obs.counts = outer(uniq.history, bayes.model.info$Data$History,
                     '==') %*% bayes.model.info$Data$count 
  
  unique.Data = NULL
  unique.Data$History = uniq.history #unique histories
  unique.Data$counts = obs.counts # total counts for related unique histories
  unique.Data$category = bayes.model.info$Data$category
  
  
  
  # negative log likelihood for each posterior (thinned) sample
  NLL.obs <- adply(thinned.beta.matrix,1,neg.log.likelihood,
               unique.Data,captureDM,thetaDM,lambdaDM,
               captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  
  
  res <- adply(thinned.beta.matrix,1,neg.likelihood.and.tukey,
                      Data,captureDM,thetaDM,lambdaDM,
                      captureOFFSET,thetaOFFSET,lambdaOFFSET)
  
  # negative log likelihood for simulated samples at each MCMC iteration
  NLL.sim <- res[,2]
  
    
  ###### Bayesian p-value using deviance #################
  deviance.p.value <- mean(2*NLL.obs[,2] < 2*NLL.sim)
  
  deviance.df <- data.frame(cbind(2*NLL.obs[,2],2*NLL.sim  ))
  names(deviance.df) <- c("NLL.obs", "NLL.sim")
  
  xymin <- 2*min(NLL.obs[,2],NLL.sim)
  xymax <- 2*max(NLL.obs[,2],NLL.sim)
  
  deviance.plot <- ggplot(deviance.df, aes(x=NLL.obs, y=NLL.sim)) +
                      geom_point()+
                      geom_abline(slope=1, intercept=0)+
                      xlab("Observed Deviance") + ylab("Simulated Deviance") + 
                      xlim(xymin, xymax)+
                      ylim(xymin, xymax)+
                      annotate("text",label = paste("Bayesian p-value = ",
                                              round(deviance.p.value,2), sep=""),
                               x =(xymin+0.2*(xymax-xymin)), hjust = 0, y = Inf, 
                               vjust = 2, color = "darkred")
  
  
  ##### Bayesian p-value using Freeman-Tukey Statistic #######
  obs.FT.statistic <- res[,3]
  sim.FT.statistic <- res[,4]
  FT.p.value <- mean( obs.FT.statistic < sim.FT.statistic )
  
  FT.statistic.df <- data.frame(cbind(obs.FT.statistic,sim.FT.statistic ))
  names(deviance.df) <- c("obs.FT.statistic", "sim.FT.statistic ")
  
  
  FT.xymin <- min(obs.FT.statistic,sim.FT.statistic)
  FT.xymax <- max(obs.FT.statistic,sim.FT.statistic)
  
  FT.plot <- ggplot(FT.statistic.df, aes(x=obs.FT.statistic, 
                                               y=sim.FT.statistic)) +
                      geom_point()+
                      geom_abline(slope=1, intercept=0)+
                      xlab("Observed FT Statistic") + ylab("Simulated FT Statistic") + 
                      xlim(FT.xymin, FT.xymax)+
                      ylim(FT.xymin, FT.xymax)+
                      annotate("text", label = paste("Bayesian p-value = ",
                                                    round(FT.p.value,2), sep=""),
                             x =(FT.xymin+0.2*(FT.xymax-FT.xymin)), hjust = 0,y = Inf, 
                             vjust = 2, color = "darkred")
  
  predictive.plots <- arrangeGrob(deviance.plot, FT.plot, nrow = 1, 
                                  main = textGrob(
                      "Bayesian p-value scatter plots using the Discrepancy functions 
                      (a) deviance  (b) Freeman-Tukey (FT) statistic",
                      just = "top", vjust = 0.75, gp = gpar(fontface = "bold")))
  
  output <- NULL
  output$deviance.p.value <- deviance.p.value
  output$FT.p.value <- FT.p.value
  output$predictive.plots <- predictive.plots
  
  return(output) 
  
} # end of bayesian.p.value
  
  

###############################################################################
###############################################################################
# Function: neg.likelihood.and.tukey calculate 
# Input : vector of parameters(x) in logit/log scale, Data,design matrices
#          and offset vectors
# Output : Negative likelihood and Freeman-Tukey statistic 
#          for a one generated data set

neg.likelihood.and.tukey <- function(x,unique.Data,
                                    captureDM,thetaDM,lambdaDM,
                                    captureOFFSET,thetaOFFSET,lambdaOFFSET ){
    
    param <- unpack.parm(x,Data,captureDM,thetaDM,lambdaDM,
                         captureOFFSET,thetaOFFSET,lambdaOFFSET)  
  
    cats = unique.Data$category
  
    cap.histories = c(paste(cats,"0",sep=""),paste(cats,cats,sep=""),
                      paste("0",cats,sep=""))
  
    gen.Data = NULL
    # all possible observable capture histories("00" is unobservable)
    gen.Data$History = c("U0", "UU", cap.histories,"0U") 
    gen.Data$category = cats
  
    # add "00" capture history to end of the vector of capture history of expected.data
    hist = c(gen.Data$History,"00")
  
    # the function "log.prob.history" is in the file "neg.log.likelihood" which
    # was used to find the log probabilities
    # call the function 'log.prob.history'. These log probabilities and the
    # counts are in the exact  order related to each capture history.
    expected.log.prob = log.prob.history(param,hist,cats)   
     
    # expected probabilities
    expected.prob = exp(expected.log.prob) - (expected.log.prob==0) 
    # expected probabilities for observable capture histories(without history"00")
    obs.hist.exp.prob  =  expected.prob[1:(length(expected.prob)-1)] 
  
    # probability of the capture history "00"
    p_00 <- sum(param$lambda*(1-param$p[,1])*(1-param$p[,2])) 
  
    # conditional probabilities for  observable capture histories
    cond.obs.hist.exp.prob <- obs.hist.exp.prob/(1-p_00)  
  
      
    sample.size <- sum(unique.Data$counts)
    gen.Data$counts <- rmultinom(1, size = sample.size, 
                                 prob = cond.obs.hist.exp.prob) 
  
    
    NLL.gen <- neg.log.likelihood(x, gen.Data,captureDM,thetaDM,lambdaDM,
                                  captureOFFSET,thetaOFFSET,lambdaOFFSET)
    
    ######### Freeman-Tukey Statistic  ########
    
    #Find expected counts for generated data
    MLE.results <- NULL
    MLE.results$param <-   param
    MLE.results$category <-   unique.Data$category
    # exp.data is a list with  expected capture histories, corresponding counts
    # and categories and variance for each expected category
    exp.data <- expected.counts(MLE.results)  
        
    ####### Freeman-Tukey Statistic  for simulated data ##################
    
    tukey.statistic.simulated.data <- sum((sqrt(gen.Data$count) - 
                                             sqrt(exp.data$count))^2)
      
    ####### Freeman-Tukey Statistic  for observed data ###################
    
    # obs.counts
    #  This produces the observed counts as same order of the expected capture
    #  history. If some histories are not present in the observed data, then 
    #  their counts will be zero
    obs.counts = outer(exp.data$History, Data$History, '==') %*% Data$count
    
    tukey.statistic.observed.data <- sum((sqrt(obs.counts) - 
                                           sqrt(exp.data$count))^2)
    
    ###########################################################################
            
    return(c(NLL.gen = NLL.gen, obs.FT.ststistic = tukey.statistic.observed.data,
                sim.FT.statistic = tukey.statistic.simulated.data))
    
}  # end of neg.likelihord.and.tukey
  
###############################################################################
###############################################################################
