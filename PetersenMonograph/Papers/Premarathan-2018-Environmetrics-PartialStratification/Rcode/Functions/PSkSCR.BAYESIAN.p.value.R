###############################################################################
######### Posterior Predictive plots and Bayesian p-value using ###############
########## two methods (likelihood and Freeman-Tukey statistic) ###############
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ######
###############################################################################
###############################################################################

# Function: PSkSCR.BAYESIAN.p.value
# inputs : model information after fitting the bayesian model
# output : Posterior Predictive plots and Bayesian p-value

# bayes.model.info=BAYESIAN_PSkSCR_model_1 ##### delete later...  

PSkSCR.BAYESIAN.p.value <- function(bayes.model.info){
  
#   # thin.beta.matrix store all the (thinned) posterior samples in all the chains 
#   thin.beta.matrix <- matrix( ncol=bayes.model.info$np, nrow =0)
#   
#   for( i in 1:bayes.model.info$n.chains ){
#     thin.beta.matrix <- rbind(thin.beta.matrix, 
#                               bayes.model.info$MCMC$thin.beta.array[ , , i])
#   }
#   

# thinned.beta.matrix store all the (thinned) posterior samples in all the chains 
  thinned.beta.matrix <- bayes.model.info$thinned.beta.matrix
  
  
  # model id and design matrices and offset vectors
  model.id <- bayes.model.info$model.id
  
  captureDM <- bayes.model.info$captureDM
  thetaDM <- bayes.model.info$thetaDM
  lambdaDM <- bayes.model.info$lambdaDM
  p_lossDM <- bayes.model.info$p_lossDM
  
  captureOFFSET <- bayes.model.info$captureOFFSET
  thetaOFFSET <- bayes.model.info$thetaOFFSET
  lambdaOFFSET <- bayes.model.info$lambdaOFFSET
  p_lossOFFSET <- bayes.model.info$p_lossOFFSET
  
  indicator <- bayes.model.info$indicator
  
  # # get unique capture histories from all the histories in the data set
  # uniq.history = unique(bayes.model.info$data$history)  
  # 
  # #sum of the counts for each of history in Data$History
  # obs.counts = outer(uniq.history, bayes.model.info$Data$History,
  #                    '==') %*% bayes.model.info$Data$count 
  # 
  # unique.Data = NULL
  # unique.Data$History = uniq.history #unique histories
  # unique.Data$counts = obs.counts # total counts for related unique histories
  # unique.Data$category = bayes.model.info$Data$category
  
  obs.data <- bayes.model.info$data
  
  # negative log likelihood for each posterior (thinned) sample
  NLL.obs <- adply(thinned.beta.matrix,1,PSkSCR.neg.log.likelihood,
                   obs.data,indicator,
                   captureDM,thetaDM,lambdaDM,p_lossDM,
                   captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
 
  # results from simulated data
  res <- adply(thinned.beta.matrix,1,PSkSCR.neg.likelihood.and.tukey,
               obs.data,indicator,
               captureDM,thetaDM,lambdaDM,p_lossDM,
               captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)
  
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
                      paste("Bayesian p-value scatter plots using the Discrepancy functions \n",
                            "(a) deviance  (b) Freeman-Tukey (FT) statistic \n", 
                            "Model: ",model.id, sep=" "),
                      just = "top", vjust = 0.75, gp = gpar(fontface = "bold")))
  
  output <- NULL
  output$deviance.p.value <- deviance.p.value #Bayesian p-value using Deviance
  output$FT.p.value <- FT.p.value #Bayesian p-value using Freeman-Tukey statistic
  output$deviance.plot <- deviance.plot
  output$FT.plot <- FT.plot
  output$predictive.plots <- predictive.plots
  
  return(output) 
  
} # end of PSkSCR.BAYESIAN.p.value
  
  
  
###############################################################################
###############################################################################


###############################################################################
###############################################################################
# Function: PSkSCR.neg.likelihood.and.tukey  for bayesian p-value calculation
# Input : vector of parameters(x) in logit/log scale, Data,design matrices
#          and offset vectors
# Output : Negative likelihood and Freeman-Tukey statistic 
#          for a one generated data set


PSkSCR.neg.likelihood.and.tukey <- function(x,obs.data,indicator,
                                     captureDM,thetaDM,lambdaDM,p_lossDM,
                                     captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET){
  
  param <- PSkSCR.unpack.parm(x,obs.data,
                              captureDM,thetaDM,lambdaDM,p_lossDM,
                              captureOFFSET,thetaOFFSET,lambdaOFFSET,p_lossOFFSET)  
  N <- param$N
  
  gen.data <- NULL
  # all possible observable capture histories("00" is unobservable)
  gen.data$history <- c(obs.data$history,"OTHER") 
  gen.data$category <- obs.data$category
  
  get.sign  <- c(sign(obs.data$counts),1) # 1 needed to for the sign of "OTHER"
  sample.size <- sum(abs(obs.data$counts))
  

  # the function "PSkSCR.log.prob.history" is in the file "PSkSCR.neg.log.likelihood" 
  # which was used to find the log probabilities
  # call the function 'PSkSCR.log.prob.history'. These log probabilities and the
  # counts are in the exact  order related to each capture history.
  expected.log.prob  <- PSkSCR.log.prob.history(param,indicator)
  
  # expected probabilities
  expected.prob <- exp(expected.log.prob) - (expected.log.prob==0) 
  
  # probability for capture history of unseen animals
  # (i.e. for two sampling occasions "00", for 3 sampling occasion "000" , ..)
  prob.hist.00 <- expected.prob[length(expected.prob)]
  
  # total probability of all other observable but not observed histories
  prob.hist.other <- 1 - ( sum(expected.prob))
  
  if(prob.hist.other < 0){
    print(" Error in probability calculations in Bayesian approach")
    stop
  }
  
  # probability for all  possible observable capture histories (without p.00)
  prob.hist.all <- c( expected.prob[1:(length(expected.prob)-1)], prob.hist.other )
  
  # conditional probabilities for  observable capture histories
  cond.obs.hist.exp.prob <- prob.hist.all/(1-prob.hist.00)  
  
  # generate data
  gen.data$counts <- rmultinom(1, size = sample.size, 
                               prob = cond.obs.hist.exp.prob) * get.sign
  
  
  ct.00 <- N-sum(abs(gen.data$counts))  # count for capture history "00"
  # counts of observable capture histories and capture history "00"
  ct <- c(abs(gen.data$counts),ct.00)
  
  # probabilities of all observable capture histories including 
  # the history "other" and capture history "00"
  prob.hist <-  c(prob.hist.all,prob.hist.00 )
  
  # log- probabilities of observable capture histories and capture history "00"
  log.prob.hist <- log(prob.hist + (prob.hist==0))
  
  #######################################
  # log-ikelihood has 2 parts
  part1 <- lgamma(N+1) - sum(lgamma(ct+1)) # factorial part
  
  part2 <- sum(ct * log.prob.hist)
  
  # log-likelihood
  LLH <-  part1 +  part2  
  
  # negative log-likelihood for generaded data
  NLL.gen <- - LLH  # -(log-likelihood)
  
  
  ######### Freeman-Tukey Statistic  ########
  
  #Find expected counts for generated data
  exp.data.counts <- sample.size * cond.obs.hist.exp.prob * get.sign
  
  ####### Freeman-Tukey Statistic  for simulated data ##################
  
  tukey.statistic.simulated.data <- sum((sqrt(abs(gen.data$counts)) - 
                                           sqrt(abs(exp.data.counts)))^2)
  
  ####### Freeman-Tukey Statistic  for observed data ###################
  
  # obs.counts
  #  This produces the observed counts as same order of the expected capture
  #  history. "0" at the end is to match the counts related to to the "OTHER"
  # in expected data
  obs.counts <- c(obs.data$counts,0)
  
  tukey.statistic.observed.data <- sum((sqrt(abs(obs.counts) )- 
                                          sqrt(abs(exp.data.counts)))^2)
  
  ###########################################################################
  
  return(c(NLL.gen = NLL.gen, obs.FT.ststistic = tukey.statistic.observed.data,
           sim.FT.statistic = tukey.statistic.simulated.data))
  
}  # end of PSkSCR.neg.likelihood.and.tukey

###############################################################################
###############################################################################

