###############################################################################
##########  Power of the study and verify the power using simulations #########
###############################################################################

# generate some data where the probabilities vary between sexes/sample times and 
# fit a model using null and alternative hypothesis. then use the power of the 
# study to detect such differences Verify the power using simulation study



###############################################################################
############################################################################### 
# Function power.study
# inputs : given.param(capture probabilities, category proportions, 
#      subsample proportions),possible categories
#      model id's and design matrices and offset vectors for unrestricted model, 
#      model id's and design matrices and offset vectors for restricted model, 
#      and sample size,number of simulations, alpha
#
# outputs : return the power of the study  for the generated data and the the
#           histogram from the  simulations to verify the power

power.study.simulation = function(given.param,categories,
                           unrest.model.id,unrest.captureDM,unrest.thetaDM,
                                  unrest.lambdaDM,unrest.captureOFFSET,
                                  unrest.thetaOFFSET,unrest.lambdaOFFSET,
                           rest.model.id,rest.captureDM,rest.thetaDM,
                                  rest.lambdaDM,rest.captureOFFSET,
                                  rest.thetaOFFSET,rest.lambdaOFFSET,
                           sample.size,n.simulation,alpha){

  cap.histories = c(paste(categories,"0",sep=""),paste(categories,categories,sep=""),
                         paste("0",categories,sep=""))

  simulate.Data = NULL
  # all possible observable capture histories("00" is unobservable)
  simulate.Data$History = c("U0", "UU", cap.histories,"0U") 
  
  simulate.Data$category = categories

  # add "00" capture history to end of the vector of capture history
  hist = c(simulate.Data$History,"00")

  # the function "log.prob.history" is in the file "neg.log.likelihood" which was
  # used to find the log probabilities of the capture histories
  sim.cats = simulate.Data$category
  log_prob_hist = log.prob.history(given.param,hist,sim.cats)  
  # These log probabilities and the counts are in the exact order related to 
  # each capture history. 

  prob_hist = exp(log_prob_hist) - (log_prob_hist==0) #probabilities for histories

  #  probabilities for observable capture histories( without history"00")
  obs.prob.hist  =  prob_hist[1:(length(prob_hist)-1)] 
  
  p_00 = prob_hist[length(prob_hist)]#  probability of the  capture history "00"
  #conditional probabilities for  observable capture histories
  cond.obs.prob.hist = obs.prob.hist/(1-p_00)  

  ## generate data as expected values. Didn't use the parameter value N    ###
  ## since it gives  different sample sizes for simulated data sets. Used  ###
  ## theconditional probabilities to generate the data. Then generate data ###
  ## as expected values                                                    ###
  
  # generate data using the sample size, and the get the expected values
  simulate.Data$counts=sample.size*cond.obs.prob.hist 
    
  ### Now analyse these expected data exactly as if they were real data  ####
  ### (to find the noncentrality parameter)                              ####
  
  ###### Fit alternative model to the simulated data ###################
  MLE.alternative.model = fit.model(model.id=unrest.model.id,simulate.Data,
                                      captureDM=unrest.captureDM,
                                        thetaDM=unrest.thetaDM,
                                        lambdaDM=unrest.lambdaDM, 
                                      captureOFFSET=unrest.captureOFFSET,
                                        thetaOFFSET=unrest.thetaOFFSET,
                                        lambdaOFFSET=unrest.lambdaOFFSET)
  
  #print results(Model info,MLEs, SE,residual plot, etc..)
  print.output(MLE.alternative.model) 
  
  #deviance part and number of parameters in the alternative  model
  alternative.model.deviance = -2* (-MLE.alternative.model$NLL)
  alternative.model.np = MLE.alternative.model$np 
  
 
  ######  Fit null hypothesis model to the simulated data ###################
  MLE.null.model = fit.model(model.id=rest.model.id,simulate.Data,
                              captureDM=rest.captureDM,thetaDM=rest.thetaDM,
                                lambdaDM=rest.lambdaDM,
                              captureOFFSET=rest.captureOFFSET,
                                thetaOFFSET=rest.thetaOFFSET,
                                lambdaOFFSET=rest.lambdaOFFSET)
  
  # print results (Model info,MLEs, SE, residual plot, etc..)
  print.output(MLE.null.model)
  
  #deviance part and number of parameters in the null model
  null.model.deviance = -2* (-MLE.null.model$NLL)
  null.model.np = MLE.null.model$np   
  
  
  ######## Deviance ################
  deviance.expected.data = null.model.deviance - alternative.model.deviance
  
  ##degree of freedom is the difference of the no:parameters of the two models 
  degree.of.freedom = alternative.model.np - null.model.np  
  
  
  ############### calculate noncentrrsality parameter ######################### 
  ############### using Devineau (2006) and Lebreton(1992) papers #############  
  # a simple way exist to compute the asymptotic value for chi-square noncentral
  # parameter under the specific alternative model The chi square "test statistic"
  # of H0(null) vs. Ha(alternative) hypothesis computed from such expected 
  # data(under Ha) is really the approximate value of noncentral parameter.
  
  noncentral.chisq.factor = deviance.expected.data # noncentrality factor 

  # threshold of the chi-square test at level alpha 
  theshold.alpha = qchisq(1-alpha,degree.of.freedom)  

  
  ############# Calculate the  power of the study #############################
  # power  = P( chi-square(L,d) > T(alpha) )
    #  L is the noncentrality factor 
    #  d is the difference of the degrees of freedom of the two models 
    #  T(alpha) is the threshold of the chi-square test at level alpha
  
  power = 1-pchisq(theshold.alpha, degree.of.freedom, ncp=noncentral.chisq.factor)   
  
  ###### Verify the power using simulation study ###############
  
  # devince for  "n.simulation" number of simulated samples
  sim.deviance = rdply(n.simulation,
                         simulate.data.deviance(simulate.Data,sample.size,
                                                cond.obs.prob.hist)) 
  names(sim.deviance) = c("sample","DEV")
  
  # calculate power from simulated data
  power.verify = mean(sim.deviance$DEV > theshold.alpha) 
  
  ################### create Histograms ######################################
  # x position to show the p-value in the graph
  x.pos = min(sim.deviance$DEV,theshold.alpha,deviance.expected.data) + 
             .3*(max(sim.deviance$DEV,theshold.alpha,deviance.expected.data) -
                    min(sim.deviance$DEV,theshold.alpha,deviance.expected.data)) 
  
  # mark the deviance for observed data
  cuts1 = data.frame(Values=c("Threshold", "Deviance(gen data)"), 
                            vals=c(theshold.alpha,deviance.expected.data)) 
  
  sim.deviance.histogram = ggplot(sim.deviance, aes(x=DEV,..density.. ))+   
                              geom_histogram(aes(y=..density.. ))+
                              geom_vline(data=cuts1,aes(xintercept=vals,
                                                        linetype=Values,colour=Values),
                                         size=1,show_guide = TRUE)+
                              xlab("Deviance") + ylab("Density")+
                              ggtitle(paste("Verify power using simulations", 
                                        "\n Model under H0 : ", rest.model.id,
                                        "\n Model under Ha : ", unrest.model.id,sep=" "))+
                              theme(axis.text=element_text(size=10,face="bold"),
                                    axis.title=element_text(size=14,face="bold"),
                                    plot.title = element_text(size = 14,face="bold"))+
                              annotate("text", label = paste("Power (generated data) = ",
                                                        round(power,3),
                                                        "\nPower (using simulation ) =",
                                                        power.verify,sep=""), 
                                       x = x.pos, hjust = 0, y = Inf, 
                                       vjust = 2, color = "darkred")
                         
  
  ## return the power of the study  for the generated data and the the histogram
  ## of the  simulation to verify the power 
  res = NULL
  res$power.Devineau = power  #power of the test for generated data(Devineau method)
  res$power.simulation = power.verify # veriry the power using simulated data
  res$sim.deviance.histogram = sim.deviance.histogram
  #res$sim.deviance.histogram.2 = sim.deviance.histogram.2
  
  
  return(res)
  
}  # end of power.study.simulation
  
  
###############################################################################
################## deviance for a one simulated sample data ###################

# Function : simulate.data.deviance

# input : Generated data set , sample size, conditional probabilities for 
#         observable capture histories to simulate data
# output : return the deviance for a simulated data set

simulate.data.deviance = function(simulate.Data, sample.size, simulate.prob){
  
  sim.data = NULL
  sim.data$History = simulate.Data$History
  sim.data$category = simulate.Data$category
  
  sim.data$counts = rmultinom(1, size = sample.size, prob = simulate.prob)
  
  ### fit the unrestricted(alternative)  model  for the simulated data set
  MLE.alternative.model.sim =  fit.model(model.id=unrest.model.id,sim.data,
                                         captureDM=unrest.captureDM,
                                           thetaDM=unrest.thetaDM,
                                           lambdaDM=unrest.lambdaDM,
                                         captureOFFSET=unrest.captureOFFSET,
                                           thetaOFFSET=unrest.thetaOFFSET,
                                           lambdaOFFSET=unrest.lambdaOFFSET)
   
  sim.alternative.model.deviance = -2* (-MLE.alternative.model.sim$NLL)
  
  
  ### fit the restricted(null)  model for the simulated data set
  MLE.null.model.sim =  fit.model(model.id=rest.model.id,sim.data,
                                  captureDM=rest.captureDM,thetaDM=rest.thetaDM,
                                     lambdaDM=rest.lambdaDM,
                                  captureOFFSET=rest.captureOFFSET,
                                     thetaOFFSET=rest.thetaOFFSET,
                                     lambdaOFFSET=rest.lambdaOFFSET)
  
  sim.null.model.deviance = -2* (-MLE.null.model.sim$NLL)
  
  sim.data.deviance = sim.null.model.deviance - sim.alternative.model.deviance
    
  return(sim.data.deviance)
      
} #end of  simulate.data.deviance



###############################################################################
###########  Test power.study function ########################################
###############################################################################
# This is not working. this was for the earliear version
# (change this according to the power.study.simulation function)
#
#  # give the  categories in the population
#  categories = c("M","F")
#     
#  # give the different  capture probabilities for p1M, p1F, p2M and p2F
#  cap.prob = c( 0.08, 0.08, 0.06, 0.06)
#     
#  # give the category proportions( total should add up to 1)
#  lambda= c(0.6, 0.4)
#     
#  # give the subsample proportions for the time 1 and 2
#  theta = c(0.8, 0.5)
#     
#  # total numbers individuals captured for the study 
#  sample.size = 10000
#     
#  n.simulation = 10 #number of simulations to verify the power of the study.
#     
#  ######  unrestricted model #####################################
#  #model identification : unrestricted model
#  unrestricted.model.id = paste("( p(c*t), theta(t), lambda(c) )")
#     
#  # give  the required design matrices
#  # Design matrix for capture recapture probabilities
#  unrestricted.captureDM = create.DM(c(1,2,3,4)) 
#
#  # Design matrix for theta(sampling(sexing) fractions)
#  thetaDM   = create.DM(c(1,2)) 
#  lambdaDM  = create.DM(c(1))  # Design matrix for lambda(Category proportion)
#     
#  #give the offset vectors(vectors of zeros should be given since no restriction)
#  captureOFFSET = c(0,0,0,0) 
#  thetaOFFSET   = c(0,0)
#  lambdaOFFSET  = c(0)
#     
#  ##### restricted model ##########################################
#     
#  #model identification : restricted model(capture probabilities equal)
#  restricted.model.id = paste("( p(.), theta(t), lambda(c) )")
#     
#  # Design matrix for capture recapture probabilities for restricted model
#  #  (capture probabilities equal)
#  restricted.captureDM = create.DM(c(1,1,1,1)) 
#     
#     
#  result = power.study(categories,cap.prob,lambda,theta,sample.size,
#                        unrestricted.model.id,unrestricted.captureDM,
#                           thetaDM,lambdaDM,captureOFFSET,thetaOFFSET,lambdaOFFSET,
#                        restricted.model.id,restricted.captureDM,n.simulation)
#     

################################################################################
################################################################################

