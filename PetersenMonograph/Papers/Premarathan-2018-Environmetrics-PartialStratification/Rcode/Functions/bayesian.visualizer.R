###############################################################################
# Following functions gives the mean and sd for the normal prior distribution
# in logit/log scale (regular form) when the prior distribution values are given
# in regular form(logit/log scale) 

###############################################################################
################## visualizer.prob.to.logit.prob ##############################
###############################################################################
# this shows plots of the prior distributions for probability with the mean and 
# the sd in regular form and logit scale 

# Function : visualizer.prob.to.logit.prob
# Input : mean and sd for the normal distribution for probability in regular form
# output : distribution plots (regular form and logit scale)

visualizer.prob.to.logit.prob <- function(priormean,priorsd ) {
  
  ######### Regular Scale ######################
  
  # generate normal random number. Since these are probabilities, consider only the 
  # between 0 and 1 when plotting the distribution and converting to logit values
  #### Warning messages will be given when there are values outside the range [0,1]###
  regular = data.frame(rnorm(20000,mean=priormean, sd=priorsd))
  names(regular) <- c("regularval")
  
  # prior distribution for the probability in regular scale values
  regular.plot <- ggplot( regular,aes(x = regularval))+
                    geom_density()+
                    ggtitle("Prior Distribution in Regular Scale")+
                    annotate("text", label = paste("Regular Values :",
                                   "\n mean = ",round(priormean,3), 
                                   "\n sd   = ",round(priorsd,3),
                                   sep=""),
                             x=-Inf,hjust = 0, 
                             y = Inf, vjust = 2, 
                             color = "darkred",fontface = "bold")+
                    xlab("Probability ") +
                    ylab("Density")+
                    xlim(0, 1)
  
    
  ######### Logit Scale  ######################
  
  logitnorm <-  data.frame(log(regular/(1-regular))) # Convert  to logit values
  names(logitnorm) <- c("logitval")
  # mean  in logit scale
  logit_prior_mean <- mean(logitnorm[,1],na.rm = TRUE)  # 
  # sd in logit scale
  logit_prior_sd  <- sd(logitnorm[,1],na.rm = TRUE)
  
  # prior distribution for the probability in logit form
  logit.plot <- ggplot( logitnorm,aes(x = logitval))+ 
                  geom_density()+
                  ggtitle("Prior Distribution in Logit Scale")+
                  annotate("text", label = paste("Logit Values :",
                                   "\n mean = ",round(logit_prior_mean,2),
                                   "\n sd   = ",round(logit_prior_sd,2),
                                   sep=""),
                            x =-Inf, hjust = 0,
                            y = Inf, vjust = 2, 
                            color = "darkred",fontface = "bold")+
                  xlab("Probability in logit scale") +
                  ylab("Density")
  
 
  ################################# 
  
  # Prior Distribution both in regular form and logit scale( 2 plots)
  prior.prob.plots <- arrangeGrob(regular.plot, logit.plot, nrow = 1, 
                        main = textGrob(
                                        paste("Prior Distribution for ",
                                              "Probability in Regular Scale and ",
                                              "Logit Scale",sep=""),
                                        just = "top", 
                                        vjust = 0.75, 
                                        gp = gpar(fontface = "bold")))
  
  return(prior.prob.plots)
  
}  # end visualizer.prob.to.logit.prob
  
###############################################################################  
################## visualizer.N.to.log.N ######################################
###############################################################################
# this shows plots of the prior distributions for population size(N) with the
# mean and the sd in regular form and log scale 

# Function : visualizer.N.to.log.N
# Input : mean and sd for normal distribution for N in regular form
# output : distribution plots (regular form and log scale)
  

visualizer.N.to.log.N <- function(priormean,priorsd ) {
  
  ######### Regular Scale ######################
  
  # generate normal random number. Since these are population sizes, consider only 
  # the values greater than 0 because since we consider the log values of the 
  # population size
  ### Warning messages will be given when there are values smaller than 0 ###
  regular = data.frame(rnorm(20000,mean=priormean, sd=priorsd))
  names(regular) <- c("regularval")
  
  # prior distribution for the N in regular form
  regualr.plot <- ggplot( regular,aes(x = regularval))+
                    geom_density()+
                    ggtitle("Prior Distribution in Regular Scale")+
                    annotate("text", label = paste("Regular Values :",
                                   "\n mean = ",prettyNum(priormean, big.mark = ","), 
                                   "\n sd   = ",prettyNum(priorsd, big.mark = ",")),
                              x =-Inf, hjust = 0,
                              y = Inf, vjust = 2, 
                              color = "darkred",fontface = "bold")+
                    xlab("N") + 
                    ylab("Density")+ 
                    xlim(max(1,min(regular[,1])),max(regular[,1]))
    
                    
  
  ######### Log Scale  ######################
  
  lognorm <-  data.frame(log(regular)) # Convert  to log values
  names(lognorm) <- c("logval")
  # mean  in log scale
  log_prior_mean <- mean(lognorm[,1],na.rm = TRUE)  # 
  # sd in log scale
  log_prior_sd  <- sd(lognorm[,1],na.rm = TRUE)
  
  # prior distribution for the probability in log scale values
  log.plot <- ggplot( lognorm,aes(x = logval))+ 
                geom_density()+
                ggtitle("Prior Distribution in Log Scale")+
                annotate("text", label = paste("Log Values :", 
                                   "\n mean = ",round(log_prior_mean,2), 
                                   "\n sd   = ",round(log_prior_sd,2), sep=""),
                           x =-Inf,hjust = 0,
                           y = Inf, vjust = 2, 
                           color = "darkred",fontface = "bold")+
                xlab("N in log scale") +
                ylab("Density")
                  
  
  #######################################
  
  # Prior Distribution both in regular form and log scale( 2 plots)
  prior.N.plots <- arrangeGrob(regualr.plot, log.plot, nrow = 1, 
                                  main = textGrob(
                                    paste("Prior Distribution for ",
                                          "Population Size(N) in Regular ", 
                                          "Scale and Log Scale",sep=""),
                                    just = "top", 
                                    vjust = 0.75, 
                                    gp = gpar(fontface = "bold")))
  
  return(prior.N.plots)
  
}  # end visualizer.N.to.log.N



###############################################################################
###################### visualizer.lambda.to.logit.lambda ######################
###############################################################################
# this shows plots of the prior distributions for category proportions (lambda)
# with the mean and the sd in regular form and logit scale 

# Function : visualizer.lambda.to.logit.lambda
# Input : vector of shape parameters for Dirichlet distribution  in regular form
#         and vector of categories
# output : distribution plots (regular form and logit scale)

visualizer.lambda.to.logit.lambda <- function(dirichlet.param,cat ) {

  if(length(dirichlet.param) != length(cat)){
    stop(" Error : Number of shape parameters for Dirichlet distribution
        are different from number of categories")
  }
  
  ######### Regular Scale ######################
  
  # get the proportion for each category
  reg.prob <- dirichlet.param/sum(dirichlet.param)
  # if we need to keep the categories as the same order as it is defined in cat
  # then use following as x=categories in qplot
  #categories <- factor(cat, levels=cat) 
  
  # generate Dirichler random variables to get the sd
  nrv <- 10000  # number of Dirichlet random variables
  dirichlet_rv <- rdirichlet(nrv,dirichlet.param)
  priormean <-colMeans(dirichlet_rv)
  priorsd <- apply(dirichlet_rv, 2, sd)
  
  regualr.lambda.plot <- qplot(x=cat,y=reg.prob, 
                               main = "Priors for lambda in Regular Scale \n(values for Dirichlet Distribution)") +
                         geom_point(colour = "darkred", size = 4)+
                         annotate("text", label = paste(
                                   "\nmean = ",round(reg.prob,2) , 
                                   "\nsd = ",round(priorsd, 3)), 
                                  x =cat,
                                  hjust = 0.5,y = max(reg.prob), 
                                  vjust = 2, color ="darkred",size=4)+                
                         xlab("Category") +
                         ylab("Proportion")

  
  regular <- as.vector(dirichlet_rv)
  Cat <- as.factor(rep(cat, each=nrv))
  regular.df <- data.frame(regular, Cat) # convert to data frame
  
  regualr.lambda.density <- ggplot(regular.df, aes(regular,colour=Cat,
                                                   linetype=Cat)) +
                              geom_density()+
                              ggtitle("Prior Distributions in Regular Scale")+
                              xlab("Probability") +
                              ylab("Density")
 

  ######### Logit Scale  ########################
  # shape parameters for Dirichlet distribution in logit scale
  logit.prob <- logit(reg.prob)
  
  # generated Dirichler random variables in logit scale
  logit_rv <- logit(dirichlet_rv)
  
  # sd in logit scale
  logit_prior_sd <- apply(logit_rv, 2, sd)
  
  # prior distribution for the probability in logit scale values
  logit.lambda.plot <- qplot(x=cat,y=logit.prob, 
                         main ="Priors for lambda in Logit Scale \n(values for Normal Distribution)") +
                       geom_point(colour = "darkred", size = 4)+
                       annotate("text", label = paste( 
                                   "\nmean = ",round(logit.prob,2), 
                                   "\nsd = ",round(logit_prior_sd,3), sep=""),
                                x =cat,
                                hjust = 0.5,y = max(logit.prob), 
                                vjust = 2, color ="darkred",size=4)+
                       xlab("Category") + 
                       ylab("Proportion in Logit Scale")
  
  
  logitvals <- as.vector(logit_rv)
  logitvals.df <- data.frame(logitvals, Cat)
  
  logit.lambda.density <- ggplot(logitvals.df, aes(logitvals,colour=Cat,
                                                   linetype=Cat)) +
                            geom_density()+
                            ggtitle("Prior Distributions in Logit Scale")+
                            xlab("Probability in logit scale") +
                            ylab("Density")
  
  
  ###########################################  
  # Prior Distribution both in regular form and log scale( 2 plots)
  prior.lambda.plots <- arrangeGrob(regualr.lambda.plot, logit.lambda.plot,
                                    regualr.lambda.density,logit.lambda.density,
                                    nrow = 2,
                                    main = textGrob(
                                       paste("Dirichlet Priors for category ",
                                             "Proportions(lambda) in Regular ", 
                                             "Scale and Logit Scale",sep=""),
                                    just = "top", 
                                    vjust = 0.75, 
                                    gp = gpar(fontface = "bold")))
  
  return(prior.lambda.plots)
  
}  # end visualizer.lambda.to.logit

###############################################################################
################ visualizer.logit.to.regular.prob #############################
###############################################################################

# this shows plots of the prior distributions for probability with the mean and 
# the sd in regular form and logit scale  when the prior is given in logit scale

# Function : visualizer.logit.to.regular.prob
# Input : mean and sd for the normal distribution for probability in logit scale
# output : distribution plots (regular form and logit scale)

visualizer.logit.to.regular.prob <- function(logit_priormean,logit_priorsd){  
  
  ######### logit Scale ######################
  
  # normal random numbers in logit form
  logitnorm = data.frame(rnorm(10000,mean=logit_priormean, sd=logit_priorsd))
  names(logitnorm) <- c("logitval")
  
  # prior distribution for the probability in logit scale values
  logit.plot <- ggplot( logitnorm,aes(x = logitval))+
                 geom_density()+
                 ggtitle("Prior Distribution in Logit Scale")+
                 annotate("text", label = paste("Logit Values :",
                                                "\n mean = ",round(logit_priormean,3), 
                                                "\n sd   = ",round(logit_priorsd,3),
                                                sep=""),
                          x=-Inf,hjust = 0,
                          y = Inf, vjust = 2, 
                          color = "darkred",fontface = "bold")+
                  xlab("Probability in logit scale") +
                  ylab("Density")
  

  ######### Regular Scale ######################
  
  regular <-  data.frame(1/( 1+ exp(-logitnorm))) # Convert  to logit values
  names(regular) <- c("regularval")
  # mean  in regular scale
  regular_prior_mean <- colMeans(regular)  # 
  # sd in regular scale
  regular_prior_sd  <- sd(regular[,1])
  
  # prior distribution for the probability in regular form
  regualr.plot <- ggplot( regular,aes(x = regularval))+ 
                    geom_density()+
                    ggtitle("Prior Distribution in Regular Scale")+
                    annotate("text", label = paste("Regular Values :",
                                                   "\n mean = ",round(regular_prior_mean,3),
                                                   "\n sd   = ",round(regular_prior_sd,3),
                                                   sep=""),
                             x =-Inf, hjust = 0,
                             y = Inf, vjust = 2, 
                             color = "darkred",fontface = "bold")+
                    xlab("Probability") +
                    ylab("Density")
    
  ####################################
  
  # Prior Distribution both in regular form and logit scale( 2 plots)
  prior.prob.plots <- arrangeGrob(logit.plot,regualr.plot, nrow = 1, 
                                  main = textGrob(
                                    paste("Prior Distribution for ",
                                          "Probability in Logit Scale and ",
                                          "Regular Scale",sep=""),
                                    just = "top", 
                                    vjust = 0.75, 
                                    gp = gpar(fontface = "bold")))
  
  return(prior.prob.plots)
  
}  # end of visualizer.logit.to.regular.prob


###############################################################################
################ visualizer.log.N.to.N ########################################
###############################################################################

# this shows plots of the prior distributions for Population size with the mean
# and the sd in regular form and log scale  when the prior is given in log scale

# Function : visualizer.log.N.to.N
# Input : mean and sd for the normal distribution for probability in log scale
# output : distribution plots (regular form and log scale)

visualizer.log.N.to.N <- function(log_priormean,log_priorsd){  
  
  ######### logit Scale ######################
  
  # normal random variables in log form
  lognorm = data.frame(rnorm(20000,mean=log_priormean, sd=log_priorsd))
  names(lognorm) <- c("logval")
   
  # prior distribution for the probability in log scale values
  log.plot <- ggplot( lognorm,aes(x = logval))+ 
                geom_density()+
                ggtitle("Prior Distribution in Log Scale") +
                annotate("text", label = paste("Log Values :", 
                                   "\n mean = ",round(log_priormean,2), 
                                   "\n sd   = ",round(log_priorsd,2), sep=""),
                         x =-Inf, hjust = 0,
                         y = Inf, vjust = 2, 
                         color = "darkred",fontface = "bold")+
                xlab("Population Size (N) in log scale") +
                ylab("Density")
  
  
  ######### Regular Scale ######################
  
  regular <-  data.frame(exp(lognorm)) # Convert  to log values
  names(regular) <- c("regularval")
  # mean  in regular scale
  regular_prior_mean <- colMeans(regular)  # 
  # sd in regular scale
  regular_prior_sd  <- sd(regular[,1])
  
  # prior distribution for the probability in regular form
  regualr.plot <- ggplot( regular,aes(x = regularval))+ 
                    geom_density()+
                    ggtitle("Prior Distribution in Regular Scale")+
                    annotate("text", label = paste("Regular Values :",
                                      "\n mean = ",prettyNum(round(regular_prior_mean),
                                                             big.mark = ","),
                                      "\n sd   = ",prettyNum(round(regular_prior_sd),
                                                             big.mark = ","), sep=""),
                              x =-Inf, hjust = 0,
                              y = Inf, vjust = 2, 
                              color = "darkred",fontface = "bold")+
                    xlab("Population Size (N)") +
                    ylab("Density")
  
  ####################################
  
  # Prior Distribution both in regular form and logit scale( 2 plots)
  prior.prob.plots <- arrangeGrob(log.plot,regualr.plot, nrow = 1, 
                                  main = textGrob(
                                    paste("Prior Distribution for ",
                                          "Population Size(N) in Log Scale and ",
                                          "Regular Scale",sep=""),
                                    just = "top", 
                                    vjust = 0.75, 
                                    gp = gpar(fontface = "bold")))
  
  return(prior.prob.plots)
  
}  # end of visualizer.log.N.to.N

###############################################################################
###############################################################################


###############################################################################
###################### Testing ################################################

# visualizer.prob.to.logit.prob(priormean=0.79,priorsd=0.2)
# visualizer.N.to.log.N (priormean=250000,priorsd=45000)
# visualizer.lambda.to.logit.lambda(c(10,5,3),c("B","A","C"))

# visualizer.logit.to.regular.prob(logit_priormean=1.32,logit_priorsd=0.12)
# visualizer.log.N.to.N(log_priormean=12.41,log_priorsd=.19)




###############################################################################
############ Following are used to find the values for the priors  ############   
############ in logit and log scale                                ############
###############################################################################
# rdirichlet(10, c(1,2))
# 
# ########################
# # logit value of regular probability
# logit(0.3298)
# 
# ## normal distribution on the logit scale can be made to fit very closely
# # to the chosen beta distribution.  When mean=0 and sd =1.78 this is 
# # approximately uniform distribution
# logitnorm = data.frame(rnorm(100000,mean=0, sd=3))
# names(logitnorm) <- c("val")
# ggplot( logitnorm,aes(x = val))+geom_density()
# 
# head(logitnorm)
# 
# regular <-  data.frame(1/( 1+ exp(-logitnorm)))
# colMeans(regular)
# sd(regular[,1])
# names(regular) <- c("val")
# ggplot( regular,aes(x = val))+ geom_density()

#########################
# 
# # log value of the regular population size
# log(205000)
# 
# ## normal distribution in log scale
# 
# lognorm = data.frame(rnorm(100000,mean=12.5, sd=0.2))
# names(lognorm) <- c("val")
# regular <-  data.frame(exp(lognorm))
# mean(regular[,1])
# sd(regular[,1])
# names(regular) <- c("val")
# 
# ggplot( regular,aes(x = val))+ geom_density()
# #########################
###############################################################################

  
  