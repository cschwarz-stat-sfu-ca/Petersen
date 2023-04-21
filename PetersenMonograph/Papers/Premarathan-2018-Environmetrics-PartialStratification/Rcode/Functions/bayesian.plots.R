#############################################################################
# Function : bayesian.plots 
# Input : results from the fit model (3 dimensional array of thinned samples,
#         and the thinned matrix and other info)
# Output : trace plots, Posterior density plots, ACF plots, and 
#          Gelman Rubin Diagnose plots with Rhat values


bayesian.plots <- function(Result){
  
  #############################################################################
  ####################### trace plots #########################################

  thinned.matrix <- Result$MCMC$summary.tables$thinned.matrix
  n.samples <- length(Result$thinned.array[,1,1])
  x= 1:n.samples  # number of posterior samples for each chain
  n.beta <- Result$np # number of beta parameters
  ncats <- length(Result$Data$category) # number of categories
  param_names <- Result$MCMC$summary.tables$param_names
  
  chains.names <- paste("chain_", (1:Result$n.chains), sep="") # chain names
  
  chain<- c() # vector named "chain" store the names name of each chain n.samples times
  for(i in 1:Result$n.chains){chain <- c(chain,
                                         rep(chains.names[i],
                                             Result$n.post.burnin/Result$n.thin))}
  
  
  sample <- rep(x,Result$n.chains ) # replicate numbers 1 to n.samples 
  thinned.df <- data.frame(cbind(thinned.matrix,sample)) # data frame 
  thinned.df<- cbind(thinned.df,chain)
  
  
  plots <- llply((1:(n.beta+3*ncats+3+ncats)),
                 function(k){qplot(x=thinned.df$sample,y=thinned.df[,k],
                                   colour=chain,
                                   geom = "line", main=param_names[k])+
                               xlab("Samples") + ylab("Value")+
                               theme(legend.position="none")}) 
  
  
  ## trace plots for beta parameters
  beta.trace.plots<- do.call(arrangeGrob, c(plots[1:n.beta], ncol=2,
                             main=paste("Trace plots for Beta parameters",
                                        "\n(values are in logit/log scale)")))
  
  ## trace plots for capture probabilities
  capture.trace.plots <- do.call(arrangeGrob, 
                                 c(plots[(n.beta+1):(n.beta+2*ncats)], 
                                   nrow=ncats,
                    main=paste("Trace plots for Capture probabilities")))
  
  ## trace plots for category proportions(lambda)
  lambda.trace.plots <- do.call(arrangeGrob, 
                                c(plots[(n.beta+1+2*ncats):(n.beta+3*ncats)], 
                                  ncol=1,
                          main=paste("Trace plots for Category Proportions")))
  
  ## trace plots for sub sample proportions(theta)
  theta.trace.plots <- do.call(arrangeGrob,
                               c(plots[(n.beta+1+3*ncats):(n.beta+3*ncats+2)], 
                                 ncol=1,
                       main=paste("Trace plots for sub-sample Proportions")))
  
  ## trace plots for sub Population size and category total
  pop.trace.plots <- do.call(arrangeGrob, 
                             c(plots[(n.beta+3*ncats+3):(n.beta+3*ncats+3+ncats)], 
                               ncol=1,
            main=paste("Trace plots for Population Size and Category Total ")))
  
  
  ######################### End of Trace plots ################################
  #############################################################################
  
  #############################################################################
  ############################### Posterior plots #############################
  
  alpha <- Result$alpha 
  CI_vals <- c(alpha/2, (1-alpha/2)) # credible interval values
  #credible intervals for all beta, and regular parameters
  CI.df <- adply(thinned.matrix,2,quantile,probs = CI_vals)
  #posterior means
  posterior.mean <- Result$MCMC$summary.tables$summary_stat_thin[,2]
  #posterior sd
  posterior.sd <- Result$MCMC$summary.tables$summary_stat_thin[,3]
  
   
  ## posterior plots for beta parameters
  density.plots<-llply(1:n.beta, 
                       function(k){qplot(x=thinned.df[,k],geom="density",
                                         main=paste("Distribution of ",
                                                    param_names[k],sep=""))+
                                    xlab(param_names[k]) + 
                                    ylab("Density")+
                                    geom_vline(xintercept=c(CI.df[k,2],
                                                            CI.df[k,3]))+
                                    annotate("text",label = paste("\nmean = ",
                                                  round(posterior.mean[k],3),
                                                  "\nsd = ",round(posterior.sd[k],3),
                                                  "\n v-lines: ",(1-alpha)*100,"% CI",
                                                  sep=""),
                                              x =min(thinned.df[,k]), hjust = 0,
                                              y = Inf,vjust = 1,size=4)})
  #stat_function(fun = dnorm, colour = "red")
  
  beta.posterior.plots <- do.call(arrangeGrob, c(density.plots, ncol=2,
                            main=paste("Posterior plots for Beta parameters")))
  
  
  ## posterior plots for Capture probabilities
  density.plots<-llply((n.beta+1):(n.beta+2*ncats), 
                       function(k){qplot(x=thinned.df[,k],geom="density",
                                         main=paste("Distribution of ",
                                              param_names[k],sep=""))+
                                    xlab(param_names[k]) + 
                                    ylab("Density")+
                                    geom_vline(xintercept=c(CI.df[k,2],
                                               CI.df[k,3]))+
                                    annotate("text", label = paste("\nmean = ",
                                                  round(posterior.mean[k],3),
                                                  "\nsd = ",round(posterior.sd[k],3),
                                                  "\n v-lines: ",(1-alpha)*100,"% CI",
                                                  sep=""),
                                              x =min(thinned.df[,k]), hjust = 0,
                                              y = Inf,vjust = 1,size=4)})
  
  capture.posterior.plots <- do.call(arrangeGrob, c(density.plots, ncol=2,
                     main=paste("Posterior plots for Capture Probabilities")))
  
  
  
  ## posterior plots for Category Proportions (lambda)
  density.plots<-llply((n.beta+2*ncats+1):(n.beta+3*ncats), 
                       function(k){qplot(x=thinned.df[,k],geom="density",
                                         main=paste("Distribution of ",
                                                    param_names[k],sep=""))+
                                     xlab(param_names[k]) + 
                                     ylab("Density")+
                                     geom_vline(xintercept=c(CI.df[k,2],
                                                             CI.df[k,3]))+
                                     annotate("text", label = paste("\n mean = ",
                                                    round(posterior.mean[k],3),
                                                    "\n sd = ",round(posterior.sd[k],3),
                                                    "\n v-lines: ",(1-alpha)*100,"% CI",
                                                    sep=""),
                                                x =min(thinned.df[,k]), hjust = 0,
                                                y = Inf,vjust = 1,size=4)})
  
  lambda.posterior.plots <- do.call(arrangeGrob, c(density.plots, ncol=1,
                      main=paste("Posterior plots for Category Proportions")))
  
  
  
  ## posterior plots for sub-sample Proportions (theta)
  density.plots<-llply((n.beta+3*ncats+1):(n.beta+3*ncats+2), 
                       function(k){qplot(x=thinned.df[,k],geom="density",
                                         main=paste("Distribution of ",
                                                    param_names[k],sep=""))+
                                     xlab(param_names[k]) + 
                                     ylab("Density")+
                                     geom_vline(xintercept=c(CI.df[k,2],
                                                             CI.df[k,3]))+
                                     annotate("text", label = paste("\nmean = ",
                                                     round(posterior.mean[k],3),
                                                     "\nsd = ",round(posterior.sd[k],3),
                                                     "\n v-lines: ",(1-alpha)*100,"% CI",
                                                     sep=""),
                                                x =min(thinned.df[,k]), hjust = 0,
                                                y = Inf,vjust = 1,size=4)})
  
  theta.posterior.plots <- do.call(arrangeGrob, c(density.plots, ncol=1,
                      main=paste("Posterior plots for sub-sample  Proportions")))
  
  
  ## posterior plots for Pupulation Size and Category Total
  density.plots<-llply((n.beta+3*ncats+3):(n.beta+3*ncats+3+ncats), 
                       function(k){qplot(x=thinned.df[,k],geom="density",
                                         main=paste("Distribution of ",
                                                    param_names[k],sep=""))+
                                     xlab(param_names[k]) + 
                                     ylab("Density")+
                                     geom_vline(xintercept=c(CI.df[k,2],
                                                             CI.df[k,3]))+
                                     annotate("text", label = paste("\nmean = ",
                                                    prettyNum(round(posterior.mean[k]),
                                                              big.mark = ","),
                                                    "\nsd = ",
                                                    prettyNum(round(posterior.sd[k]),
                                                              big.mark = ","),
                                                    "\n v-lines: ",(1-alpha)*100,"% CI",
                                                    sep=""),
                                              x =min(thinned.df[,k]), hjust = 0,
                                              y = Inf,vjust = 1,size=4)})
  
  pop.posterior.plots <- do.call(arrangeGrob, c(density.plots, ncol=1,
         main=paste("Posterior plots for Population Size and Category Total")))
  
  ######################## End of Posterior plots #############################
  #############################################################################
  
  #############################################################################
  ############################### ACF plots ###################################
  
  # function to create a ACF plot
  # Input: column number of the data frame for the thnned.matrix
  # output : ACF plot
  acf.function <- function(k){
    acf.out <- acf(thinned.df[, k], plot = FALSE)  # compute ACF without plot
    acf.df <- with(acf.out, data.frame(lag, acf))  # turn into data frame
    acf.plot <- qplot( x = acf.df$lag, y = acf.df$acf) +
      geom_segment(xend = acf.df$lag, yend = 0, size = 1)+
      ggtitle(param_names[k])+
      xlab("Lag") +ylab("ACF")
    return(list(acf.plot))
  }
  
  acf.plots <- list()
  for(i in 1:(n.beta+3*ncats+3+ncats)){
    acf.plots[i] <- acf.function(i)  # create list of ACF plots
  }
  
  
  ## ACF plots for beta parameters
  beta.acf.plots<- do.call(arrangeGrob, c(acf.plots[1:n.beta], ncol=2,
                        main=paste("Autocorrelation for Beta parameters")))
  
  ## ACF plots for Capture Probabilities
  if(ncol(Result$captureDM) > 0){
    capture.acf.plots<- do.call(arrangeGrob, c(acf.plots[(n.beta+1):(n.beta+2*ncats)],
                                            ncol=2,
                          main=paste("Autocorrelation for Capture Probabilities")))
  }
  
  ## ACF plots for category proportions(lambda)
  if(ncol(Result$lambdaDM) > 0){
     lambda.acf.plots<- do.call(arrangeGrob, c(acf.plots[(n.beta+2*ncats+1):(n.beta+3*ncats)],
                                               ncol=1,
                  main=paste("Autocorrealtion for Category Proportions(lambda)")))
  }
  
  ## ACF plots for sub-sample proportions(theta)
  if(ncol(Result$thetaDM) > 0){
    theta.acf.plots<- do.call(arrangeGrob, c(acf.plots[(n.beta+3*ncats+1):(n.beta+3*ncats+2)],
                                              ncol=1,
                      main=paste("Autocorrelation for Sub-Sample Proportions(theta)")))
  }
  
  ## ACF plots for Population size(N) and category total
  pop.acf.plots<- do.call(arrangeGrob, c(acf.plots[(n.beta+3*ncats+3):(n.beta+3*ncats+3+ncats)],
                                           ncol=1,
        main=paste("Autocorrealtion for Population size(N) and Category Total")))
  
  
  ########################## End of  ACF plots ################################
  #############################################################################
  
  #############################################################################
  ################ Gelman Rubin Diagnose  #####################################
                
  # converting results from Markov chain simulations, that might not be from BUGS,
  # to bugs object. Used mainly to display results with plot.bugs.
  mcmcm_object <- as.bugs.array(Result$thinned.array)
  
  plot(mcmcm_object)
  
  posterior.summary <- mcmcm_object$summary
  
  # create a list of mcmc objects for all the chains ( for beta parameters)
  mcmclist <- list()
  for(i in 1:Result$n.chains){
    mcmclist[[i]] <- as.mcmc(Result$thinned.array[ , i,1:n.beta])
  }
  
  
  #Potential scale reduction factors(Rhat)
  GRpsrf <- gelman.diag(mcmclist, confidence = 0.95, transform=FALSE, 
                        autoburnin=FALSE, multivariate=TRUE)
  psrf.point.est <- GRpsrf$psrf[,1] # point estimates for Rhat
  
  # Gelman and Rubin's shrink factor as the number of iterations increases.
  GR.shrink<-gelman.plot(mcmclist)
  
  # Gelman Rubin Diagnosis Plots with Point Estimates for Potential 
  # Scale Reduction Factors (Rhat)
  psrf.plots <-llply(1:n.beta,
                     function(k){qplot(x=GR.shrink$last.iter,
                                       y=GR.shrink$shrink[,k ,1], 
                                       geom="line",main=param_names[k])+
                                   geom_line(y=GR.shrink$shrink[,k ,2],
                                             linetype=2, colour="red") +
                                   ylim(min(GR.shrink$shrink[,k ,1],
                                            GR.shrink$shrink[,k ,2]), 
                                        max(GR.shrink$shrink[,k ,1],
                                            GR.shrink$shrink[,k ,2]))+
                                   xlab("last iteration in chain") +
                                   ylab("shrink factor")+
                                   annotate("text", 
                                            label = paste("\nRhat = ",
                                                          round(psrf.point.est[k],2),
                                                          "\nsolid-line : median ",
                                                          "\ndashed-lines : 97.5%",
                                                          sep=""),
                                            x =0.45*(max(GR.shrink$last.iter)-
                                                        min(GR.shrink$last.iter)), 
                                            hjust = 0,y = Inf,vjust = 1,size=3)})
    
  
  ## Gelman Rubin Diagnose Plots for beta parameters
  BGR.diag.plots<- do.call(arrangeGrob, c(psrf.plots, ncol=2,
              main=paste("Gelman Rubin Diagnosis Plots with Point Estimates", 
                         "\nfor Potential Scale Reduction Factors (Rhat)")))
  
  ############### End of  Gelman Rubin Diagnosis ##############################
  #############################################################################
  res <- NULL
  
  res$beta.trace.plots <- beta.trace.plots
  res$capture.trace.plots <- capture.trace.plots
  res$lambda.trace.plots <- lambda.trace.plots
  res$theta.trace.plots <- theta.trace.plots
  res$pop.trace.plots <- pop.trace.plots
  
  res$beta.posterior.plots <- beta.posterior.plots
  res$capture.posterior.plots <- capture.posterior.plots
  res$lambda.posterior.plots <- lambda.posterior.plots
  res$theta.posterior.plots <- theta.posterior.plots
  res$pop.posterior.plots <- pop.posterior.plots
  
  res$beta.acf.plots <- beta.acf.plots
  if(ncol(Result$captureDM) > 0){res$capture.acf.plots <- capture.acf.plots}
  if(ncol(Result$lambdaDM) > 0){res$lambda.acf.plots <- lambda.acf.plots}
  if(ncol(Result$thetaDM) > 0){res$theta.acf.plots <- theta.acf.plots}
  res$pop.acf.plots <- pop.acf.plots

  res$BGR.diag.plots <- BGR.diag.plots
  res$posterior.summary <-posterior.summary
  
  
  return(res)

} # end of the function bayesian.plots

###############################################################################
###############################################################################

