############################################################################### 
## predictive Plots of Capture Probabilities (logit scale and regular scale) ##
#################### for the data with individual covariates ##################
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
###############################################################################

# Function "ic.fit.model" 
#  Inputs: optimized values for parameters(logit form), Data( here the length and 
#           lengthsq are standardized separately), true_length, capture formula, 
#           design matrices and offset vectors for category proportions and
#           sub-sample proportions
#  Outputs: Plot of Capture Probabilities (logit scale and regular scale)

ic.plots <- function (logit.est ,Data,true_length,
                       captureformula, thetaDM,lambdaDM,
                       thetaOFFSET,lambdaOFFSET){
  
  true.length <- true_length # true length
  ncats <- length(Data$category) # number of categories
  
  # create values of length to use for predictions
  min.length <- min(Data$length)
  max.length <- max(Data$length)
  
  # create sequence of points to plot the capture probabilities
  seq_length   <- seq(from=min.length,to=max.length,length=100)
  seq_lengthsq <- seq_length^2
  
  # these sequence of points are in standardized scale. Now convert them into 
  # the regular scale. 
  seq_length_regular <- seq_length + mean(true.length)
  

  # data for prediction with standardized lengths 
  pred_data <- Data
  pred_data$length <- seq_length
  pred_data$lengthsq <- seq_lengthsq 
  pred_data$History  <- rep("pred",length(seq_length) )
  
  
  # Design matrix for created values
  pred.captureDM <- ic.create.DM(pred_data,captureformula)
  
  new.captureOFFSET = matrix(0, nrow=(nrow(pred.captureDM)/(length(Data$category)*2)),
                             ncol=(length(Data$category)*2)) 
  
  cap.prob <- as.vector(ic.unpack.parm(logit.est,  pred_data ,
                                       captureDM=pred.captureDM,thetaDM,lambdaDM,
                                       captureOFFSET=new.captureOFFSET,thetaOFFSET,
                                       lambdaOFFSET)$p)

  
  ########## Plot of Capture Probabilities (logit scale) ###################### 
  ncats <- length(Data$category)
  pred.beta.data.matrix <- matrix(c(round(logit(cap.prob),4), 
                                    rep(round(seq_length,2),ncats*2)),
                                  ncol = 2, 
                                  nrow = (length(seq_length)*ncats * 2),
                                  byrow=FALSE)
  
  pred.logit.data.frame <- data.frame(pred.beta.data.matrix)
  
  pred.logit.data.frame <- cbind(pred.logit.data.frame ,
                                 c(rep(c("Male_t1", "Female_t1"),
                                       each=length(seq_length)),
                                   rep(c("Male_t2", "Female_t2"),
                                       each=length(seq_length))))
  
  colnames(pred.logit.data.frame) <- c("probability","length", "category_time")
  
  pred.logit.ggplot <- ggplot(pred.logit.data.frame,
                              aes(x=length,y=probability,colour=category_time))+
                          geom_line(size=1, aes(linetype =category_time)) +
                          ggtitle(paste("Plot of Capture Probabilities (logit scale) with individual covariate 'length' \n Capture Formula : ",
                                  paste(captureformula)[1], paste(captureformula)[2],sep=" " )) +
                          xlab("Length") + 
                          ylab("Probability ( in logit scale)")+
                          theme(axis.text=element_text(size=12,face="bold"),
                                axis.title=element_text(size=12),
                                plot.title = element_text(size = 12,face="bold"),
                                legend.title = element_text(size = 12,face="bold"),
                                legend.text = element_text(size = 12) ) 
  
    
  
  ######### Plot of Capture Probabilities (regular scale) #####################
  pred.data.matrix <- matrix(c(round(cap.prob,4),
                               rep(round(seq_length,2),ncats*2)),
                             ncol = 2, 
                             nrow = (length(seq_length)*ncats * 2),
                             byrow=FALSE)
  
  pred.data.frame <- data.frame(pred.data.matrix)
  
  pred.data.frame <- cbind(pred.data.frame ,
                           c(rep(c("Male_t1", "Female_t1"),
                                 each=length(seq_length)),
                             rep(c("Male_t2", "Female_t2"),
                                 each=length(seq_length))))
  
  colnames(pred.data.frame) <- c("probability","length", "category_time")
  
  pred.ggplot <- ggplot(pred.data.frame, 
                        aes(x=length, y=probability,colour=category_time)) +
                    geom_line(size=1, aes(linetype =category_time)) +
                    ggtitle(paste("Plot of Estimated Capture Probabilities with individual covariate 'length' \n Capture Formula : ",
                            paste(captureformula)[1], paste(captureformula)[2],sep=" " )) +
                    xlab("Standardized Length") + 
                    ylab("Probability")+
                    theme(axis.text=element_text(size=12,face="bold"),
                          axis.title=element_text(size=12),
                          plot.title = element_text(size = 12,face="bold"),
                          legend.title = element_text(size = 12,face="bold"),
                          legend.text = element_text(size = 12) ) 
  
  ##############################
  
  # Plot of Capture Probabilities (regular scale)  with length in regular scale
  pred.data.matrix_reg <- matrix(c(round(cap.prob,4),
                                 rep(round(seq_length_regular,2),ncats*2)),
                                 ncol = 2, 
                                 nrow = (length(seq_length_regular)*ncats * 2),
                                 byrow=FALSE)
  
  pred.data.frame_reg <- data.frame(pred.data.matrix_reg)
  
  pred.data.frame_reg <- cbind(pred.data.frame_reg ,
                               c(rep(c("Male_t1", "Female_t1"),
                                     each=length(seq_length)),
                                 rep(c("Male_t2", "Female_t2"),
                                      each=length(seq_length)) ) )
  
  colnames(pred.data.frame_reg) <- c("probability","length", "category_time")
  
  pred.ggplot_reg <- ggplot(pred.data.frame_reg, 
                            aes(x=length, y=probability,colour=category_time)) +
                      geom_line(size=1, aes(linetype =category_time)) +
                      ggtitle(paste("Plot of Estimated Capture Probabilities with individual covariate 'length' \n Capture Formula : ",
                              paste(captureformula)[1], paste(captureformula)[2],sep=" " )) +
                      xlab("Length") + 
                      ylab("Probability")+
                      theme(axis.text=element_text(size=12,face="bold"),
                      axis.title=element_text(size=12),
                      plot.title = element_text(size = 12,face="bold"),
                      legend.title = element_text(size = 12,face="bold"),
                      legend.text = element_text(size = 12) ) 
  
  ###########
  
  out<- NULL
  out$pred.logit.ggplot <- pred.logit.ggplot
  out$pred.ggplot <- pred.ggplot
  out$pred.ggplot_reg <- pred.ggplot_reg
  
  return(out)
  
} # end of ic.plots



###############################################################################
###############################################################################
## plots for the generated data without standardizing the length and lengthsq #
###############################################################################

## this function only use to see the plots for capture probabilities with 
## individual covariates with generated lengths. Here lengths are not standardise.
## However the shape of the graphs from the analysis( after standardizing the
## generated length and lengthsquare separately) should be closer to the graphs
## creating from the following code even though the x-axis scales are different

## ***here,  Data$length and Data$lengthsq are not standardized ******** ###
  
ic.generated.data.plots <- function (logit.est ,Data,
                        captureformula, thetaDM,lambdaDM,
                        thetaOFFSET,lambdaOFFSET){
    
    ncats <- length(Data$category) # number of categories
    
    # create values of length to use for predictions
    min.length <- min(Data$length)
    max.length <- max(Data$length)
    length.values <- seq(from=min.length,to=max.length,length=50)
    
    pred_data <- Data
    pred_data$length <- length.values
    pred_data$lengthsq <- length.values^2
    pred_data$History  <- rep("pred",length(length.values) )
    
    # Design matrix for crated values
    pred.captureDM <- ic.create.DM(pred_data,captureformula)
    
    new.captureOFFSET = matrix(0, nrow=(nrow(pred.captureDM)/(length(Data$category)*2)),
                               ncol=(length(Data$category)*2)) 
    
    cap.prob <- as.vector(ic.unpack.parm(logit.est, pred_data,
                                         captureDM=pred.captureDM,thetaDM,lambdaDM,
                                         captureOFFSET=new.captureOFFSET,thetaOFFSET,
                                         lambdaOFFSET)$p)
    
    
    ########## Plot of Capture Probabilities (logit scale) ###################### 
    ncats <- length(Data$category)
    pred.beta.data.matrix <- matrix(c(round(logit(cap.prob),4), 
                                      rep(round(length.values,2),ncats*2)),
                                    ncol = 2, 
                                    nrow = (length(length.values)*ncats * 2),
                                    byrow=FALSE)
    
    pred.logit.data.frame <- data.frame(pred.beta.data.matrix)
    
    pred.logit.data.frame <- cbind(pred.logit.data.frame ,
                                   c(rep(c("Male_t1", "Female_t1"),
                                         each=length(length.values)),
                                     rep(c("Male_t2", "Female_t2"),
                                         each=length(length.values))))
    
    colnames(pred.logit.data.frame) <- c("probability","length", "category_time")
    
    pred.logit.ggplot <- ggplot(pred.logit.data.frame,
                                aes(x=length,y=probability,colour=category_time))+
                            geom_line(size=1, aes(linetype =category_time)) +
                            ggtitle(paste("Plot of Capture Probabilities (logit scale) with Individual Covariate 'length' \n For the Generated lengths \n Capture Formula : ",
                                          paste(captureformula)[1], paste(captureformula)[2],sep=" " )) +
                            xlab("Length") + 
                            ylab("Probability ( in logit scale)")+
                            theme(axis.text=element_text(size=12,face="bold"),
                                 axis.title=element_text(size=12),
                                 plot.title = element_text(size = 12,face="bold"),
                                 legend.title = element_text(size = 12,face="bold"),
                                 legend.text = element_text(size = 12) ) 
    
    
    
    ######### Plot of Capture Probabilities (regular scale) #####################
    pred.data.matrix <- matrix(c(round(cap.prob,4),
                                 rep(round(length.values,2),ncats*2)),
                               ncol = 2, 
                               nrow = (length(length.values)*ncats * 2),
                               byrow=FALSE)
    
    pred.data.frame <- data.frame(pred.data.matrix)
    
    pred.data.frame <- cbind(pred.data.frame ,
                             c(rep(c("Male_t1", "Female_t1"),
                                   each=length(length.values)),
                               rep(c("Male_t2", "Female_t2"),
                                   each=length(length.values))))
    
    colnames(pred.data.frame) <- c("probability","length", "category_time")
    
    pred.ggplot <- ggplot(pred.data.frame, 
                          aes(x=length, y=probability,colour=category_time)) +
      geom_line(size=1, aes(linetype =category_time)) +
      ggtitle(paste("Plot of Estimated Capture Probabilities with Individual Covariate 'length' \n For the Generated lengths \n Capture Formula : ",
                    paste(captureformula)[1], paste(captureformula)[2],sep=" " )) +
      xlab("Length") + 
      ylab("Probability")+
      theme(axis.text=element_text(size=12,face="bold"),
            axis.title=element_text(size=12),
            plot.title = element_text(size = 12,face="bold"),
            legend.title = element_text(size = 12,face="bold"),
            legend.text = element_text(size = 12) ) 
    
    #################
    
    out<- NULL
    out$pred.logit.ggplot <- pred.logit.ggplot
    out$pred.ggplot <- pred.ggplot
    
    return(out)
    
  } # end of ic.generated.data.plots

###############################################################################