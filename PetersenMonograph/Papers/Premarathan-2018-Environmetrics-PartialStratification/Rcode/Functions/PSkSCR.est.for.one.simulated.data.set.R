############################################################################### 
###########to check the simulated data : PSkSCR experiments ###################
### Simulate a one data set and estimate the parameters from that data set ####
###############################################################################

### this function return the estimate of the parameters for one simulated data set  ##

# Function : PSkSCR.est.for.one.simulated.data.set
#  Input : parameter values to simulate a data set and
#          model id, design matrices and offset vectors
#  Output : estimates of for the parameters for the simulated data set

PSkSCR.est.for.one.simulated.data.set <- function(N.gen, st.gen, theta.gen,p_loss.gen,
                                                  cat1, M_lambda, p_M,
                                                  cat2, F_lambda, p_F,
                                                  model.id,
                                                  captureDM,thetaDM,lambdaDM,p_lossDM,
                                                  captureOFFSET,thetaOFFSET,
                                                  lambdaOFFSET,p_lossOFFSET){
  #create data for cat1
  cat1_data <- PSkSCR.data.generate(st=st.gen, N=N.gen, theta=theta.gen, category=cat1, 
                                    lambda=M_lambda, p=p_M, p_loss=p_loss.gen)
  
  #create data for cat2
  cat2_data <- PSkSCR.data.generate(st=st.gen, N=N.gen, theta=theta.gen, category=cat2, 
                                    lambda=F_lambda, p=p_F, p_loss=p_loss.gen)
  
  ## combine created data
  gen.data <- rbind( cat1_data ,cat2_data )
  
  # data to be used to analysis 
  data <- NULL
  data$history <- as.character(gen.data$history)
  data$counts <- gen.data$counts
  data$category <- c(cat1,cat2)
  str(data)
  
  MLE_PSkSCR_model_full_estimates <- PSkSCR.fit.model(model.id=model.id,
                                                      data=data,
                                                      captureDM=captureDM,
                                                      thetaDM=thetaDM,
                                                      lambdaDM=lambdaDM,
                                                      p_lossDM=p_lossDM,
                                                      captureOFFSET=captureOFFSET,
                                                      thetaOFFSET=thetaOFFSET,
                                                      lambdaOFFSET=lambdaOFFSET,
                                                      p_lossOFFSET=p_lossOFFSET)$est$full
  return(MLE_PSkSCR_model_full_estimates)
  
} # end of the function " PSkSCR.est.for.one.simulated.data.set" 

###############################################################################
###############################################################################