###############################################################################
########## Generate data for Partial Stratification in  #######################
##########    k-Sample Capture-Recapture Experiments    #######################
##########             with known dead removals         ####################### 
##########                    AND                       #######################
##########            Simulation  Study                 #######################
###############################################################################

# setwd("U:\\Lasantha/Research/Rcode")
# setwd("c:\\Lasantha/SFU/Research/Rcode")

source("load.R")  # load required functions and packages
set.seed(123)


######################################################## 
#### consider 2 categories and 3 sampling occasions ####
########################################################
N.gen <- 300000 # population size
st.gen <- 3 # number of sampling occasions
theta.gen <- c(0.8,0.3, 0.4) # sub sample proportions


# p_loss =  probability that animal is dead at the sampling occasion  j wherej = 1, ... , t
# this is used to get the know removals in each sampling time
p_loss.gen <- c(0.005, 0.001, 0.004)



###################################
########## Category 1 - Male ######
###################################
cat1 <- "M" # category
M_lambda <- 0.7 # category proportion for males
p_M <- c( 0.015, 0.009, 0.005) # capture probabilities for each sampling time #

#create data for cat1
cat1_data <- PSkSCR.data.generate(st=st.gen, N=N.gen, theta=theta.gen, category=cat1, 
                                  lambda=M_lambda, p=p_M, p_loss=p_loss.gen)


###################################
########## Category 1 - Female ####
###################################
cat2 <- "F" # category
F_lambda <- 0.3 # category proportion for males
p_F <- c( 0.005, 0.01, 0.007) # capture probabilities for each sampling time  

#create data for cat2
cat2_data <- PSkSCR.data.generate(st=st.gen, N=N.gen, theta=theta.gen, category=cat2, 
                                  lambda=F_lambda, p=p_F, p_loss=p_loss.gen)

###################################
## combine created data
gen.data <- rbind( cat1_data ,cat2_data )

# data to be used to analysis 
data <- NULL
data$history <- as.character(gen.data$history)
data$counts <- gen.data$counts
data$category <- c(cat1,cat2)
str(data)

data.frame(data$history,data$counts)

## The following data set is used for the PSkSCR Analysis 
#write.csv(data, file = "PSkSCR_Test_data_1.csv")


###############################################################################
###############################################################################


###############################################################################
###############################################################################
### to check the simulated data : PSkSCR experiments         ##################
### Simulate "n.samples" number of samples  and estimate the parameters    ####
### from those simulated samples and and create histograms for the         ####
### parameter estimates of those simulated samples                         ####
###############################################################################


### parameter values given to simulate data sets ##############################

N.gen <- 300000 # population size
st.gen <- 3 # number of sampling occasions
theta.gen <- c(0.8,0.3, 0.4) # sub sample proportions

# p_loss =  probability that animal is dead at the sampling occasion  j wherej = 1, ... , t
# this is used to get the know removals in each sampling time
p_loss.gen <- c(0.005, 0.001, 0.004)

cat1 <- "M" # category
M_lambda <- 0.7 # category proportion for males
p_M <- c( 0.015, 0.009, 0.005) # capture probabilities for each sampling time #

cat2 <- "F" # category
F_lambda <- 0.3 # category proportion for males
p_F <- c( 0.005, 0.01, 0.007) # capture probabilities for each sampling time  

### parameters used  are as a vector ##############
cap.probability <- as.vector(rbind(p_M,p_F))
cat.proportion <- c(M_lambda,F_lambda)
parameters <- c(cap.probability,cat.proportion,theta.gen,p_loss.gen,N.gen )

###############################################################################

#model identification : unrestricted model
model.id <- paste("{ p(c*t), theta(t), lambda(c) }")

# give  the required design matrices for capture recapture probabilities,
# theta(sampling(sexing) fractions ) and lambda(Category proportion)
captureDM <- create.DM(c(1,2,3,4,5,6)) 
thetaDM   <- create.DM(c(1,2,3))
lambdaDM  <- create.DM(c(1)) 
p_lossDM  <- create.DM(c(1,2,3))

#give the offset vectors(vectors of zero's should be given since no restriction)
captureOFFSET <- c(0,0,0,0,0,0) 
thetaOFFSET   <- c(0,0,0)
lambdaOFFSET  <- c(0)
p_lossOFFSET  <- c(0,0,0)

n.samples=1000 # number of simulated data sets

# estimates of the parameters for n.samples  number of simulated sample
est.n.samples <- rdply(n.samples,
                       PSkSCR.est.for.one.simulated.data.set(N.gen, st.gen, theta.gen,p_loss.gen,
                                                             cat1, M_lambda, p_M,
                                                             cat2, F_lambda, p_F,
                                                             model.id,
                                                             captureDM,thetaDM,lambdaDM, p_lossDM,
                                                             captureOFFSET,thetaOFFSET,
                                                             lambdaOFFSET,p_lossOFFSET)) 
 
cats <- c(cat1,cat2 )
names(est.n.samples)<- c("Sample", paste(rep(paste("p_",rep(1:st.gen,each=length(cats)),
                                                   sep="")), rep(cats,st.gen),sep=""),
                         paste("lambda_",cats,sep=""),
                         paste("theta_",c(1:st.gen), sep=""), 
                         paste("nu_",c(1:st.gen), sep=""),"N")

## means of the estimates of the simulated samples
mean.est.simulated <- round(apply(est.n.samples, 2, mean)[2:ncol(est.n.samples)],5)

## sd of the estimates of the simulated samples
sd <- round(apply(est.n.samples, 2, sd)[2:ncol(est.n.samples)],5)

## initial parameter values used
names(parameters) <- names(est.n.samples)[2:ncol(est.n.samples)]
data.frame( parameters, mean.est.simulated,sd)


###### estimates from the data set that was used for the #####
###### k-Sample Capture-Recapture Experiments analysis   #####
# get  data that ws used for the analysis
data <- PSkSCR.get.data("PSkSCR_Test_data_1.csv")
str(data)
## fit the model and get estimates
MLE_PSkSCR_best_model <- PSkSCR.fit.model(model.id,data,
                                          captureDM,thetaDM,lambdaDM,p_lossDM,
                                          captureOFFSET,thetaOFFSET,
                                          lambdaOFFSET,p_lossOFFSET)
est.best.model <- MLE_PSkSCR_best_model$est$full


# plot the histograms of the estimates of the simulated data sets
# Solid lines represent the parameter values used for data simulations
# dashed lines represent the estimates foud using the data set which was
# used for PSkSCR Analysis
plots.hist <- llply(2:ncol(est.n.samples),
               function(k){qplot(x=est.n.samples[,k],
                                 geom = "histogram", 
                                 main=paste(names(est.n.samples)[k],sep="") )+
                   theme(plot.title = element_text(size = rel(2)))+
                   xlab(names(est.n.samples)[k]) + ylab("Frequancy")+
                   geom_vline(xintercept=parameters[k-1],colour="red", size=1.1)+
                   geom_vline(xintercept=est.best.model[k-1],colour="red",
                              linetype = "longdash",size=1)})
                                           

## all the hitograms of the estimates  parameters in a one plot
all.plots <- do.call(arrangeGrob, c(plots.hist, ncol=3,
                      main=paste("Histograms of the Estimates of the Parameters")))
                                   
                                                     
###############################################################################
############################################################################### 






############################################################################### delete later
# library(grid)
# grid.arrange(rectGrob(), rectGrob())
# 
# ## Not run:
# library(ggplot2)
# pl <- lapply(1:11, function(.x) qplot(1:10,rnorm(10), main=paste("plot",.x)))
# ml <- marrangeGrob(pl, nrow=2, ncol=2)
# ## interactive use; open new devices
# ml
# ## non-interactive use, multipage pdf
# ggsave("multipage.pdf", ml)
# ## End(Not run)


# 
# grobframe <- arrangeGrob(p1, p2, ncol=1, nrow=2,
#                          main = textGrob("\nArrangeGrob Test", gp = gpar(fontsize=18, fontface="bold.italic", fontsize=18)),
#                          sub = textGrob("*subtitle location*", x=0, hjust=-0.5, vjust=0.1, gp = gpar(fontface = "italic", fontsize = 15)))
# 
# 
# aaaaaaaaa <- do.call(arrangeGrob, c(plots.hist, ncol=3,
#                                     main = textGrob("Histograms of the Estimates of the Parameters", 
#                                                     gp = gpar(fontsize=18, fontface="bold", fontsize=18))))
# 
# 
# marrangeGrob(plots.hist, nrow=4, ncol=3, main = textGrob("Histograms of the Estimates of the Parameters", 
#                                                          gp = gpar(fontsize=18, fontface="bold", fontsize=18)))
