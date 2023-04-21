###############################################################################
###############################################################################
####################  calculation of effort  ##################################
##### ( non linear optimization with a linear constraints for total cost)  ####
#####              and creating Contour plots                              ####
###############################################################################

#setwd("D:\\SFU/Research/Rcode/Rcode")
#setwd("U:\\Lasantha/Research/Rcode")
#setwd("c:\\Lasantha/SFU/Research/Rcode")


source("load.R")  # load required functions

set.seed(123)
Data <- get.data("Mille_Lacs_Walleye.csv")
str(Data)

###############################################################################
############ Optimization under linear equality Constraint  ###################
###                                                                        ####
###  c0 + n1*c1 +  n1.star*c1.star + n2*c2 + n2.star*c2.star = total.cost  ####
###                                                                        ####
###############################################################################

# Optimal allocation of the sample size and the sub-sample size at both
# sampling times are calculated for given guestimates and for given cost
# by minimizingf the standard error of the population size N.

# The output is the optimized values for the counts of the capture histories
# U0, UU, 0U, C0, CC, 0C
#    and
# n1 - sample size at time 1
# n1.star - subsample size at time 1
# n2 - sample size at time 2
# n2.star - sub-sample size at time 2
#    and
# standard error of population size N at the optimized values
#
# n1, n1.star, n2, and n2.star are found after optimising the counts of
# the capture histories U0, UU, 0U, C0, CC, 0C. Have to do this way because
# sample size at time 2 depends on the capture histories CC and UU that is also
# related to sample size at time 1. and also there is a restriction such that
# n2.star <= n2 - E(n_UU + n_CC)


###############################################################################
##  guesstimate  for population size (N), category proportions (lambda)
##  and ratio of capture probabilities of  male to female at each sampling times

         N <- 205500  # guesstimate  for N
cat_prop <- c(0.33,0.67) # c(0.3298384, 0.6701616)   # guesstimate  for lambda
      r1 <- 6.5  # 0.07580455/0.01157313     # guesstimate  of ratio of p1m/p1f
      r2 <- 0.4 #  0.007892997/0.020999333   # guesstimate  of ratio of p2m/p2f

category <- c("M", "F") # consider two categories in the populations


################ Costs related to capture and stratification ##################
## these are the actual cost for the Mille Lacs Walleye data
#
#      c0 <- 0           # fixed cost for the study
#      c1 <- (10*60)/38  # cost for  sampling at time 1
# c1.star <- 1/12        # cost for  sub-sampling and processing  at time 1
#      c2 <- (10*60)/10  # cost for  sampling at time 2
# c2.star <- 2           # cost for  sub-sampling and processing  at time 2
# total.cost <- 313028   # total available cost
#
# cost.vec <- c( c0, c1, c1.star, c2, c2.star)
#####################
c1 <- 4 #2   # cost for  sampling at time 1
c2 <- 6 # 3   # cost for  sampling at time 2
c1.star <-0.4 # .2   # cost for  sub-sampling and processing  at time 1
c2.star <-0.4 # 0.2  # cost for  sub-sampling and processing  at time 2
c0 <- 0            # fixed cost for the study
total.cost <- 90000  #45000      # total available cost

cost.vec <- c( c0, c1, c1.star, c2, c2.star)
names(cost.vec) <- c("c0", "c1","c1.star", "c2", "c2.star")

####### lower bound and the upper bounds for the optimization routine #########
lower.b <- rep(0,6)
upper.b <- rep(15000,6)


##############################################################################
#### Initial values for the capture histories U0, UU, C0, CC, 0U, and 0C #####

## starting parameter vector should be  in the following order
##  ( n_U0, n_UU, n_0U, n_C0, n_CC, n_0C )
## give the values for starting parameter vector or get it from
## the Mille Lacs Walleye dataset

#x0 <- c(50,100,2000,5000,200,500)
x0 <- c(10,10,3000,4000,30,1000)
names(x0) <- c("n_U0", "n_UU", "n_0U", "n_C0", "n_CC", "n_0C")

#x0 <- starting.param.optimal.allocation(Data) #from the Mille Lacs Walleye dataset
#x0<- c(10,20,1000,6000,30,3000)
##############################################################################

#model identification"
model.id <- paste("{ p(c*t), theta(t), lambda(.) }")

#get the required design matrices
captureDM <- create.DM(c(1,2,3,4)) # Design matrix for capture recapture probabilities
thetaDM   <- create.DM(c(1,2))     # Design matrix for theta(sampling(sexing) fractions)
lambdaDM  <- matrix(0, ncol=0,nrow=1)       # Design matrix for lambda(Category proportion)

#give the offset vectors(vectors of zeros should be given since no restriction)
captureOFFSET <- c(0,0,0,0)
thetaOFFSET   <- c(0,0)
lambdaOFFSET  <- (logit(cat_prop[1]))


##Optimal allocation under linear equality Constraint
OA.eqal.constraint <- planning.tool.optimal.allocation( N, cat_prop,r1,r2,
                                                        category,
                                                        cost.vec,total.cost,x0,
                                                        lower.b, upper.b,
                                                        model.id,
                                                        captureDM,thetaDM, lambdaDM,
                                                        captureOFFSET, thetaOFFSET,
                                                        lambdaOFFSET)

# print output
print.optimal.allocation.equality.constraint(OA.eqal.constraint)

###############################################################################
###############################################################################



###############################################################################
############ Optimal allocation under linear inequality Constraint ############
##                                                                            #
## LB <= c0 + n1*c1 +  n1.star*c1.star + n2*c2 + n2.star*c2.star <= UB        #
##                                                                            #
### where LB is lower bound  and UB LB is upper boundfor the available total ##
### cost, so that linear cost function is between LB and UB (total.cost)     ##
###############################################################################
                                                                              #
inequalityLB <- 80000                                                         #
inequalityUB <- total.cost                                                    #
                                                                              #
OA.inq.cont <- planning.tool.optimal.allocation.inequality.cost(N,cat_prop,   #
                                                                r1,r2,        #
                                                                category,     #
                                                                cost.vec,     #
                                                                inequalityLB, #
                                                                inequalityUB, #
                                                                x0,           #
                                                                lower.b,      #
                                                                upper.b,      #
                                                                model.id,     #
                                                                captureDM,    #
                                                                thetaDM,      #
                                                                lambdaDM,     #
                                                                captureOFFSET,#
                                                                thetaOFFSET,  #
                                                                lambdaOFFSET) #
                                                                              #
# print output                                                                #
print.optimal.alloc.inequality.constraint(OA.inq.cont)                        #
                                                                              #
###############################################################################
###############################################################################

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
####################### End of Optimal Allocation  ############################
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##




######################## CONTOUR PLOT #########################################
###                                                                         ###
###  Contour plot for standard error of N for given n1 and n2 whre n1.star  ###
###  are and n2.star fixed at the optimised values                          ###
###                                                                         ###
###############################################################################


# Optimized values for n1, n1.star, n2, and n2.star
     n1 <- 9012
n1.star <- 8991
     n2 <- 8297
n2.star <- 1451


## create a grid of values for n1 and n2
## n1.star and n2.star are fixed at the optimal value
n1grid <- seq(n1,n1+3000,length.out=5)
n2grid <- seq(n2-6000,n2+3000,length.out=5)
grid.n1.n2 <- expand.grid(n1grid  ,n1.star,  n2grid,  n2.star)


## create a grid of values for n1.star and n2.star
## n1 and n2  are fixed at the optimal value
n1.stargrid <- seq(n1.star-7000,n1.star+20,length.out=30)
n2.stargrid <- seq(n2.star-1450,n2.star+6000,length.out=30)

grid.n1.star.n2.star <- expand.grid(n1  ,n1.stargrid,  n2,  n2.stargrid)



######## n1.star and n2.star are fixed at optimized values ###########
N.SE.grid.n1.n2  <- adply( grid.n1.n2, 1,function(x){
                      contour.plot.SE.N(unlist(x),N=N, cat_prop=cat_prop,r1=r1,r2=r2,
                                model.id=model.id,
                                captureDM=captureDM,thetaDM=thetaDM, lambdaDM=lambdaDM,
                                captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
                                lambdaOFFSET=lambdaOFFSET)})
names(N.SE.grid.n1.n2)=c("n1","n1.star","n2","n2.star","SE_N")
nrow(N.SE.grid.n1.n2)
head(N.SE.grid.n1.n2)
summary(N.SE.grid.n1.n2$SE_N)


contour.plot.n1.n2 <- ggplot(N.SE.grid.n1.n2, aes(n1, n2, z = SE_N))+
                      stat_contour(aes(colour = ..level..),size = 1) +
                      ggtitle(paste("Contour plot of SE of N when n1.star and n2.star are fixed" )) +
                      theme(plot.title = element_text(family = "Trebuchet MS",face="bold", size=18)) +
                      theme(axis.title = element_text(family = "Trebuchet MS",face="bold", size=16))

direct.label(contour.plot.n1.n2) # this produce contours with values



#total cost calcualtion for each row of values in the N.SE.grid.n1.n2
N.SE.grid.n1.n2$TotalCost <- as.matrix(cbind(rep(1,nrow(N.SE.grid.n1.n2)),N.SE.grid.n1.n2[,c("n1","n1.star","n2","n2.star")])) %*% cost.vec
#keep the rows in N.SE.grid.n1.n2 if the cost is less than  or equal "total.cost"
reduced.N.SE.grid.n1.n2 <- N.SE.grid.n1.n2[N.SE.grid.n1.n2$TotalCost <= total.cost,]
nrow(reduced.N.SE.grid.n1.n2)
summary(reduced.N.SE.grid.n1.n2$SE_N)

#contour plot for reduced df (i.e.reduced.N.SE.grid.n1.n2 )
reduced.contour.plot.n1.n2 <- ggplot(reduced.N.SE.grid.n1.n2, aes(n1, n2, z = SE_N))+
                              stat_contour(aes(colour = ..level.. ),size = 1) +
                              ggtitle(paste("Contour plot of SE of N when n1.star and n2.star are fixed" )) +
                              theme(plot.title = element_text(family = "Trebuchet MS",face="bold", size=18)) +
                              theme(axis.title = element_text(family = "Trebuchet MS",face="bold", size=16))

direct.label(reduced.contour.plot.n1.n2) # this produce contours with values



######## n1 and n2 are fixed at optimized values ###########
N.SE.grid.n1.star.n2.star  <- adply( grid.n1.star.n2.star , 1,function(x){c
                               contour.plot.SE.N(unlist(x),N=N, cat_prop=cat_prop,r1=r1,r2=r2,
                                 model.id=model.id,
                                 captureDM=captureDM,thetaDM=thetaDM, lambdaDM=lambdaDM,
                                 captureOFFSET=captureOFFSET, thetaOFFSET=thetaOFFSET,
                                 lambdaOFFSET=lambdaOFFSET)})
names(N.SE.grid.n1.star.n2.star)=c("n1","n1.star","n2","n2.star","SE_N")
nrow(N.SE.grid.n1.star.n2.star)
head(N.SE.grid.n1.star.n2.star)
summary(N.SE.grid.n1.star.n2.star$SE_N)

contour.plot.n1.star.n2.star <- ggplot(N.SE.grid.n1.star.n2.star, aes(n1.star, n2.star, z = SE_N))+
                                stat_contour(aes(colour = ..level..),size = 1) +
                                ggtitle(paste("Contour plot of SE of N when n1 and n2 are fixed" )) +
                                theme(plot.title = element_text(family = "Trebuchet MS",face="bold", size=18)) +
                                theme(axis.title = element_text(family = "Trebuchet MS",face="bold", size=16))

direct.label(contour.plot.n1.star.n2.star)



#total cost calcualtion for each row of values in the N.SE.grid.n1.star.n2.star
N.SE.grid.n1.star.n2.star$TotalCost = as.matrix(cbind(rep(1,nrow(N.SE.grid.n1.star.n2.star)),
                                                      N.SE.grid.n1.star.n2.star[,c("n1","n1.star",
                                                                                   "n2","n2.star")])) %*% cost.vec
#keep the rows in N.SE.grid.n1.n2 if the cost is less than  or equal "total.cost"
Reduced.N.SE.grid.n1.star.n2.star <- N.SE.grid.n1.star.n2.star[N.SE.grid.n1.star.n2.star$TotalCost <= total.cost,]
nrow(Reduced.N.SE.grid.n1.star.n2.star)
summary(Reduced.N.SE.grid.n1.star.n2.star$SE_N)

reduced.contour.plot.n1.star.n2.star <- ggplot(Reduced.N.SE.grid.n1.star.n2.star, aes(n1.star, n2.star, z = SE_N))+
                  stat_contour(aes(colour = ..level..),size = 1) +
                  ggtitle(paste("Contour plot of SE of N when n1 and n2 are fixed" ))+
                  theme(plot.title = element_text(family = "Trebuchet MS",face="bold", size=18)) +
                  theme(axis.title = element_text(family = "Trebuchet MS",face="bold", size=16))

direct.label(reduced.contour.plot.n1.star.n2.star)



##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#######################  End of Contour Plots  ################################
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### delete the following later...
# ggplot(Reduced.N.SE.grid.n1.star.n2.star, aes(n1.star, n2.star, z = SE_N))+
#   geom_raster(aes(fill = SE_N))  +
#   stat_contour(colour = "white",size = 1)+
#   ggtitle(paste("Contour plot of SE of N when n1 and n2 are fixed" ))+
#   theme(plot.title = element_text(family = "Trebuchet MS",face="bold", size=18)) +
#   theme(axis.title = element_text(family = "Trebuchet MS",face="bold", size=16))

