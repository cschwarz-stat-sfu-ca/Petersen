###############################################################################
#########    Create design matrix for the  capture probabilities       ######## 
#########                for data with individual covariates           ########
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
############################################################################### 

# Function :ic.create.DM 

# First create a data frame to create a design matrix for capture probabilities
# dimension of the  data frame : 
#  number of rows = (n.rows x n.cat x 2)
#                           # 2 for time 1 and time 2
#  number of columns = (number of individual covariates + 2)
#                       # +2 represent a column for "category" and column for "time"
#
# if there are two categories(male, female) then there are 4 blocks such as
#  male at time 1,female at time 1, male at time 2, female at time 2 respectively
#  so that total number of rows in the big data frame is  length(Data$History)x4
#
#
#Function :ic.create.DM
# Input : Data and formula
#         ( formula should contain the exact names of the individual covariates and
#           if the "category" and "time" is in the formula that should taken as
#           reserved words)
#
# Output : Design matrix 
#  number of rows = (n.rows x n.cat x 2)
#                           # 2 for time 1 and time 2
#  number of columns = depends on the formula; first column is for the intercept
#

ic.create.DM <- function(Data,captureformula){
  
  n.rows <- length(Data$History) # number of rows of data
  n.cat <- length(Data$category) # number of categories
  
  ind <- c(which(names(Data)=="History"), which(names(Data)=="counts"),
           which(names(Data)=="category") )
  
  # get only the individual covariates  from Data. These are used to create 
  # the data frame and then create the design matrix according to the formula
  ind_cov_data <- Data[-ind] 
  
  new.data.matrix <- matrix(c( c(do.call(cbind,rep(ind_cov_data,each = (n.cat*2) ))),
                               rep(0:(n.cat-1),each=n.rows,len=n.rows*n.cat*2), 
                               rep( 0:1, each=(n.rows*n.cat)) ),
                            ncol = (length(names(ind_cov_data))+2),
                            nrow = (n.rows*n.cat * 2),
                            byrow=FALSE)
  new.data.frame <- data.frame(new.data.matrix)
  colnames(new.data.frame) <- c(names(ind_cov_data), "category","time")
  
  captureDM <- model.matrix(captureformula,new.data.frame)
  
  return(captureDM)
  
} # end of the function ic.create.DM



###############################################################################
################## Test the above functions ###################################
###############################################################################
# 
# ## test data
# 
# testdata<- NULL
# testdata$History <- c("MM", "MM", "FF", "FF", "UU", "U0","0U")
# testdata$counts <- c(1,1,1,1,1,1,1)
# testdata$length <- c(1,2,3,4,5,6,7)
# testdata$lengthsq <- testdata$length^2
# testdata$weight <- c(22,23,24,25,26,27,28)
# testdata$category <- c("A", "B")
# str(testdata)
# 
# ##########  Test the function "ic.create.DM" ############################
# captureformula1 <- ~length+weight+category +time
# ic.create.DM(testdata,captureformula1)
# 
# captureformula2 <- ~length+category +time
# ic.create.DM(testdata,captureformula2)
# 
# captureformula3 <- ~weight+category +time
# ic.create.DM(testdata,captureformula3)
# 
# captureformula4 <- ~length+weight
# ic.create.DM(testdata,captureformula4)
# 
# captureformula5 <- ~length+weight+category*time
# ic.create.DM(testdata,captureformula5)
# 
# captureformula6 <- ~length+I(length^2) + weight+category+time
# ic.create.DM(testdata,captureformula6)
# 
# 
# captureformula7 <- ~length+I(length^2) + weight + I(weight^2)+category*time
# ic.create.DM(testdata,captureformula7)
# 
# captureformula8 <- ~length+ lengthsq+category+time
# ic.create.DM(testdata,captureformula8)
# 
# captureformula9 <- ~length+ lengthsq+category*time
# ic.create.DM(testdata,captureformula9)
# 
# captureformula10 <- ~length+ lengthsq+weight +category*time
# ic.create.DM(testdata,captureformula10)
#
# captureformula11 <- ~length*category*time+ lengthsq*category*time+ category*time
# ic.create.DM(testdata,captureformula11)
###############################################################################
################## End of test the function ic.create.DM ######################
###############################################################################


