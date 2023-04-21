###############################################################################
########### Read data set with individual covariates ########################## 
####  Partial Stratification in two-Sample Capture-Recapture Experiments   ####
###############################################################################

# Function "ic.get.data"
#   Input: Data Set in Raw CSV form. 
#   Output: Data Set(as a list) ready for analysis. 

ic.get.data = function(testdata.csv){ 
  
  datacsv = read.csv(testdata.csv, header = TRUE, as.is=TRUE, strip.white=TRUE)
  
  data = NULL
  data$History   <- datacsv$history  # capture histories of the raw data
  data$counts    <- datacsv$counts  # counts of the corresponding capture histories
  data$length    <- datacsv$length # length for each individual
  data$lengthsq  <- datacsv$lengthsq # square of length for each individual
  data$category  <- datacsv$category[datacsv$category!=""] # categories in the data set
  
  return(data)
}

############################################################################### 
###############################################################################
########### Test the function "ic.get.data"  ##################################
# 
# test.ic.get.dat = function(){
#   Data <- ic.get.data("ic.test.data.csv") 
#   return(Data)
# }
# 
# test.ic.get.dat()

###############################################################################