###############################################################################
#####                           get data for                             ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ######
#####                     with known dead removals                       ###### 
###############################################################################

# Function "PSkSCR.get.data"
#   Input: Data Set in Raw CSV form. 
#   Ouput: Data Set(as a list) ready for analysis. 

PSkSCR.get.data = function(testdata.csv){
  
  datacsv = read.csv(testdata.csv, header = TRUE, as.is=TRUE, strip.white=TRUE)
  
  data = NULL
  data$history   = datacsv$history  # capture histories of the raw data
  data$counts    = datacsv$counts  # counts of the corresponding capture histories
  data$category  = datacsv$category[datacsv$category!=""] # categories in the data set
  
  return(data)
  
} # End of "PSkSCR.get.data"

###############################################################################