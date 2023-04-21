# Function "get.data"
#   Input: Data Set in Raw CSV form. 
#   Ouput: Data Set(as a list) ready for analysis. 

get.data = function(testdata.csv){
  
  datacsv = read.csv(testdata.csv, header = TRUE, as.is=TRUE, strip.white=TRUE)
  
  Data = NULL
  Data$History   = datacsv$CaptureHistory  # capture histories of the raw data
  Data$counts    = datacsv$Counts  # counts of the corresponding capture histories
  Data$category  = datacsv$Category[datacsv$Category!=""] # categories in the data set
  
  return(Data)
}


