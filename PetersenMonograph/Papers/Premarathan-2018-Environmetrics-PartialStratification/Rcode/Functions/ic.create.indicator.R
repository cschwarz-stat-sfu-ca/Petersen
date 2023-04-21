###############################################################################
#####                    Create indicator variables for                  ######
##### Partial Stratification in two-Sample Capture-Recapture Experiments ######
###############################################################################

# Function : ic.create.indicator
#  Input : data
#  Output : indicator variables

ic.create.indicator <- function(Data){
  
  History <- Data$History
  counts  <- Data$counts
  cats    <- Data$category
  
  ##number of sampling occasions
  st <-  2
  
  ### create an  indicator variable "h" to identify the capture history
  ###  1 - captured
  ###  0 - Not not captured
  h <- matrix(0, nrow=length(History), ncol=st) 
  for(i in 1:st){
    h[,i] <- substr(History,i,i)!=0
  }
  
  
  ### creating indicator variable "cat.ind" to identify the category
  ### 0 - indicates not captured i.e. "0" 
  ### 1 - represents cats[1], i.e "M" 
  ### 2 - represents cats[2], i.e. "F" and so on..
  ### 
  ### largest number represents "U"
  
  cat.ind <- matrix(0, nrow=length(History), ncol=st) 
  for(i in 1:st){
    cat.ind[,i] <- match(substr(History,i,i),c(cats,"U"), nomatch = 0 )
  }
  
  
  ## create an indicator variable "s" to identify whether the animal is 
  ## stratified (sub-sampled) or not at the current sampling time.
  ###  1 - stratified
  ###  0 - Not stratified
  s <- matrix(0, nrow=length(History), ncol=st) 
  for(j in 1:length(History)){
    temp <- 1
    while( temp <= st){
      if((cat.ind[j,temp]!=0) && (cat.ind[j,temp]!= (length(cats)+1))){
        s[j,temp] <- 1
        break
      }
      temp <- temp+1
    }
  }
  
  
  ## create an indicator variable "cs" to identify whether the animal is 
  ## captured in the sample or not at the current sampling time from the 
  ## animals not marked earlier. 
  ###  1 - captured in the sample ( i.e. "U", "M", "F" )
  ###  0 - Not captured  in the sample ( "0")
  cs <- matrix(0, nrow=length(History), ncol=st) 
  for(j in 1:length(History)){
    temp <- 1
    while( temp <= st){
      if(cat.ind[j,temp]!=0){
        cs[j,temp] <- 1
        break
      }
      temp <- temp+1
    }
  }
  
  
  ## Create an indicator variable "s.ind" to indicate whether the captured animal  
  ## is in particular category or "U"
  ## "s.ind"  is a matrix with length(cats) columns and length(history rows)
  #### column 1 indicates the captures animal is in cats[1] or "U" by 1 and 
  ####    otherwise 0 ( ie. "M" or "U")
  #### column 2 indicates the captures animal is in cats[2] or "U" by 1 and 
  ####    otherwise 0( ie. "F" or "U")
  ####
  ####  and so on.....
  
  s.ind <- matrix(0, nrow=length(History), ncol=length(cats)) 
  for(j in 1:length(History)){
    for(i in 1:length(cats)){
      temp <- 1
      while( temp <= st){
        if(cat.ind[j,temp]==i ||cat.ind[j,temp]==3){
          s.ind[j,i] <- 1
          break
        }
        temp <- temp+1
      }
    }
  }
  
  

  indicator <- NULL
  indicator $h <- h
  indicator$cat.ind <- cat.ind
  indicator$s <- s
  indicator$cs <- cs
  indicator$s.ind <- s.ind

    return(indicator)
  
} # end of the function "ic.create.indicator"

###############################################################################
###############################################################################


###############################################################################
#####  function to test "ic.create.indicator" #################################
# 
# test.ic.create.indicator <- function(){
# 
#    Data <- NULL
#    Data$History <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U","MM","M0","0F","U0","F0")
#    Data$counts   <- c(41,9,16,4,23,7,12,25,63,2,1,3,8,2)
#    Data$category <- c("M","F")
#    str(Data)
#    print(Data)
#    
#    
#    out <- ic.create.indicator(Data)
#    
#    print(out)
# }
# 
# test.ic.create.indicator()

######################### End of test functions## #############################
###############################################################################