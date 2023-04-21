###############################################################################
#####                    Create indicator variables for                  ######
##### Partial Stratification in k-Sample Capture-Recapture Experiments   ###### 
#####                     with known dead removals                       ######
###############################################################################

# Function : PSkSCR.create.indicator
#  Input : data
#  Output : indicator variables

PSkSCR.create.indicator <- function(data){
  
  history <- data$history
  count  <- data$count
  cats    <- data$category
  
  ##number of sampling occasions
  st <-  nchar(history[1])
  
  
  ### create an  indicator variable "h" to identify the capture history
  ###  1 - captured
  ###  0 - Not not captured
  h <- matrix(0, nrow=length(history), ncol=st) 
  for(i in 1:st){
    h[,i] <- substr(history,i,i)!=0
  }
  
  
  ## creat an indicator variable "z" to identify whether the animal is 
  ## available to capture at the current sampling time. animal may 
  ## alive of dead at the current sampling time
  ###  1 - available to capture 
  ###  0 - Not available to capture
  z <- matrix(1, nrow=length(history), ncol=st) 
  for(j in 1:length(history)){
    temp <- st
    if(count[j] < 0){
      while(h[j,temp] !=1){
        z[j,temp] <- 0
        temp <- temp-1
      }
    }
  }
  
  
  ### creating indicator variable "cat.ind" to identify the category
  ### 0 - indicates not stratified i.e. "0" or "U"
  ### 1 - represents cats[1], i.e "M" 
  ### 2 - represents cats[2], i.e. "F" and so on..
  ### 
  ### largest number represents "U"
  
  cat.ind <- matrix(0, nrow=length(data$history), ncol=st) 
  for(i in 1:st){
    cat.ind[,i] <- match(substr(data$history,i,i),c(cats,"U"), nomatch = 0 )
  }
  
  
  ## create an indicator variable "ss" to identify whether the animal is 
  ## stratified (sub-sampled) or not at the current sampling time.
  ###  1 - stratified
  ###  0 - Not stratified
  s <- matrix(0, nrow=length(history), ncol=st) 
  for(j in 1:length(history)){
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
  cs <- matrix(0, nrow=length(history), ncol=st) 
  for(j in 1:length(history)){
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
  #### and so on...
  ##
  s.ind <- matrix(0, nrow=length(history), ncol=length(cats)) 
  for(j in 1:length(history)){
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
  
  
  ## creat an indicator variable "dead" to identify whether the 
  ## captured animal is alive or dead
  ###  1 - animal is alive
  ###  0 -animal is dead
  dead <- matrix(1, nrow=length(data$history), ncol=st) 
  for(j in 1:length(data$history)){
    temp <- st
    if(data$counts[j] < 0){
      while(h[j,temp] !=1){
        dead[j,temp] <- 0
        temp <- temp-1
      }
      dead[j,temp] <- 0
    }
  }
  
  
  
  indicator <- NULL
  indicator $h <- h
  indicator$z <- z
  indicator$cat.ind <- cat.ind
  indicator$s <- s
  indicator$cs <- cs
  indicator$s.ind <- s.ind
  indicator$dead <- dead
  
  return(indicator)
  
} # end of the function "PSkSCR.create.indicator"

###############################################################################
###############################################################################


###############################################################################
#####  function to test "PSkSCR.create.indicator" #############################

# test.PSkSCR.create.indicator <- function(){
# 
#    data <- NULL
#    data$history <- c("U0","UU","M0","MM","F0","FF","0M","0F","0U","MM","M0","0F","U0","F0")
#    data$count   <- c(41,9,16,4,23,7,12,25,63,-2,-1,-3,-8,-2)
#    data$category <- c("M","F")
#    str(data)
#    print(data)
#    
#    
#    out <- PSkSCR.create.indicator(data)
#    
#    print(out)
# }
# 
# test.PSkSCR.create.indicator()

######################### End of test functions## #############################
###############################################################################
