# Create the various design matrices.


# The following function creates a design matrix using forcing vector informations; 
#
#   For Example: Assume Forcing Vector is 1 ,1, 3, 4,
#   The Design Matrix becomes:
#     
#        1     0     0
#        1     0     0
#        0     1     0
#        0     0     1
# 
# crate.DM.
# Input: Vector of length(2*categories)
# Ouput: Matrix of dimention :  length(2*categories)x(Unique numbers in input vector) 


create.DM = function(info){
  # Create the design matrix for row or column pooling
  un  <- unique(info)
  m   <- as.matrix(diag(length(un)))
  res <- m[ match(info, un), ,drop=FALSE]
  return(res)
}

######################################################
#Testing above functions
create.DM(c(1,1,2,3))
create.DM(c(1,1,2,2))
create.DM(c(1,1,1,1))
create.DM(c(1,2))
