fun_get_freq <- function(Enc.Mat,data.type="LR"){
  #This function calculates an initial frequency count from the feasible set of latent history frequencies (x_0 in paper) by assuming a 1-to-1 mapping of the recorded histories to latent histories.
  # For the LR and LRB data types, the index of latent history w_i is i = 1+ sum_t (w_t*4^(T-t)) (w_t being a 0, 1, 2, or 3 to represent
  # a non-encounter, a left-sided encounter, a right-sided encounter, and both-sided encounter, respectively).
  # For the LRAB data type, the index of latent history w_i is i = 1+ sum_t (w_t*5^(T-t)) (w_t being a 0, 1, 2, 3, or 4 to represent
  # a non-encounter, a left-sided encounter, a right-sided encounter, left- and right-sided encounter, and simultaneous both-sided encounter, respectively).
  
  # Arguments: 
  # Enc.Mat is a matrix of recorded histories with each row representing an individual and each column representing a sampling occasion.
  #   Entries for Enc.Mat are 0, 1, or 2 for the LR data type and 0, 1, 2, or 3 for the LRB data type.
  # data.type is the data type that determines mapping of recorded histories to latent histories (see Table 1 in paper). 
  #   LR data type indicated by "LR"; LRB data type indicated by "LRB"; LRAB data type indicated by "LRAB".
  
  if(data.type=="LR"){
    ind=4
    if(any(!dunif(Enc.Mat,0,2))) stop("LR encounter histories can only include 0, 1, and 2 entries")
  } else if(data.type=="LRB"){
    ind=4
    if(any(!dunif(Enc.Mat,0,3))) stop("LRB encounter histories can only include 0, 1, 2, and 3 entries")
    if(!any(Enc.Mat==3)) cat("Warning: LRB recorded histories contain no both-sided encounters -- should you be using the LR data type?")
  } else if(data.type=="LRAB"){
    ind=5
    if(any(!dunif(Enc.Mat,0,4))) stop("LRAB encounter histories can only include 0, 1, 2, 3, and 4 entries")
    if(!any(Enc.Mat==4)) cat("Warning: LRAB recorded histories contain no both-sided encounters -- should you be using the LR data type?")
  }
  temp.noccas=ncol(Enc.Mat)
  Hist.num=rep(0,nrow(Enc.Mat))
  Value=ind^c((temp.noccas-1):0)
  for(i in 1:nrow(Enc.Mat)){
    Hist.num[i]=1+Enc.Mat[i,] %*% Value      #assign a history index i to each latent history, where i = 1+ \sum_{t=1}^T \omega_t ind^{T-t}
  }
  #frequency counts
  temp.Freq=rep(0,ind^temp.noccas)
  for(i in 1:nrow(Enc.Mat)){
    temp.Freq[Hist.num[i]]=temp.Freq[Hist.num[i]]+1
  }
  temp.Freq
}

fun_get_A <- function(noccas,data.type="LR"){
  #This function constructs A, for f=A'x.  Function returns A matrix, latent histories (All.hists), and indices for the latent histories spawning >1 recorded history (ivect).
  
  # Arguments: 
  # noccas = number of sampling occasions (T in paper)
  # data.type = data type that determines mapping of recorded histories to latent histories (see Table 1 in paper). 
  #   LR data type indicated by "LR"; LRB data type indicated by "LRB"; LRAB data type indicated by "LRAB"
  
  if(data.type=="LR"){
    ind=4
    r=4^noccas-2*(2^noccas-1)          #r is the number of free variables (i.e., the number of basis vectors)
  } else if(data.type=="LRB"){
    ind=4
    r=3^noccas-2^(noccas+1)+2
  } else if(data.type=="LRAB"){
    ind=5
    r=4^noccas-2*(2^noccas-1) 
  } else{stop("data type must be LR, LRB, or LRAB")}
  
  cat("Dimension of null space is ",ind^noccas," x ",ind^noccas-r," = ",ind^noccas*(ind^noccas-r))

  
  #first, create correponding matrix w/all possible latent histories using recursive algorithm
  dimAllhists=(ind^noccas)*noccas
  nonzeroAllhists=ind^(noccas-1)*noccas*(ind-1)
  indexAllhists=integer(nonzeroAllhists)
  valueAllhists=integer(nonzeroAllhists)
  
  indexcount=0
  for(i in 1:noccas){
    if(data.type=="LRAB"){
      series=rep(c(rep(0,5^(noccas-i)),rep(1,5^(noccas-i)),rep(2,5^(noccas-i)),rep(3,5^(noccas-i)),rep(4,5^(noccas-i))),5^(i-1))      
    } else {
      series=rep(c(rep(0,4^(noccas-i)),rep(1,4^(noccas-i)),rep(2,4^(noccas-i)),rep(3,4^(noccas-i))),4^(i-1))
    }
    seriesind=which(series>0)
    indexAllhists[indexcount+1:length(seriesind)]=seriesind+(i-1)*dimAllhists/noccas;
    valueAllhists[indexcount+1:length(seriesind)]=series[seriesind];
    indexcount=indexcount+length(seriesind)
  }
  
  All.hists=Matrix(0,nrow=(ind^noccas),ncol=noccas)
  All.hists[indexAllhists]=valueAllhists
  
  Value=ind^c((noccas-1):0)
  
  # Construct A matrix   
  ivect=which((rowSums(All.hists==1)>0 & rowSums(All.hists==2)>0 & (rowSums(All.hists==3)==0 & rowSums(All.hists==4)==0)) | (data.type!="LRB")*(rowSums(All.hists==3)>0 & rowSums(All.hists==4)==0))
  temp.hist=All.hists[ivect,]
  temp.R=Matrix(0,nrow=length(ivect),ncol=noccas)
  temp.L=Matrix(0,nrow=length(ivect),ncol=noccas)
  temp.R[which(temp.hist==2 | ((temp.hist==3)*(data.type!="LRB")))]=2
  temp.L[which(temp.hist==1 | ((temp.hist==3)*(data.type!="LRB")))]=1
  placesR=as.vector(temp.R %*% Value)
  placesL=as.vector(temp.L %*% Value) 
  
  A=sparseMatrix(i=c(ivect,ivect),j=c(placesL+1,placesR+1),dims=c(ind^noccas,ind^noccas),x=1)
  diag(A)[-ivect]=1
  A=A[,-1]
  A=A[,-which(colSums(A)==0)]
  
  out=list(A=A,All.hists=All.hists,ivect=ivect)
}

fun_get_basis_vectors <- function(noccas,A,ivect,Freq,div=3,data.type="LR"){
  #This function caculates basis vectors based on data type (data.type) and latent frequencies (Freq).  
  #Function returns a matrix of the relevant basis vectors for the null space of A'.
  
  # Arguments: 
  # noccas = number of sampling occasions (T in paper)
  # A and ivect are list objects returned by "fun_get_A" above.
  # Freq = initial frequencies for latent histories (x_0 in paper) return by "fun_get_freq" above.
  # div = integer for dividing matrix manipulation workload for very large null spaces. If experiencing matrix size issues, try larger values for div.
  # data.type = data type that determines mapping of recorded histories to latent histories (see Table 1 in paper). LR data type indicated by "LR"; LRB data type indicated by "LRB"; LRAB data type indicated by "LRAB".
  
  #determine basis vectors
  if(data.type=="LR"){
    ind=4
    r=4^noccas-2*(2^noccas-1)                   #r is the number of free variables (i.e., the number of basis vectors)
  } else if(data.type=="LRB"){
    ind=4
    r=3^noccas-2^(noccas+1)+2
  } else if(data.type=="LRAB"){
    ind=5
    r=4^noccas-2*(2^noccas-1) 
  }
  
  if(noccas>7){
    cat("This might be a while...that's a darned big null space!")
  }
  
  Basis=Matrix(0,ind^noccas,r)
  
  free=c(1,ivect)   # indices for the free variables (i.e., the latent histories that spawn >1 recorded history)
  
  bound=seq(1:(ind^noccas))[-free]  # indices for the bound variables (i.e., the latent histories that spawn only 1 recorded history)
  
  k=length(bound)
  
  div_r = ceiling(seq(2,r,length=ifelse(noccas>2,div,r-1)))
  temp = t(A)
  for(i in 1:(length(div_r)-1)){
    Basis[bound,(div_r[i]):(div_r[i+1])] = -temp[,free[(div_r[i]):(div_r[i+1])]]   
  }
  
  freemat=Matrix(0,ind^noccas,ind^noccas)
  diag(freemat)=1
  Basis=Basis+freemat[,-bound]
  
  #Given the initial frequency vector 'Freqs', eliminate Basis vectors which always produce negative frequencies.  
  #For example, with T=3, if latent history \omega_2 = '00L' has frequency x_2 = 0, then any basis vector with a '-1' in the second row will produce negative frequencies and can be eliminated. For the LR and LRAB data types, one would eliminate basis vectors 2, 5, 6, 23, 24, 29, and 30. For the LRB data type, one would eliminate basis vectors 3, 9, and 13.
  Freqs=c(1,which(Freq[-1]>0)+1)
  temp=sparseMatrix(i=rep(Freqs,r),j=rep(seq(1,r),each=length(Freqs)),dims=c(ind^noccas,r))
  temp=temp+Basis
  temp=which(colSums(temp<0)<1)
  Basis=Basis[,temp]
  return(Basis)
}

fun_sim_data <- function(N=100,noccas=5,p=0.4,delta_R=0.4,delta_L=0.4,alpha=0.5,data.type="LR"){
  #This function generates ecounter histories for bilateral photo-ID data  
  
  # Arguments: 
  # N = population size
  # noccas = number of sampling occasions
  # p = probability of detection
  # delta_R = probability of right-side encounter, given detection
  # delta_L = probability of left-side encounter, given detection
  # alpha = probability of both-sided detection, given both sides were encountered (alpha=0 for LR data type and alpha=1 for LRB data type)
  # data.type = data type that determines mapping of recorded histories to latent histories (see Table 1 in paper). LR data type indicated by "LR"; LRB data type indicated by "LRB"; LRAB data type indicated by "LRAB".
  
  delta_B=1-delta_R-delta_L
  
  if(data.type=="LR"){
    alpha=0
  } else if(data.type=="LRB") {
    alpha=1
  }
  
  tEnc.Mat=matrix(rbinom(N*noccas,1,p),N,noccas)        #"true" latent histories
  Rand.Mat=matrix(runif(N*noccas,0,1),N,noccas)
  tEnc.Mat[which(tEnc.Mat==1 & Rand.Mat<delta_R)]=2      # Right side encounters
  tEnc.Mat[which(tEnc.Mat==1 & Rand.Mat>(1-delta_B))]=4  # Both sided encounters
  tEnc.Mat[which(tEnc.Mat==4)]=tEnc.Mat[which(tEnc.Mat==4)]-(runif(sum(tEnc.Mat==4))<(1-alpha))   # unobserved both sided encounters
  
  Enc.Mat=tEnc.Mat
  z.Enc.Mat=which(rowSums(Enc.Mat)==0)
  if(length(z.Enc.Mat)) {
    Enc.Mat=Enc.Mat[-z.Enc.Mat,]        #remove all zero simulated histories
  }
  
  if(data.type=="LR"){
    for(i in 1:nrow(Enc.Mat)){      #histories with R and L or B result in 2 'ghost' histories; change R's to 0 and B's to 1 to get L ghosts, then change B's to 2 and append with R ghosts to end of array
      if((sum(Enc.Mat[i,]==1)>0 & sum(Enc.Mat[i,]==2)>0) | sum(Enc.Mat[i,]==3)>0) {
        curR=which(Enc.Mat[i,]==2)
        curB1=which(Enc.Mat[i,]==3)
        Enc.Mat[i,curR]=0
        Enc.Mat[i,curB1]=1
        new.hist=rep(0,noccas)
        new.hist[curR]=2
        new.hist[curB1]=2
        Enc.Mat=rbind(Enc.Mat,new.hist)
      }
    }
  } else if(data.type=="LRB"){
    for(i in 1:nrow(Enc.Mat)){      #histories with R and L but no B result in 2 'ghost' histories; change R's to 0 to get L ghosts and append R ghosts to end of array
      if(sum(Enc.Mat[i,]==1)>0 & sum(Enc.Mat[i,]==2)>0 & sum(Enc.Mat[i,]==4)==0) {
        curR=which(Enc.Mat[i,]==2)
        curA=which(Enc.Mat[i,]==3)
        Enc.Mat[i,curR]=0
        Enc.Mat[i,curA]=1
        new.hist=rep(0,noccas)
        new.hist[curR]=2
        new.hist[curA]=2
        Enc.Mat=rbind(Enc.Mat,new.hist)
      }
    }
    Enc.Mat[which(Enc.Mat==4)]=Enc.Mat[which(Enc.Mat==4)]-1
  } else if(data.type=="LRAB"){
    for(i in 1:nrow(Enc.Mat)){      #histories with R and L or A but no B result in 2 'ghost' histories; change R's to 0 and A's to 1 to get L ghosts, then change A's to 2 and append with R ghosts to end of array
      if(sum(Enc.Mat[i,]==1)>0 & sum(Enc.Mat[i,]==2)>0 & sum(Enc.Mat[i,]==4)==0 | (sum(Enc.Mat[i,]==3)>0 & sum(Enc.Mat[i,]==4)==0)) {
        curR=which(Enc.Mat[i,]==2)
        curA=which(Enc.Mat[i,]==3)
        Enc.Mat[i,curR]=0
        Enc.Mat[i,curA]=1
        new.hist=rep(0,noccas)
        new.hist[curR]=2
        new.hist[curA]=2
        Enc.Mat=rbind(Enc.Mat,new.hist)
      }
    }
  }
  return(Enc.Mat)
}