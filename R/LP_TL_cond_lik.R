# This function should not be called by the user.
# The function is called twice. First to fit the model for the capture probabilities
# and then to use the fitted model to compute the abundance estimates


#' Compute the conditional likelihood and eventual abundance estimates for LP estimatior
#' ACCOUNTING FOR TAG LOSS
#' using condition likelihood (conditional on being captured at time 2)
#' @param all_beta Initial or final all_beta  values for modeling p1 (prob capture at t1) and for modelling rho (tag retention)
#' @template param.data
#' @template param.dt_type
#' @template param.p_model
#' @template param.N_hat
#' @param all_beta_vcov Variance matrix for the all_beta from the fit
#' @param what.return Indicate if return the negative conditional log-likelihood,
#' or estimates of model and abundance
#' @template param.conf_level
#' @param trace If debugging should be done.
#'
#' @importFrom stats model.matrix coef qnorm
#' @importFrom plyr ddply llply
#' @importFrom tidyr pivot_wider

#' @noRd




LP_TL_cond_lik <- function(all_beta, data, dt_type, p_model=NULL, rho_model=NULL,
                        N_hat=~1,
                        all_beta_vcov=NULL,
                        what.return=c("negcll","estimates")[1],
                        conf_level=0.95,
                        trace=FALSE){
  # compute the conditional likelihood under tag loss. We condition on
  # recoveries at time 2 only and so can eliminate histories that have
  # indicate no capture at time 2.
  N_hat_f = N_hat
  # null some "hidden" variables to avoid R CMD check errors. See
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  #cat("Call with ", p_beta, "    ", Sys.time(), "\n")
  if(trace)browser()
  # all_beta is on the logit scale

  # add an index to the capture histories
  data$..index <- 1:nrow(data)

  # extract the capture history vector from the data$cap_hist
  cap_hist <- cbind(t1=substr(data$cap_hist, 1,2), t2=substr(data$cap_hist, 3,4)) # this splits the tag history into two elements

  # we need to paste columns 1 and colu
  n.cap = ncol(cap_hist)
  if(n.cap !=2)stop("Only use for 2 sampling events")

  # only want histories where captured at t2.. we condition here
  data.cond   <- data[ cap_hist[,2,drop=TRUE] %in% c("11","10","01","1X", "1P","0P"),]
  data.expand <- data.cond[ rep(1:nrow(data.cond), each=n.cap),]
  data.expand$..tag   <- as.character(1:2) # tag number

  # separate the two parts of the beta matrix (p_beta , followed by rho_beta)
  p1.model.matrix <- model.matrix(p_model, data=data.cond)
  p1_beta <- all_beta[1:ncol(p1.model.matrix)] # extract the beta values for p1

  # get the estimates of p1
  p1.logit <- as.vector(p1.model.matrix %*% p1_beta)
  p1       <- expit(pmax(-10, pmin(10,p1.logit)))

  # get the estimates of rho (tag retention probability)
  rho.model.matrix <- model.matrix(rho_model, data=data.expand)
  rho_beta  <- all_beta[-(1:length(p1_beta))]
  rho.logit <- as.vector(rho.model.matrix %*% rho_beta)
  data.expand$rho       <- expit(pmax(-10,pmin(10,rho.logit)))
  if(dt_type == valid_dt_type()[3]){  # tag2 is permanent
    data.expand$rho[ data.expand$..tag==2] <- 1
  }
  rho <- tidyr::pivot_wider(data.expand,
                            id_cols="..index",
                            names_from="..tag",
                            names_prefix="tag",
                            values_from="rho")
  rho <- rho[ order(rho$..index),] # be sure it is sorted properly

  # compute fraction that are single tagged. This is assumed to be fixed and known
  pST   <- sum(data$freq[ cap_hist[,1,drop=TRUE] %in% c("10","01","0P")])/
           sum(data$freq[ cap_hist[,1,drop=TRUE] %in% c("10","11","01","1P","0P")])
  pST.1 <- sum(data$freq[ cap_hist[,1,drop=TRUE] %in% c("10")])      / sum(data$freq[ cap_hist[,1,drop=TRUE] %in% c("10","01","0P")])
  pST.2 <- sum(data$freq[ cap_hist[,1,drop=TRUE] %in% c("01","0P")]) / sum(data$freq[ cap_hist[,1,drop=TRUE] %in% c("10","01","0P")])
  if(pST==0){
    pST.1 <-0
    pST.2 <-0
  }
  # we now get the conditional probabilities of A., .A, BB, B., .B, U  using notation from Petersen Monograph

  # pA. = probability single tagged recovered at t2 (ignoring the p2 term which cancels when we condition)
  pA.  <- p1 * pST * pST.1 * rho[,"tag1",drop=TRUE]
  # p.A = probability single tagged recovered (tagged at position 2) at t2 (ignoring the p2 term which cancels when we condition)
  p.A  <- p1 * pST * pST.2 * rho[,"tag2",drop=TRUE]
  # pBB= probability double tagged fish recovered at t2 with both tags (ingoring the p2 term)
  pBB <- p1 * (1-pST) * rho[,"tag1",drop=TRUE] * rho[,"tag2",drop=TRUE]
  # pB. = probability double tagged fish lost second tag and recovered at t2 (ignoring the p2 term)
  pB. <- p1 * (1-pST) * rho[,"tag1",drop=TRUE] * (1-rho[,"tag2",drop=TRUE])
  # p.B = probability double tagged fish lost first tag and recovered at t2 (ignoring the p2 term)
  p.B <- p1 * (1-pST) * (1-rho[,"tag1",drop=TRUE]) * rho[,"tag2",drop=TRUE]
  # pBX = probability that double tagged fish lost either t1 or t2 but cannot tell
  pBX = pB. + p.B
  # p.U = prob apparently unmarked fish captured at t2
  pU  <- p1 * pST * pST.1 *(1-rho[,"tag1",drop=TRUE])+ # single tagged, lost tag, recaptured
         p1 * pST * pST.2 *(1-rho[,"tag2",drop=TRUE])+
         p1 * (1-pST) *(1-rho[,"tag1",drop=TRUE])*(1-rho[,"tag2",drop=TRUE])+ # double tagged lost both tags
         (1-p1)

  # force to sum to 1
  sum.p <- pA. + p.A + pBB + pB. + p.B + pU
  pA. <- pA.  / sum.p
  p.A <- p.A  / sum.p
  pBB <- pBB  / sum.p
  pB. <- pB.  / sum.p
  p.B <- p.B  / sum.p
  pBX <- pBX  / sum.p
  pU  <- pU   / sum.p

  #browser()
  # compute the conditional log likelihood
  cll <- 0
  select <- data.cond$cap_hist == "1010" # pA.
  cll <- cll + sum(data.cond$freq[select] * log(pA.[select]), na.rm=TRUE)
  select <- data.cond$cap_hist %in% c("0101", "0P0P") # p.A
  cll <- cll + sum(data.cond$freq[select] * log(p.A[select]), na.rm=TRUE)
  select <- data.cond$cap_hist %in% c("1111","1P1P") # pBB
  cll <- cll + sum(data.cond$freq[select] * log(pBB[select]), na.rm=TRUE)
  select <- data.cond$cap_hist %in% c("1101", "1P0P") # p.B
  cll <- cll + sum(data.cond$freq[select] * log(p.B[select]), na.rm=TRUE)
  select <- data.cond$cap_hist == "1110" # pB.
  cll <- cll + sum(data.cond$freq[select] * log(pB.[select]), na.rm=TRUE)
  select <- data.cond$cap_hist == "111X" # pBX
  cll <- cll + sum(data.cond$freq[select] * log(pBX[select]), na.rm=TRUE)
  select <- data.cond$cap_hist == "0010" # pU
  cll <- cll + sum(data.cond$freq[select] * log(pU[select]), na.rm=TRUE)

  # finally the conditional log-likelihood
  neg_cll <- - cll
  if(what.return=="negcll")return(neg_cll)

  # The abundance is computed as sum( freq/ p1) so 1/p1 is the expansion factor
  p1.model.matrix <- model.matrix(p_model, data=data)

  # get the estimates of p1 and se
  p1.logit <- as.vector(p1.model.matrix %*% p1_beta)
  p1       <- expit(p1.logit)

  p1_beta_vcov <- all_beta_vcov[1:length(p1_beta), 1:length(p1_beta)]
  # https://stackoverflow.com/questions/21708489/compute-only-diagonals-of-matrix-multiplication-in-r
  p1.logit.var  <- rowSums((p1.model.matrix %*% p1_beta_vcov)  * p1.model.matrix)
  p1.var        <- p1.logit.var * p1^2 * (1-p1)^2
  p1.se         <- sqrt(p1.var)

  data<- cbind(data, p1=p1, p1.se=p1.se)
  data$K    <- data$p1
  data$..EF <- 1/data$p1

  #browser()
  # get se of estimates of rho in data.expand
  rho_beta_vcov  <- all_beta_vcov[-(1:length(p1_beta)), -(1:length(p1_beta))]
  rho.logit.var  <- rowSums((rho.model.matrix %*% rho_beta_vcov)  * rho.model.matrix)
  rho.var        <- rho.logit.var * data.expand$rho^2 * (1-data.expand$rho)^2
  data.expand$rho.se <- sqrt(rho.var)


  # restrict the data to having being captured at t1
  data.cond <- data[ substr(data$cap_hist,1,2) != "00",]

  # compute parts of the second term of the variance
  # we have p(i,t)=expit( p*=model.matrix(p) %*% p1_beta)
  # So we need the partial of p wrt to p1_beta
  # This is computed in two parts
  #.     dp/dp1_beta = dp/dp* dp*/dp1_beta
  # The first term  is
  #.   d expit(p*)= d 1/(1+exp(-p*))= exp(-p*/(1+exp(-p*))^2)= p(1-p).
  #.The second term is simply the model matrix

  # We need to avoid creating huge matrices (e.g., via the diagonal term)
  # see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
  # m * rep(v, rep(nrow(m),length(v))),
  if(trace)browser()
  #partial.wrt.p1_beta.old = diag(as.vector(data.expand$p *(1-data.expand$p))) %*%
  #                     model.matrix(p_model, data=data.expand)

  partial.wrt.p1_beta =  model.matrix(p_model, data=data.cond)
  temp <- as.vector(data.cond$p1 *(1-data.cond$p1))
  temp <- rep(temp, length.out=length(partial.wrt.p1_beta))
  partial.wrt.p1_beta <- partial.wrt.p1_beta * temp


  # multiply by frequency
  #partial.wrt.p1_beta <- diag(as.vector(data$freq)) %*%
  #  diag(as.vector(-1/data.collapse$K^2)) %*%
  #  partial.wrt.p1_beta
  temp <- as.vector(data.cond$freq * -1/data.cond$K^2)
  temp <- rep(temp, length.out=length(partial.wrt.p1_beta))
  partial.wrt.p1_beta <- partial.wrt.p1_beta * temp



  # get the estimates of N as defined by the formula for N
    # HT estimates of N
    terms <- model.matrix(N_hat_f, data.cond)
    N_hat <- t(terms) %*% (data.cond$freq * data.cond$..EF)

    # see Huggings equation for s2 on p.135
    # this is the first component of variance that assumes that the HT denominator are known exactly
    N_hat_var1 <- diag(as.vector(t(terms) %*% ( data.cond$freq * 1/data.cond$K^2 *(1-data.cond$K))), ncol=ncol(terms))
    if(trace)browser()
    temp <- t(terms) %*% partial.wrt.p1_beta
    N_hat_var2 <- temp %*% p1_beta_vcov %*% t(temp)
    #browser()
    N_hat_vcov <- N_hat_var1 + N_hat_var2

    N_hat_SE = sqrt(diag(N_hat_vcov))

    # find the confidence interval using a log() transformation
    z <- qnorm(1-(1-conf_level)/2)

    log_N_hat    = log(N_hat)
    log_N_hat_SE = N_hat_SE / N_hat
    lcl = exp( log_N_hat -z*log_N_hat_SE)
    ucl = exp( log_N_hat +z*log_N_hat_SE)

    N_hat_res <- list(N_hat_f = N_hat_f,
         N_hat=N_hat, N_hat_vcov=N_hat_vcov, N_hat_SE=N_hat_SE,
         N_hat_conf_level  = conf_level,
         N_hat_conf_method = "logN",
         N_hat_LCL         = lcl,
         N_hat_UCL         = ucl
    )

  res<- list(data        = data,
             data.expand = data.expand,
             p_model     = p_model,
             rho_model    = rho_model,
             cond.ll     = cll,   # conditional log likelihood
             n.parms     = length(p1_beta)+ length(rho_beta), # number of parameters
             N_hat       =N_hat_res,
             method      ="LP_TS CondLik",
             DateTime    =Sys.time())
  return(res)

}



