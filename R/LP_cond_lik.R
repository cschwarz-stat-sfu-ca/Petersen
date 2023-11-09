# This function should not be called by the user.
# The function is called twice. First to fit the model for the capture probabilities
# and then to use the fitted model to compute the abundance estimates


#' Compute the conditional likelihood and eventual abundance estimates for LP estimatior
#' using condition likelihood
#' @param p_beta Initial or final p_beta  values
#' @template param.data
#' @template param.p_model
#' @template param.N_hat
#' @param p_beta_vcov Variance matrix for the p_beta from the fit
#' @param what.return Indicate if return the negative conditional log-likelihood,
#' or estimates of model and abundance
#' @template param.conf_level
#' @param trace If debugging should be done.
#'
#' @importFrom stats model.matrix coef qnorm
#' @importFrom plyr ddply llply


#' @noRd




LP_cond_lik <- function(p_beta, data, p_model=NULL,
                        N_hat=~1,
                        p_beta_vcov=NULL,
                        what.return=c("negcll","estimates")[1],
                        conf_level=0.95,
                         trace=FALSE){
  # compute the conditional likelihood. See Alho (1990) for details
  N_hat_f = N_hat
  # null some "hidden" variables to avoid R CMD check errors. See
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  #cat("Call with ", p_beta, "    ", Sys.time(), "\n")
  p <- sum.p <- prod.p <- p.hist <- K <- freq <- NULL
  if(trace)browser()
  # p_beta is on the logit scale

  # add an index to the capture histories
  data$..index <- 1:nrow(data)

  # extract the capture history vector from the data$cap_hist
  cap_hist <- split_cap_hist(data$cap_hist)
  n.cap = ncol(cap_hist)
  if(n.cap !=2)stop("Only use for 2 sampling events")

  # expand the data data.frame with each row repeated length(p_beta) times
  data.expand <- data[ rep(1:nrow(data), each=n.cap),]
  data.expand$..time         <- as.character(1:n.cap) # ensure that treated as factor
  data.expand$cap_hist_indiv <- as.vector(t(cap_hist))

  # get the model matrix for the time points and apply the p_beta
  data.expand$p.logit <- as.vector(model.matrix(p_model, data=data.expand) %*% p_beta)
  data.expand$p.odds  <- exp(data.expand$p.logit)
  data.expand$p       <- expit(data.expand$p.logit)
  data.expand$p.hist  <- (data.expand$cap_hist_indiv=='1')*(  data.expand$p) +
    (data.expand$cap_hist_indiv=='0')*(1-data.expand$p)

  # get the (1-p_i) terms for each time point
  p1.index <- seq(1, nrow(data.expand),2)
  data.expand$..1mp <- 1
  data.expand$..1mp[  p1.index] <- 1- data.expand$p[1+p1.index]
  data.expand$..1mp[1+p1.index] <- 1- data.expand$p[  p1.index]

  data.collapse <- data.frame(
                      sum.p = data.expand$p[p1.index] + data.expand$p[1+p1.index],
                      prod.p= data.expand$p[p1.index] * data.expand$p[1+p1.index])
  data.collapse$K          = data.collapse$sum.p - data.collapse$prod.p
  data.collapse$p.hist.cond= data.expand$p.hist[p1.index]*data.expand$p.hist[1+p1.index]/ data.collapse$K
  data.collapse$freq       = data.expand$freq  [p1.index]

  # find the conditional factor (K(theta)) in Alho's paper
  # now to compute the conditional log-likelihood (see Alho p.626)
  cll <- sum(data.collapse$freq*log(data.collapse$p.hist.cond))
  # finally the conditional log-likelihood
  neg_cll <- - cll
  if(what.return=="negcll")return(neg_cll)


  # compute parts of the second term of the variance
  # we have p(i,t)=expit( p*=model.matrix(p) %*% p_beta)
  # So we need the partial of p wrt to p_beta
  # This is computed in two parts
  #.     dp/dp_beta = dp/dp* dp*/dp_beta
  # The first term  is
  #.   d expit(p*)= d 1/(1+exp(-p*))= exp(-p*/(1+exp(-p*))^2)= p(1-p).
  #.The second term is simply the model matrix

  # We need to avoid creating huge matrices (e.g., via the diagonal term)
  # see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
  # m * rep(v, rep(nrow(m),length(v))),
  if(trace)browser()
  #partial.wrt.p_beta.old = diag(as.vector(data.expand$p *(1-data.expand$p))) %*%
  #                     model.matrix(p_model, data=data.expand)

  partial.wrt.p_beta =  model.matrix(p_model, data=data.expand)
  temp <- as.vector(data.expand$p *(1-data.expand$p))
  temp <- rep(temp, length.out=length(partial.wrt.p_beta))
  partial.wrt.p_beta <- partial.wrt.p_beta * temp

  # now for the (1-p of other terms) part of equation on bottom of page 6 of Lee paper
  #partial.wrt.p_beta = diag(as.vector(data.expand$..1mp)) %*% partial.wrt.p_beta

  temp <- as.vector(data.expand$..1mp)
  temp <- rep(temp, length.out=length(partial.wrt.p_beta))
  partial.wrt.p_beta <- partial.wrt.p_beta * temp

  # add up the two partials for p1 or p2
  #temp <- diag(nrow(data))
  #temp <- temp[rep(1:nrow(data),each=n.cap),]
  #partial.wrt.p_beta  <- t(temp) %*% partial.wrt.p_beta  # add up the two capture times
  partial.wrt.p_beta <- partial.wrt.p_beta[p1.index,] + partial.wrt.p_beta[1+p1.index,]


  # multiply by frequency
  #partial.wrt.p_beta <- diag(as.vector(data$freq)) %*%
  #  diag(as.vector(-1/data.collapse$K^2)) %*%
  #  partial.wrt.p_beta
  temp <- as.vector(data$freq * -1/data.collapse$K^2)
  temp <- rep(temp, length.out=length(partial.wrt.p_beta))
  partial.wrt.p_beta <- partial.wrt.p_beta * temp

  data.collapse$..EF <- 1/data.collapse$K
  data <- cbind(data, data.collapse[,c("p.hist.cond","K", "..EF")])

  # get the estimates of N as defined by the formula for N
    # HT estimates of N
    terms <- model.matrix(N_hat_f, data)
    N_hat <- t(terms) %*% (data$freq * data$..EF)

    # see Huggings equation for s2 on p.135
    # this is the first component of variance that assumes that the HT denominator are known exactly
    N_hat_var1 <- diag(as.vector(t(terms) %*% ( data$freq * 1/data$K^2 *(1-data$K))), ncol=ncol(terms))
    if(trace)browser()
    temp <- t(terms) %*% partial.wrt.p_beta
    N_hat_var2 <- temp %*% p_beta_vcov %*% t(temp)
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
             cond.ll     = cll,   # conditional log likelihood
             n.parms     = length(p_beta), # number of parameters
             N_hat       =N_hat_res,
             method      ="CondLik",
             DateTime    =Sys.time())
  return(res)

}
