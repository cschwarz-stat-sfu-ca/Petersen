#' Simulate data from a Lincoln-Petersen Model with Tag Loss
#'
#' This function creates simulated capture histories for the Lincoln-Petersen
#' model with tag loss.
#'

#' @template param.dt_type
#' @param N Population size
#' @param cov1 Function to generate first covariate for each member of population as function of \code{N}
#' @param cov2 Function to generate second covariate for each member of population as function of \code{cov1}.
#' @param p1 Function to generate P(capture) at event 1 for each member of population as function of \code{cov1,cov2}.
#' @param pST Function to generate P(single tag) if captured at event 1 as function of \code{cov1,cov2}.
#' @param pST.1 Function to generate p(apply single tag to first position at event 1) as function of \code{cov1,cov2}.
#' @param rho1 Function to generate P(tag1 retained) as function of \code{cov1,cov2}.
#' @param rho2 Function to generate P(tag2 retained) as function of \code{cov1,cov2}.
#' @param p2 Function to generate P(capture) at event 2 for each member of population as function of \code{cov1,cov2}.
#' @param trace Trace flag to help debug if things fail.
#' @param seed Initial value of random seed

#' @details
#' The \code{cov1} function takes the value \code{N} and returns N covariate values. For example these could be
#' simulated length, or sex of each fish. The \code{cov2} function takes the \code{cov1} values and generates
#' a second covariate. Two covariates should be sufficient for most capture-recapture simulations.
#' If generating continuous covariates, you should round the covariate to
#' about 100 distinct values to speed up your simulation.
#'
#' The remaining functions take the two covariate values and generate capture probabilities, single tag probabilities,
#' placing tags on fish, and tag retention probabilities. These should all be in the range of 0 to 1.
#'
#' After generating capture histories for the entire population, animals never seen are "discarded" and the
#' data set is compress to unique combinations of the two covariates and the capture history with the frequency
#' variable set accordingly.
#'
#' @return Data frame with observed capture histories
#'
#'
#' @importFrom formula.tools is.one.sided
#' @importFrom plyr is.formula
#' @importFrom stats as.formula model.matrix coef runif
#' @importFrom bbmle mle2

#' @examples
#'
#' sim_data <-LP_TL_simulate(
#'       dt_type="t2perm",  # permanent tag
#'       N=1000,
#'       cov1=function(N)         {rep(1,N)},
#'       cov2=function(cov1)      {rep(1,  length(cov1))},
#'       p1  =function(cov1, cov2){rep(.1, length(cov1))},
#'       pST =function(cov1, cov2){rep(.25,length(cov1))},
#'       rho1=function(cov1, cov2){rep(.70,length(cov1))},
#'       rho2=function(cov1, cov2){rep(1,  length(cov1))},  # permanent second tag
#'       p2  =function(cov1, cov2){rep(.1, length(cov1))},
#'       seed=round(1000000*runif(1)))
#' sim_data
#'
#' @export LP_TL_simulate
#'

LP_TL_simulate <- function(
                     dt_type=NULL,
                     N=1000,
                     cov1 = function(N){rep(1,N)},
                     cov2 = function(cov1){rep(1,length(cov1))},
                     p1   = function(cov1, cov2){rep(.1,length(cov1))},
                     pST  = function(cov1, cov2){rep(.5,length(cov1))},
                     pST.1= function(cov1, cov2){rep(1 ,length(cov1))},
                     rho1 = function(cov1, cov2){rep(.8,length(cov1))},
                     rho2 = function(cov1, cov2){rep(.8,length(cov1))},
                     p2   = function(cov1, cov2){rep(.1,length(cov1))},
                     seed=round(100000000*runif(1)),
                     trace=FALSE){

  # check that arguments are ok
  check.numeric(seed, min.value=1, check.whole=TRUE)
  set.seed(seed)

  if(trace)browser()
  # generate population with two covariates
  check.numeric(N, min.value=100, check.whole=TRUE)
  if(N > 500000)warning("N is 500,000+. Are you sure!")

  check.func <- function(func, func.name="", req.n.arg=2, arg1=pop$cov1, arg2=pop$cov2){
    if(!is.function(func))              stop(func.name, " must be a function")
    if(length(formals(func))!=req.n.arg)stop(func.name, "  must have ", req.n.arg," argument(s)")
    #browser()
    temp <-tryCatch( do.call(func, plyr::compact(list(arg1, arg2))),
        error=function(cond) {
            stop(func.name, " function did not evaluate properly")
            return(NA)
        },
        warning=function(cond) {
            stop(func.name, " function did not evaluate properly")
            return(NULL)
        }
    )
    #browser()
    if(length(temp) != N)stop(func.name, " did not return N values. It returned ", length(temp), " values.")
    temp
  }

  # generate initial covariate
  temp <- check.func(cov1, "cov1", 1, arg1=N, arg2=NULL)
  if(length(unique(temp))>100)warning("More than 100 unique values for cov1. Your simulation will run faster if you round this covariate")
  pop <- data.frame(cov1=temp)

  temp <- check.func(cov2, "cov2", 1, arg1=pop$cov1, arg2=NULL)
  if(length(unique(temp))>100)warning("More than 100 unique values for cov2. Your simulation will run faster if you round this covariate")
  pop$cov2 <- temp

  # generate p1 and p2
  temp <- check.func(p1, "p1", 2, arg1=pop$cov1, arg2=pop$cov2)
  if(!all(temp >=0 & temp <=1))stop("p1 must be between 0 and 1. Your range is ", range(temp))
  pop$p1 <- temp

  temp <- check.func(p2, "p2", 2, arg1=pop$cov1, arg2=pop$cov2)
  if(!all(temp >=0 & temp <=1))stop("p2 must be between 0 and 1. Your range is ", range(temp))
  pop$p2 <- temp

  # generate p(single tag)
  temp <- check.func(pST, "pST", 2, arg1=pop$cov1, arg2=pop$cov2)
  if(!all(temp >=0 & temp <=1))stop("pST must be between 0 and 1. Your range is ", range(temp))
  pop$pST <- temp

  # generate p(single tag) is in first position
  temp <- check.func(pST.1, "pST.1", 2, arg1=pop$cov1, arg2=pop$cov2)
  if(!all(temp >=0 & temp <=1))stop("pST.1 must be between 0 and 1. Your range is ", range(temp))
  pop$pST.1 <- temp

  # generate tag retention
  temp <- check.func(rho1, "rho1", 2, arg1=pop$cov1, arg2=pop$cov2)
  if(!all(temp >=0 & temp <=1))stop("rho1 must be between 0 and 1. Your range is ", range(temp))
  pop$rho1<- temp

  temp <- check.func(rho2, "rho2", 2, arg1=pop$cov1, arg2=pop$cov2)
  if(!all(temp >=0 & temp <=1))stop("rho2 must be between 0 and 1. Your range is ", range(temp))
  pop$rho2<- temp

  # find population that is tagged with one or two tags
  # event1, tag1, event1, tag2 etc
  # these are the actual status of the tags
  tag        <- runif(N) <= pop$p1
  single.tag <- runif(N) <= pop$pST
  single.tag.1 <- runif(N) <= pop$pST.1

  pop$e1t1 <- as.numeric( tag & (single.tag &  single.tag.1  | !single.tag))
  pop$e1t2 <- as.numeric( tag & (single.tag & !single.tag.1  | !single.tag))

  # are tags retained to event2
  pop$e2t1 <- pop$e1t1 * as.numeric( runif(N) <= pop$rho1)
  pop$e2t2 <- pop$e1t2 * as.numeric( runif(N) <= pop$rho2)

  # animals recaptured at t2 regardless of tag status
  recap <- runif(N) <= pop$p2

  # if not recaptured, then both tags not observed.
  pop$e2t1[ !recap] <- 0  # tag not observed
  pop$e2t2[ !recap] <- 0  # tag not observed

  # if not tagged at event 1 and recaptured, then set tag1 to 1
  pop$e2t1[ recap & pop$e1t1==0 & pop$e1t2==0] <- 1

  # if previously tagged, and all tags lost, and recaptured, then create a "recycled" fish
  # that looks like a 0010
  select <- (pop$e1t1 !=0 | pop$e1t2 !=0) & (pop$e2t1==0 & pop$e2t2==0) & recap
  new.pop <- pop[select,]
  if(sum(select)>0){  # if any recycled fish, add them to end of population
     new.pop$e1t1 <- 0
     new.pop$e1t2 <- 0
     new.pop$e2t1 <- 1
     new.pop$e2t2 <- 0
     pop <- rbind(pop, new.pop)
  }

  # create capture history
  pop$cap_hist=apply(pop[,c("e1t1","e1t2","e2t1","e2t2")], 1, FUN=paste0, collapse="")

  # change history for type of double tag
  if(dt_type == valid_dt_type()[1]){ # two indistinguishable double tags
     pop$cap_hist[ pop$cap_hist %in% c("1110","1101")] <- "111X"
  }
  if(dt_type == valid_dt_type()[3]){ # permanent second tag
     # this does not check that permanent tag was lost. User is responsible for checking this.
     select <- substr(pop$cap_hist,2,2) == "1"
     substr(pop$cap_hist[ select],2,2) <- "P"
     select <- substr(pop$cap_hist,4,4) == "1"
     substr(pop$cap_hist[ select],4,4) <- "P"
  }

  # remove animals never seen
  pop <- pop[ pop$cap_hist != '0000',]

  # collapse to capture histories
  data <- plyr::ddply(pop, c("cov1","cov2","cap_hist"), plyr::summarize,
                freq=length(cov1))
  data

}
