
#' Wrapper (*_fit) to call the Time Stratified Petersen Estimator
#' with Diagonal Entries function in BTSPAS.
#'
#' Takes the data structure as described below, and uses Bayesian methods to fit a fit a spline
#' through the population numbers and a hierarchical model for the trap
#' efficiencies over time.  An MCMC object
#' is also created with samples from the posterior.
#'
#' Use the \code{Petersen::LP_BTSPAS_NonDiag_fit} function for cases
#' where recaptures take place outside the stratum of release.
#'
#' @template param.data
#' @template param.p_model
#' @param p_model_cov Data frame with covariates for the model for prob capture at second sampling event. If this
#' data frame is given, it requires one line for each of the temporal strata at the second sampling event (even
#' if missing in the \code{data} that has the capture histories) with one variable being \code{..time}
#'  to represent
#' the second temporal stratum.
#' @template param.BTSPAS.jump.after
#' @template param.logitP.fixed
#' @param InitialSeed Numeric value used to initialize the random numbers used
#' in the MCMC iterations.
#' @param trace  Internal tracing flag.


#' @details
#'
#' The frequency variable (\code{freq} in the \code{data} argument) is the number of animals with the corresponding capture history.
#'
#' Capture histories (\code{cap_hist} in the \code{data} argument) are character values of the format
#' \code{xx..yy} is a capture_history where \code{xx} and \code{yy} are the temporal stratum
#' (e.g., julian week) and \code{'..'} separates
#' the two temporal strata.
#'  If a fish is released in temporal stratum and never captured again, then \code{yy} is set to 0;
#'  if a fish is newly captured in temporal stratum \code{yy}, then \code{xx} is set to zero.
#'  For example, a capture history of \code{23..23} indicates animals released in temporal stratum
#'  23 and recaptured in temporal stratum 23; a capture history of \code{23..00}
#'   indicates animals released in temporal stratum
#'  23 and never seen again; a capture history of \code{00..23}
#'   indicates animals newly captured in temporal stratum
#'  23 at the second sampling event.
#'
#'  In the diagonal case, no fish should move between temporal strata.
#'
#'  It is not necessary to label the temporal strata starting at 1; BTSPAS will treat the smallest
#'  value of the temporal strata seen as the first stratum and will interpolate for temporal strata
#'  without any data. Temporal strata labels should be numeric, i.e., do NOT use A, B, C etc.

#'
#
#'
#' @return An MCMC object with samples from the posterior distribution. A
#' series of graphs and text file are also created in the results along with
#' a summary object

#' @importFrom stats runif var sd
#' @importFrom BTSPAS TimeStratPetersenDiagError_fit
#'
#' @export LP_BTSPAS_Diag_fit
#'
#'

LP_BTSPAS_Diag_fit <- function(
     data,
     p_model=~1,
     p_model_cov=NULL,
     jump.after=NULL,
     logitP.fixed=NULL, logitP.fixed.values=NULL,
     InitialSeed=ceiling(stats::runif(1,min=0, max=1000000)),
     trace=FALSE){

  # some basic data checking
  check.cap_hist_temporal.df(data)

  # check that p_model_cov is a data frame
  if(!is.null(p_model_cov) && !is.data.frame(p_model_cov))stop("p_model_cov must be a data frame")

  # check the model for p. Must be a formula and all the variables must be in data frame p_model_cov
  if(!plyr::is.formula(p_model))stop("p_model must be a formula")
  if(!formula.tools::is.one.sided(p_model))stop("p_model must be one sided formula")
  temp <- as.character(p_model)
  if(length(temp)>1)stop("p_model must be length 1")
  p_model_vars <-all.vars(p_model)
  if(length(p_model_vars)>0){
     p_model_vars <- c(p_model_vars)
     temp <- p_model_vars %in% c(names(p_model_cov))
     if(!all(temp))stop("p_model refers to variables not in data :",
                     paste(p_model_vars[!temp],sep=", ", collapse=""))
  }

  nmu <- cap_hist_to_n_m_u(data)
  if(trace)browser()
  # get the model matrix for p
  if(length(all.vars(p_model))==0) logitP.cov <- rep(1, max(nmu$..ts, na.rm=TRUE)-min(nmu$..ts, na.rm=TRUE)+1)
  if(length(all.vars(p_model)) >0) logitP.cov <- model.matrix(p_model, data=p_model_cov)

  # now make the call
  res <- BTSPAS::TimeStratPetersenDiagError_fit(
           title="", prefix="TSPDE-",
           time      = nmu$..ts,
           n1        = nmu$n1,
           m2        = as.vector(nmu$m2),
           u2        = nmu$u2,
           jump.after= jump.after,
           bad.n1=c(), bad.m2=c(), bad.u2=c(),
           logitP.cov= logitP.cov,
           logitP.fixed=logitP.fixed, logitP.fixed.values=logitP.fixed.values,
           n.chains=3, n.iter=200000, n.burnin=100000, n.sims=2000,
           tauU.alpha=1, tauU.beta=.05, taueU.alpha=1, taueU.beta=.05,
           prior.beta.logitP.mean = c(logit(sum(nmu$m2,na.rm=TRUE)/sum(nmu$n1,na.rm=TRUE)),
                                    rep(0,  ncol(as.matrix(logitP.cov))-1)),
           prior.beta.logitP.sd   = c(stats::sd(logit((nmu$m2+.5)/(nmu$n1+1)),na.rm=TRUE),
                                    rep(10, ncol(as.matrix(logitP.cov))-1)),
           tauP.alpha=.001, tauP.beta=.001,
           run.prob=seq(0,1,.1),  # what percentiles of run timing are wanted
           debug=FALSE, debug2=FALSE,
           InitialSeed=InitialSeed,
           save.output.to.files=FALSE,
           trunc.logitP=15
  )

  summary <- data.frame(
     p_model    = paste0(as.character(p_model), collapse=""),
     name_model = paste0("p: ", toString(p_model), collapse=""),
     cond.ll    = NA,
     n.parms    = NA,
     nobs       = sum(data$freq),
     method     = "BTSPAS-Diag"
  )

  res <- list(summary    = summary,
              data       = data,
              p_model    = p_model,
              p_model_cov= p_model_cov,
              jump.after = jump.after,
              InitialSeed= InitialSeed,
              name_model=paste("p: ", toString(p_model), collapse=""),
              fit=res,
              datetime=Sys.time())
  class(res) <- "BTSPAS-Diag-fit"

  res
}

