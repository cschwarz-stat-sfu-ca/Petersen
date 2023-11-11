#' Extract estimates of abundance after BTSPAS fit
#'
#' This will take a previous fit and return estimates of abundance.
#'

#' @param LP_BTSPAS_fit A result of an call to fitting at BTSPAS object.
#' @template param.conf_level
#' @param trace If trace flag is set in call when estimating functions
#' @param parm Which parameter from the BTSPAS fix is to be extracted?
#'
#' @examples
#' \donttest{
#' # NOTE. To keep execution time to a small value as required by CRAN
#' # I've made a very small example.
#' # Additionally, I've set the number of MCMC chains, iterations, burning, simulation to save to
#' # small values. Proper mixing may not have occurred yet.
#' # When using this routine, you likely want to the use the default values
#' # for these MCMC parameters.
#'
#' data(data_btspas_diag1)

#' # extract the strata of interest
#' temp<- cbind(data_btspas_diag1,
#'              split_cap_hist( data_btspas_diag1$cap_hist,
#'                              sep="..", make.numeric=TRUE))

#' # only use data up to week 10 to keep example small
#' temp <- temp[ temp$t1 %in% 0:10 & temp$t2 %in% 0:10,]
#'
#' fit <- Petersen::LP_BTSPAS_fit_Diag(
#'   temp,
#'   p_model=~1,
#'   InitialSeed=23943242,
#'   # the number of chains and iterations are too small to be useful
#'   # they are set to a small number to pare execution time to <5 seconds for an example
#'   n.chains=2, n.iter=20000, n.burnin=1000, n.sims=100,
#'   quietly=TRUE
#' )
#' fit$summary
#'
#' # now get the estimates of abundance
#' est <-  Petersen::LP_BTSPAS_est (fit)
#' est$summary
#' }
#'
#'
#' @returns An list object of class *LP_BTSPAS_est* with the following elements
#' * **summary** A data frame  with the estimates of abundance, SE, and CI
#' * **datetime** Date and time the fit was done

#' @template author
#'
#' @importFrom formula.tools is.one.sided
#' @importFrom plyr is.formula
#' @importFrom stats as.formula model.matrix coef
#' @importFrom bbmle mle2

#' @export LP_BTSPAS_est
#'

LP_BTSPAS_est <- function(LP_BTSPAS_fit, parm="Ntot", conf_level=0.95, trace=FALSE){
  # After the Fit of model for the capture probability using BTSPAS

  # check the fitted objects
  if(!inherits(LP_BTSPAS_fit, c("LP_BTSPAS_fit_Diag","LP_BTSPAS_fit_NonDiag")))
    stop("LP_BTSPAS argument must be the results of a call to fitting BTSPAS objectl")

  # check the confidence level
  check.conf_level(conf_level)

  post <- extract_posterior(LP_BTSPAS_fit$fit, parm)  # allow to extract other parameters as well

  if(trace)browser()
  summary <-
       data.frame(N_hat_f = NA,
                  N_hat_rn= NA,
                  N_hat   = mean(post$value),
                  N_hat_SE= sd  (post$value),
                  N_hat_conf_level =conf_level,
                  N_hat_conf_method="Posterior",
                  N_hat_LCL = quantile(post$value, prob=  (1-conf_level)/2),
                  N_hat_UCL = quantile(post$value, prob=1-(1-conf_level)/2)
       )

  #browser()
  summary$p_model    = toString(LP_BTSPAS_fit$p_model)
  summary$name_model = paste0("p: ", toString(LP_BTSPAS_fit$p_model), collaspe="")
  summary$cond.ll    = LP_BTSPAS_fit$summary$cond.ll
  summary$n.parms    = LP_BTSPAS_fit$summary$n.parms
  summary$nobs       = sum(LP_BTSPAS_fit$data$freq)
  summary$method     = LP_BTSPAS_fit$method

  rownames(summary) <- NULL
  res <- list(summary=summary,
              datetime=Sys.time()
              )
  class(res)<- "LP_BTSPAS_est"
  res

}
