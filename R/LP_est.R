#' Estimate abundance after the LP conditional likelihood fit.
#'
#' This will take a previous fit and return estimates of abundance.
#' The population abundance is estimated using
#' a Horvitz-Thompson type estimator and the user can request abundance
#' estimates for sub-sets of the population.
#'

#' @param LP_fit A result of an LP_fit() call.
#' @template param.N_hat
#' @template param.conf_level
#' @param trace If trace flag is set in call when estimating functions
#'
#'
#' @returns An list object with abundance estimates and other information with the following elements
#' * **summary** Data frame with abundance estimates, their SE, and CIs as requested
#' * **detail** List with many components, including the rawdata, model fitting information, observed and expected values, residual plot, etc
#' * **datetime** Date and time the estimation was done from the fit.

#' @template author
#'
#' @importFrom formula.tools is.one.sided
#' @importFrom plyr is.formula
#' @importFrom stats as.formula model.matrix coef
#' @importFrom bbmle mle2

#' @examples
#'
#' data(data_rodli)
#' fit <- Petersen::LP_fit(data=data_rodli, p_model=~..time)
#' fit$summary
#' est <- Petersen::LP_est(fit, N_hat=~1)
#' est$summary

#' @export LP_est
#'

LP_est <- function(LP_fit, N_hat=~1, conf_level=0.95, trace=FALSE){
  # After the Fit of model for the capture probability using conditional likelihood
  # and use a Horvitz-Thompson estimator to estimate abundance
  N_hat_f = N_hat # save the formula for N_hat

  # check that LP_fit is a fitted objects
  if(!inherits(LP_fit, "LP_fit"))stop("LP_fit argument must be the results of an LP_fit call")

  # check the formula for N_hat
  if(!is.formula(N_hat))stop("N_hat argument must be a formula. You have ", N_hat)
  vars <- all.vars(N_hat) # extract the variables and allow access to the Expansion factor as well
  if(!all(vars %in% c(names(LP_fit$data),"..EF")))stop("Invalid variables in formula for N_hat")

  # check the confidence level
  check.conf_level(conf_level)

  # now to find the estimates of abundance est
  est <- LP_cond_lik(LP_fit$fit@coef,
                     data       =LP_fit$data,
                     p_model    =LP_fit$p_model,
                     p_beta_vcov=LP_fit$fit@vcov,
                     what.return="estimates",
                     N_hat = N_hat,
                     conf_level=conf_level,
                     trace=trace)
  #extract the N_hats if there are more than 1 and their respective SE
  summary <-
       data.frame(N_hat_f = paste0(as.character(N_hat_f),collapse=""),
                  N_hat_rn= row.names(est$N_hat$N_hat),
                  N_hat   = as.vector(est$N_hat$N_hat),
                  N_hat_SE= as.vector(est$N_hat$N_hat_SE),
                  N_hat_conf_level =as.vector(est$N_hat$N_hat_conf_level),
                  N_hat_conf_method=as.vector(est$N_hat$N_hat_conf_method),
                  N_hat_LCL = as.vector(est$N_hat$N_hat_LCL),
                  N_hat_UCL = as.vector(est$N_hat$N_hat_UCL)
       )

  #browser()
  summary$p_model    = toString(est$p_model)
  summary$name_model = paste0("p: ", toString(est$p_model), collaspe="")
  summary$cond.ll    = est$cond.ll
  summary$n.parms    = est$n.parms
  summary$nobs       = sum(est$data$freq)
  summary$method     = est$method

  est$fit = LP_fit

  res <- list(summary=summary,
              detail =est,
              datetime=Sys.time()
              )

  res

}
