#' Estimate abundance after the LP_IS conditional likelihood fit.
#'
#' This will take a previous fit and return estimates of abundance.
#' The population abundance is estimated using
#' a Horvitz-Thompson type estimator and the user can request abundance
#' estimates for sub-sets of the population
#'

#' @param LP_IS_fit A result of an LP_IS_fit() call.
#' @template param.N_hat
#' @template param.conf_level
#' @param trace If trace flag is set in call when estimating functions
#'
#' @returns An list object with abundance estimates and other information with the following elements
#'  * **summary** Data frame with abundance estimates, their SE, and CIs as requested
#'  * **detail** List with many components, including the rawdata, model fitting information, observed and expected values, residual plot, etc
#' * **datetime** Date and time the estimation was done from the fit.
#'

#' @template author
#'
#' @importFrom formula.tools is.one.sided
#' @importFrom plyr is.formula
#' @importFrom stats as.formula model.matrix coef

#' @examples
#'
#' data(data_wae_is_short)
#' fit <- Petersen::LP_IS_fit(data=data_wae_is_short, p_model=~..time)
#' fit$summary
#' est <- LP_IS_est(fit, N_hat=~1)
#' est$summary

#' @export LP_IS_est
#'

LP_IS_est <- function(LP_IS_fit, N_hat=~1, conf_level=0.95, trace=FALSE){
  # After the incomplete stratification Fit, extract the estimates
  N_hat_f = N_hat # save the formula for N_hat

  # check that LP_IS_fit is a fitted objects
  if(!inherits(LP_IS_fit, "LP_IS_fit"))stop("LP_IS_fit argument must be the results of an LP_IS_fit call")

  # check the formula for N_hat
  # we currently only allow ~1 and ~..cat in formula for N_ht
  if(!plyr::is.formula(N_hat))stop("N_hat argument must be a formula. You have ", N_hat)
  N_hat.vars <- all.vars(N_hat) # extract the variables and allow access to the Expansion factor as well
  if(length(N_hat.vars)>0){
    if(!all(N_hat.vars %in% c("..cat")))stop("Invalid variables in formula for N_hat. Only ~1 and ..cat allowed")
  }

  # check the confidence level
  check.conf_level(conf_level)

  est <- LP_IS_fit$fit$est
  se  <- LP_IS_fit$fit$se

  #browser()

  #extract the N_hats if there are more than 1 and their respective SE
  summary <-
       data.frame(N_hat_f = toString(N_hat),
                  N_hat_conf_level =conf_level,
                  N_hat_conf_method="logN")
  if(length(N_hat.vars)==0){  # estimate of total abundance
     summary$N_hat_rn <- "All"
     summary$N_hat    <- est$N
     summary$N_hat_SE <- se $N
  }
  if(length(N_hat.vars)>0){  # estimate of individual stratum abundances
     summary <- summary[rep(1,length(rownames(est$p))),]
     summary$N_hat_rn <- rownames(est$p)
     summary$N_hat    <- est$N_lambda
     summary$N_hat_SE <- se $N_lambda
  }

  summary$N_hat_LCL = exp(log(summary$N_hat) - qnorm(1-(1-conf_level)/2)*summary$N_hat_SE/summary$N_hat)
  summary$N_hat_UCL = exp(log(summary$N_hat) + qnorm(1-(1-conf_level)/2)*summary$N_hat_SE/summary$N_hat)

  #browser()
  summary$p_model    = LP_IS_fit$summary$p_model
  summary$theta_model= LP_IS_fit$summary$theta_model
  summary$lambda_mode= LP_IS_fit$summary$lambda_model
  summary$name_model = LP_IS_fit$summary$name_model
  summary$cond.ll    = LP_IS_fit$summary$cond.ll
  summary$n.parms    = LP_IS_fit$summary$n.parms
  summary$nobs       = LP_IS_fit$summary$nobs
  summary$method     = LP_IS_fit$summary$method

  res <- list(summary=summary,
              detail =LP_IS_fit$fit,
              datetime=Sys.time()
              )

  res

}
