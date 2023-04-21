#' Extract estimates of abundance after SPAS fit
#'
#' This will take a previous fit and return estimates of abundance.
#'

#' @param LP_SPAS_fit A result of an call to fitting at SPAS object.
#' @template param.conf_level
#' @param trace If trace flag is set in call when estimating functions
#'
#' @return An list object with abundance estimates
#' @template author
#'
#' @importFrom formula.tools is.one.sided
#' @importFrom plyr is.formula
#' @importFrom stats as.formula model.matrix coef
#' @importFrom bbmle mle2

#' @examples
#'
#' #data(data_rodli)
#' #rodli.fit <- Petersen::LP_fit(data=data_rodli, p_model=~..time)
#' #Petersen::LP_est(rodli.fit, N_hat=~1)

#' @export LP_SPAS_est
#'

LP_SPAS_est <- function(LP_SPAS_fit, conf_level=0.95, trace=FALSE){
  # After the Fit of model for the capture probability using SPAS

  # check the fitted objects
  if(!inherits(LP_SPAS_fit, c("LP_SPAS_fit")))
    stop("LP_SPAS argument must be the results of a call to fitting SPAS objectl")

  # check the confidence level
  check.conf_level(conf_level)


  if(trace)browser()
  summary <-
       data.frame(N_hat_f = NA,
                  N_hat_rn= NA,
                  N_hat   = LP_SPAS_fit$fit$est$real$N,
                  N_hat_SE= LP_SPAS_fit$fit$se $real$N,
                  N_hat_conf_level =conf_level,
                  N_hat_conf_method="Large sample",
                  N_hat_LCL = LP_SPAS_fit$fit$est$real$N - qnorm(1-(1-conf_level)/2)*LP_SPAS_fit$fit$se $real$N,
                  N_hat_UCL = LP_SPAS_fit$fit$est$real$N + qnorm(1-(1-conf_level)/2)*LP_SPAS_fit$fit$se $real$N
       )


  #browser()
  summary$p_model    = LP_SPAS_fit$summary$p_model
  summary$name_model = LP_SPAS_fit$summary$name_model
  summary$cond.ll    = LP_SPAS_fit$summary$cond.ll
  summary$n.parms    = LP_SPAS_fit$summary$n.parms
  summary$nobs       = LP_SPAS_fit$summary$nobs
  summary$method     = LP_SPAS_fit$summary$method
  summary$cond.factor= LP_SPAS_fit$summary$cond.factor

  rownames(summary) <- NULL
  res <- list(summary=summary,
              datetime=Sys.time()
              )
  class(res) <- "LP_SPAS_est"
  res

}
