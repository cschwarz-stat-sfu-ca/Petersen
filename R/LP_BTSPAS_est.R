#' Extract estimates of abundance after BTSPAS fit
#'
#' This will take a previous fit and return estimates of abundance.
#'

#' @param LP_BTSPAS_fit A result of an call to fitting at BTSPAS object.
#' @template param.conf_level
#' @param trace If trace flag is set in call when estimating functions
#' @param parm Which parameter from the BTSPAS fix is to be extracted?
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

#' @export LP_BTSPAS_est
#'

LP_BTSPAS_est <- function(LP_BTSPAS_fit, parm="Ntot", conf_level=0.95, trace=FALSE){
  # After the Fit of model for the capture probability using BTSPAS

  # check the fitted objects
  if(!inherits(LP_BTSPAS_fit, c("BTSPAS-Diag-fit","BTSPAS-NonDiag-fit")))
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
