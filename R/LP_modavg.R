#' Create an table of individual estimates and the model averaged values
#'
#' This will take a series of LP fits and computes the model averages for each set of N_hat
#'

#' @template param.fit
#' @template param.N_hat
#' @template param.conf_level

#' @return An data frame with model averaged values for abundance
#'
#' @importFrom plyr ldply l_ply
#' @importFrom AICcmodavg modavgCustom

#' @examples
#'
#' data(data_rodli)
#' mt <- Petersen::LP_fit(data=data_rodli, p_model=~..time)
#' m0 <- Petersen::LP_fit(data=data_rodli, p_model=~1)
#' Petersen::LP_modavg(m0,mt)
#'
#' @export LP_modavg
#'

LP_modavg <- function(..., N_hat=~1,conf_level=0.95){
  # Take a series of LP_fits

  # get all of the fits
  fits <-  c(as.list(environment()), list(...))

  # exclude the N_hat
  fits <- fits[ !names(fits) %in% c("N_hat","conf_level")]
  #browser()

  valid_models <- c("LP_fit",make_class_LP_TL(valid_dt_type()),"LP_SPAS_fit","LP_IS_fit")
  plyr::l_ply(fits, function(x){
    if(!inherits(x,valid_models))stop("All objects must be one of ", paste(valid_models, collapse=", "))
  })

  # check that all of the model are the same type, i.e. don't mix models
  first.class <-class(fits[[1]])
  plyr::l_ply(fits, function(x){
     if( !inherits(x , first.class))stop("Not all model of the same class  ")
  })

  if(!is.formula(N_hat))stop("N_hat argument must be a formula. You have ", N_hat)
  vars <- all.vars(N_hat) # extract the variables and allow access to the Expansion factor as well
  if(!all(vars %in% c(names(fits[[1]]$data),"..EF","..cat")))stop("Invalid variables in formula for N_hat :",vars)

  # check the confidence level
  check.conf_level(conf_level)

  # get the aic table
  aic.table <- LP_AICc(...)

  # get a table of estimates.
  est <- plyr::ldply(fits, function(x){
     if(inherits(x,"LP_fit")){
       est <- LP_est(x, N_hat=N_hat, conf_level=conf_level)
     }
     if(inherits(x,make_class_LP_TL(valid_dt_type()))){
       est <- LP_TL_est(x, N_hat=N_hat, conf_level=conf_level)
     }
     if(inherits(x,"LP_SPAS_fit")){
       est <- LP_SPAS_est(x)
     }
     if(inherits(x,"LP_IS_fit")){
       est <- LP_IS_est(x, N_hat=N_hat, conf_level=conf_level)
     }
     est$summary
  })
  #browser()

  # append the aic stufff
  # get the model averages for each type of estimate
  ma.est <- plyr::ddply(est, c(".id","N_hat_f","N_hat_rn"), function(x){
      ma <- AICcmodavg::modavgCustom( logL=x$cond.ll, K=x$n.parms, nobs=x$nobs, modnames=x$name_model,
                                      estimate=x$N_hat, se=x$N_hat_SE, conf.level=conf_level)
      #browser()
      res <- ma$Mod.avg.table
      #browser()

      res <- merge(res, x[,c("name_model","N_hat_f","N_hat_rn","N_hat_conf_level","N_hat_conf_method","N_hat_LCL","N_hat_UCL")],
                   by.x="Modnames", by.y="name_model")
      res <- res[ order(res$AICcWt, decreasing=TRUE),]
      res <- plyr::rbind.fill(res,
                              data.frame(
                                Modnames="Model averaged",
                                Estimate=ma$Mod.avg.est,
                                SE      =ma$Uncond.SE,
                                N_hat_conf_level=ma$Conf.level,
                                N_hat_LCL     =ma$Lower.CL,
                                N_hat_UCL     =ma$Upper.CL))
      res
  })
  ma.est
}

#LP_AICc(rodli.fit.m0, rodli.fit.mt)
