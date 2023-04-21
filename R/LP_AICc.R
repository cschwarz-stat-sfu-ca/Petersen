#' Create an AIC table comparing multiple LP fits
#'
#' This will take a series of LP fits and computes the usual AICc table and model weights
#'

#' @template param.fit
#'
#' @return An data frame with an AICc table and model weights etc
#'
#' @importFrom plyr ldply l_ply
#' @importFrom AICcmodavg AICcCustom

#' @examples
#'
#' data(data_rodli)
#' mt <- Petersen::LP_fit(data=data_rodli, p_model=~..time)
#' m0 <- Petersen::LP_fit(data=data_rodli, p_model=~1)
#' Petersen::LP_AICc(m0,mt)
#'
#' @export LP_AICc
#'

LP_AICc <- function(...){
  # Take a series of LP_fits
  # get all of the fits
  fits <-  c(as.list(environment()), list(...))

  # check that these are all LP_fit's.
  valid_models <- c("LP_fit",make_class_LP_TL(valid_dt_type()),"LP_SPAS_fit","LP_IS_fit")
  plyr::l_ply(fits, function(x){
     if( !inherits(x , valid_models))stop("Not all arguments in call are one of ", paste(valid_models, collapse=", "))
  })

  # check that all of the model are the same type, i.e. don't mix models
  first.class <-class(fits[[1]])
  #browser()
  plyr::l_ply(fits, function(x){
     if( !inherits(x , first.class))stop("Not all model of the same class  ")
  })

  #browser()
  # extract the summary from each of the fits
  aic.table <- plyr::ldply(fits, function(x){
    res <- x$summary
    res
  })

  #browser()
  aic.table$AICc  <- AICcCustom(logL =aic.table$cond.ll,
                                K    =aic.table$n.parms,
                                nobs =aic.table$nobs)
  aic.table       <- aic.table[ order(aic.table$AICc),]
  aic.table$Delta <- aic.table$AICc - min(aic.table$AICc)
  aic.table$AICcWt <- exp(-aic.table$Delta/2)
  aic.table$AICcWt <- aic.table$AICcWt / sum(aic.table$AICcWt)

  aic.table <- aic.table[,c("name_model","cond.ll","n.parms","nobs",'AICc',"Delta","AICcWt")]

  aic.table

}

#LP_AICc(rodli.fit.m0, rodli.fit.mt)
