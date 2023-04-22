#' Test for equal recapture probability in LP experiment
#'
#' This function takes the capture histories and stratification variable and computes the test
#' for equal recapture probabilities.
#'

#' @template param.data
#' @template param.strat_var
#' @param do.fisher.test Should a fisher exact test be done?
#' @importFrom stats na.pass xtabs

#' @return List containing the contingency table and the chi-square test and fisher-exact test
#'

#' @examples
#'
#' data(data_NorthernPike)
#' LP_test_equal_recap(data_NorthernPike, "Sex")
#'
#' @export LP_test_equal_recap
#'

LP_test_equal_recap <- function(data, strat_var, do.fisher.test=FALSE){
  # Contingency table test for equal recapture
  # check the data frame
  check.cap_hist.df(data)

  # check the do.fisher.tes
  if(!is.logical(do.fisher.test))stop("do.fisher.test must be logical of length 1")
  if(length(do.fisher.test)!=1)  stop("do.fisher.test must be logical of length 1")

  # check the stratification variable
  if(!is.character(strat_var))stop("Stratification variable must be character string")
  if(length(strat_var)!=1)    stop("Stratification variable must be length 1")
  if(!strat_var %in% names(data))stop("Stratification variable not in data frame")

  temp <- data[ data$cap_hist %in% c("11","10"),]
  temp$..strat_var <- temp[,strat_var]
  tab      <- xtabs(freq~cap_hist + ..strat_var, data=temp, exclude=NULL, na.action=na.pass)
  tab.prop <- prop.table(tab, margin=2)
  chisq.test  <- chisq.test(tab)
  if(do.fisher.test)fisher.test <- fisher.test(tab, simulate.p.value=TRUE)
  #browser()
  res <- list(examine ='equal recap',
       strat_var=strat_var,
       table    =tab,
       table.prop=tab.prop,
       chisq.test=chisq.test)
  if(do.fisher.test)res <- c(res, list(fisher.test=fisher.test))
  #browser()
  res
}

#data(NorthernPike)
#LP_test_equal_recap(NorthernPike, "Sex")
