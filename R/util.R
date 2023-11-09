# Utility functions used in this package.
# These functions will not be visible to the user.


#' @param variable Variable to be checked if numeric
#' @param min.value Minimum value for variable
#' @param max.value Maximum value for variable
#' @param req.length Required length for variable
#' @param check.whole Check if all values in variable are whole numbers
#' @param data A  data frame
#' @param dt_type Type of double tags. See valid_dt_type()
#' @param sep. Separator for capture histories

#'
#' @importFrom reshape2 melt
#' @noRd



#######################################################################3
# Check that argument is numeric, specified length, and all values are between min /max
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol}  # taken from examples of is.integer


#######################################################################3

check.numeric <- function(variable, min.value=0, max.value=Inf, req.length=1, check.whole=TRUE){
  var.name <- deparse(substitute(variable))
  if(!is.numeric(min.value) | !is.numeric(max.value) | !is.numeric(req.length))
    stop("Bad input values to check.numeric: ", paste(c(var.name,"; ", min.value,"; ", max.value,"; ", req.length), collapse=", "))
  if(length(req.length) !=1)stop("Bad required length for check.numeric. You had ", req.length)
  if(length(min.value) !=1 & length(min.value) != req.length)
    stop("Min values wrong length in check.numeric. You have: ",paste(min.value, collapse=", "))
  if(length(max.value) !=1 & length(max.value) != req.length)
    stop("Max values wrong length in check.numeric. You have: ",paste(max.value, collapse=", "))
  if(!is.numeric(variable))stop("'",var.name, "' must be numeric")
  if( any(is.na(variable)))stop("'",var.name, "' cannot be missing. You have ", paste(variable, collapse=", "))
  if( length(variable) != req.length)stop("'",var.name, "' must have length ", req.length)
  if( any(variable < min.value) | any(variable > max.value))
    stop("'",var.name, "' must between (",paste(min.value, collapse=", "),") and (",
         paste(max.value, collapse=", "),
         ') (inclusive). You have ',paste(variable, collapse=", "))
  # check if values are all whole numbers
  if(check.whole){
    if(any(!is.wholenumber(variable)))stop("'",var.name, "' must be integers. You had ", paste(variable, collapse=", "))
  }
  invisible()
}


#######################################################################3
# check the capture history data frame for validity
check.cap_hist.df <- function(data, type=c("LP","LP_TL")[1], dt_type=""){
  if(!is.character(type))stop("Type of data set must be character")
  if(!all(type %in% c("LP","LP_TL")))stop("data type must be LP or LP_TL")

  if(!is.data.frame(data))stop("Capture history data must be a data frame")
  if(!all(c("cap_hist","freq")%in% names(data)))stop("Capture history data frame must contain a 'cap_hist' and 'freq' variable")
  if("..time" %in% names(data))warning("*** Capture-history data frame contains a variable '..time'  that may be overwritten ")
  if("..tag"  %in% names(data))warning("*** Capture-history data frame contains a variable '..tag'   that may be overwritten ")
  if("..EF"   %in% names(data))warning("*** Capture-history data frame contains a variable '..EF'    that may be overwritten ")

  if(any(is.na(data$freq)))stop("All frequencies must be present and non-negative")
  check.numeric(data$freq, min.value=0, req.length=nrow(data), check.whole=FALSE)

  if(any(is.na(data$cap_hist)))stop("All capture-histories 'cap_hist' must be non-missing")
  if(!is.character(data$cap_hist))stop("data$cap_hist must be character vector")

  if(!is.character(dt_type) | length(dt_type)!=1)stop("dt_type set must be character vector of length 1")
  if(type=="LP_TL" & dt_type=="")stop("Specify dt_type when type=LP_TL")
  if(type=="LP_TL" & !dt_type %in% valid_dt_type())stop("Invalid dt_type ", dt_type)

  # check for valid CH values
  if(type=="LP"  )  valid.cap_hist <- c("10","01","11")
  if(type=="LP_TL" & dt_type ==valid_dt_type()[1] )  valid.cap_hist <- c("1000","1010","1100","111X","1111","0010")
  if(type=="LP_TL" & dt_type ==valid_dt_type()[2] )  valid.cap_hist <- c("1000","1010","0100","0101","1100","1101","1110","1111","0010")
  if(type=="LP_TL" & dt_type ==valid_dt_type()[3] )  valid.cap_hist <- c("1000","1010","0P00","0P0P","1P00","1P0P",       "1P1P","0010")
  if(!all(data$cap_hist %in% valid.cap_hist))stop("All capture.histories must be one of ",paste(valid.cap_hist, sep="", collapse=", "))

  invisible()
}

check.cap_hist_temporal.df <- function(data, sep=".."){
  if(!is.data.frame(data))stop("Capture history data must be a data frame")
  if(!all(c("cap_hist","freq")%in% names(data)))stop("Capture history data frame must contain a 'cap_hist' and 'freq' variable")
  if("..time" %in% names(data))warning("*** Capture-history data frame contains a variable '..time'  that may be overwritten ")
  if("..tag"  %in% names(data))warning("*** Capture-history data frame contains a variable '..tag'   that may be overwritten ")
  if("..EF"   %in% names(data))warning("*** Capture-history data frame contains a variable '..EF'    that may be overwritten ")

  if(any(is.na(data$freq)))stop("All frequencies must be present and non-negative")
  check.numeric(data$freq, min.value=0, req.length=nrow(data), check.whole=FALSE)

  if(any(is.na(data$cap_hist)))stop("All capture-histories 'cap_hist' must be non-missing")
  if(!is.character(data$cap_hist))stop("data$cap_hist must be character vector")

  # now to try and split the capture histories into the temporal strata and see if ok
  strata <- tryCatch({split_cap_hist(data$cap_hist, sep=sep, make.numeric=TRUE, prefix="ts")
                     },
        error=function(cond) {
             stop("Unable to separate the capture histories. Likely invalid format??")
        },
        warning=function(cond) {
            stop(cond)
        }
  )
  # check if any missing values in temporal strata
  if(any(is.na(strata[,"ts1"])))stop("Some of the temporal strata were invalid ", paste(strata[,"ts1"],collapse=", "))
  if(any(is.na(strata[,"ts2"])))stop("Some of the temporal strata were invalid ", paste(strata[,"ts2"],collapse=", "))
  if(any(strata<0))stop("All temporal strata must be non-negative")

  invisible()
}

# check the data frame for geographic stratification
check.cap_hist_geographic.df <- function(data, sep=".."){
  if(!is.data.frame(data))stop("Capture history data must be a data frame")
  if(!all(c("cap_hist","freq")%in% names(data)))stop("Capture history data frame must contain a 'cap_hist' and 'freq' variable")
  if("..time" %in% names(data))warning("*** Capture-history data frame contains a variable '..time'  that may be overwritten ")
  if("..tag"  %in% names(data))warning("*** Capture-history data frame contains a variable '..tag'   that may be overwritten ")
  if("..EF"   %in% names(data))warning("*** Capture-history data frame contains a variable '..EF'    that may be overwritten ")

  if(any(is.na(data$freq)))stop("All frequencies must be present and non-negative")
  check.numeric(data$freq, min.value=0, req.length=nrow(data), check.whole=FALSE)

  if(any(is.na(data$cap_hist)))stop("All capture-histories 'cap_hist' must be non-missing")
  if(!is.character(data$cap_hist))stop("data$cap_hist must be character vector")

  # now to try and split the capture histories into the temporal strata and see if ok
  strata <- tryCatch({split_cap_hist(data$cap_hist, sep=sep, make.numeric=FALSE, prefix="gs")
  },
  error=function(cond) {
    stop("Unable to separate the capture histories. Likely invalid format??")
  },
  warning=function(cond) {
    stop(cond)
  }
  )
  # check if any missing values in geographic strata
  if(any(is.na(strata[,"gs1"])))stop("Some of the geographic strata were invalid ", paste(strata[,"gs1"],collapse=", "))
  if(any(is.na(strata[,"gs2"])))stop("Some of the geographic strata were invalid ", paste(strata[,"gs2"],collapse=", "))
  if(any(strata==""))stop("All geographic strata must be non-blank")

  invisible()
}


# check the capture history data frame for validity for incomplete stratification
check.cap_hist_IS.df <- function(data){

  if(!is.data.frame(data))stop("Capture history data must be a data frame")
  if(!all(c("cap_hist","freq")%in% names(data)))stop("Capture history data frame must contain a 'cap_hist' and 'freq' variable")
  if("..time" %in% names(data))warning("*** Capture-history data frame contains a variable '..time'  that may be overwritten ")
  if("..tag"  %in% names(data))warning("*** Capture-history data frame contains a variable '..tag'   that may be overwritten ")
  if("..EF"   %in% names(data))warning("*** Capture-history data frame contains a variable '..EF'    that may be overwritten ")
  if("..cat"  %in% names(data))warning("*** Capture-history data frame contains a variable '..cat'   that may be overwritten ")

  if(any(is.na(data$freq)))stop("All frequencies must be present and non-negative")
  check.numeric(data$freq, min.value=0, req.length=nrow(data), check.whole=FALSE)

  if(any(is.na(data$cap_hist)))stop("All capture-histories 'cap_hist' must be non-missing")
  if(!is.character(data$cap_hist))stop("data$cap_hist must be character vector")

  # check for valid CH values
  # first extract the categories that appear in the capture histories
  data <- cbind(data, split_cap_hist(data$cap_hist, sep=""))
  strata <- setdiff(unique(c(data$t1, data$t2)), "0") # ignore the 0 code
  valid.cap_hist <- c( paste0(strata,"0"), paste0("0",strata),paste0(strata,strata))

  if(!all(data$cap_hist %in% valid.cap_hist))stop("All capture.histories must be one of ",paste(valid.cap_hist, sep="", collapse=", "))

  invisible()
}


#######################################################################3
# check the confidence level
check.conf_level <- function(conf_level){
  if(!is.numeric(conf_level))stop("Confidence level must be numeric")
  if(length(conf_level) != 1)stop("Confidence level must be vector of length 1")
  if(conf_level <= 0 | conf_level >= 1)stop("Confidence level must be between 0 and 1")
  if(conf_level < .80)warning("Confidence level is less than .80. Are you sure?")

  invisible()

}


#######################################################################3
# what are valid dt_types (types of double tags)
# inorder, 2 tags, not distinguishable; two tags, distinguisable; second tag is permanent
valid_dt_type <- function(){c('notD','twoD','t2perm')}


#######################################################################3
# create class variables from LP_TL_fit and dt_type
make_class_LP_TL <- function(dt_type){paste0("LP_TL_fit","-",dt_type)}


#######################################################################3
extract_posterior <- function(fit, effect.name, source="Bayesian"){
  # extract the posterior for an effect in the long-format
  # index.type must be one of "year", "site", or "year.site"
  #browser()
  select <- grepl(effect.name, colnames(fit$sims.matrix) )
  post   <- as.data.frame(fit$sims.matrix[,select,drop=FALSE])

  post$sim <- 1:nrow(post)
  post.long <- reshape2::melt(post,
                              id.vars="sim",
                              value.name="value")
  #browser()
  post.long$Source <- source
  post.long$variable <- as.character(post.long$variable)
  post.long
}

#######################################################################3
# suppress cat() messages from function call
# see https://stackoverflow.com/questions/34208564/how-to-hide-or-disable-in-function-printed-message
quiet.eval <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

