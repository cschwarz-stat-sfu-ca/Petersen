#' Fit a Lincoln-Petersen Model with Tag Loss using conditional likelihood
#'
#' This will take a data frame of capture histories, frequencies, and
#' additional covariates (e.g., strata and/or continuous covariates) and the model
#' for p1 and the tag retention probabilities
#' and will use conditional likelihood (conditional on capture at time 2)
#' to fit the model. The population abundance is estimated using
#' a Horvitz-Thompson type estimator and the user can request abundance
#' estimates for sub-sets of the population.
#'

#' @template param.data
#' @template param.dt_type
#' @template param.p_model
#' @param rho_model Model for retention probabilities
#' @param all_beta.start Initial values for call to optimization routine for the beta parameters (on the logit scale).
#' The values will be replicated to match
#' the number of initial beta parameters needed. Some care is needed here since the parameter order are for the p1 probabilities
#' and then for the rho probabilities
#' @param trace If trace flag is set in call when estimating functions

#' @details
#'
#' The frequency variable (\code{freq} in the \code{data} argument) is the number of animals with the corresponding capture history.
#'
#' Capture histories (\code{cap_hist} in the \code{data} argument) are character values of length 4.
#'
#' If the tag loss model is two indistinguishable tags (\code{dt_type="notD"}), then valid capture histories are:
#' \itemize{
#'   \item \strong{1100}  Animals double tagged but never seen again.
#'   \item \strong{111X}  Animals double tagged, but only 1 tag was present when animal recaptured at second event.
#'   \item \strong{1111}  Animals double tagged and both tags present when animal recaptured at second event.
#'   \item \strong{1000}  Animals single tagged and never seen again.
#'   \item \strong{0100}  Animals single tagged and never seen again.
#'   \item \strong{1010}  Animals single tagged and recaptured with the single tag.
#'   \item \strong{0101}  Animals single tagged and recaptured with the single tag.
#'   \item \strong{0010}  Animals APPARENTLY captured for the first time at event 2. This includes animals that are
#'   newly captured, plus fish that were tagged and lost all their tags, and were captured again
#' }
#'
#' If the tag loss model is two distinguishable tags (\code{dt_type="twoD"}), then valid capture histories are the same
#' as above except the history \code{111X} is replaced by:
#' \itemize{
#'   \item \strong{1110} Animals double tagged, but only the first of the double tags applied  was present when animal recaptured at event 2,
#'   \item \strong{1101} Animals double tagged, but only the second of the double tags applied was present when animal recaptured at event 2.
#'  }
#'
#' If the second tag is a permanent batch mark (\code{dt_type="t2perm"}), then valid capture histories are:
#' \itemize{
#'   \item \strong{1P00}  Animals double tagged but never seen again.
#'   \item \strong{1P0P}  Animals double tagged,but non-permanent tag missing when animal recaptured at second event.
#'   \item \strong{1P1P}  Animals double tagged and both tags present when animal recaptured at second event.
#'   \item \strong{1000}  Animals single tagged and never seen again.
#'   \item \strong{0P00}  Animals single tagged with a permanent batch mark only and never seen again.
#'   \item \strong{1010}  Animals single tagged and recaptured with the single tag.
#'   \item \strong{0P0P}  Animals single tagged with the permanent batch mark and recaptured with the permanent tag.
#'   \item \strong{0010}  Animals APPARENTLY captured for the first time at event 2. This includes animals that are
#'   newly captured, plus fish that were tagged and lost all their tags, and were captured again
#' }
#'
#' @return An list object with the fit information
#' @template author
#'
#' @importFrom formula.tools is.one.sided
#' @importFrom plyr is.formula
#' @importFrom stats as.formula model.matrix coef
#' @importFrom bbmle mle2

#' @examples
#'
#' data(data_kokanee_tagloss)
#' Petersen::LP_TL_fit(data=data_kokanee_tagloss, p_model=~1, rho_model=~1, dt_type="notD")
#'
#' @export LP_TL_fit
#'

LP_TL_fit <- function(data,
                     dt_type=NULL,
                     p_model, rho_model, all_beta.start=NULL, trace=FALSE){
  # Fit a model for the Peterson model with tagloss using conditional likelihood


  # check that double tag type
  if(is.null(dt_type))stop("dt_type must be specified. Valid values are ", paste(valid_dt_type(),collape=", "))
  if(!is.character(dt_type) | length(dt_type) !=1)stop("dt_type must be character of length 1")
  valid.values <- valid_dt_type()
  if(!all(dt_type %in% valid.values))stop("dt_type must be one of ", paste(valid.values, collapse=", "))

  # check the data frame to ensure that it is valid for the data type
  check.cap_hist.df(data, type="LP_TL", dt_type=dt_type)

  # check the model for p1. Must be a formula and all the variables must be in data frame
  if(!plyr::is.formula(p_model))stop("p_model must be a formula")
  if(!formula.tools::is.one.sided(p_model))stop("p_model must be one sided formula")
  temp <- as.character(p_model)
  if(length(temp)>1)stop("p_model must be length 1")
  # we need to convert : to * to properly extract variables
  p_model_vars <- all.vars(p_model)
  temp <- p_model_vars %in% c(names(data))
  if(!all(temp))stop("p_model refers to variables not in data :",
                     paste(p_model_vars[!temp],sep=", ", collapse=""))

  if(!plyr::is.formula(rho_model))stop("rho_model must be a formula")
  if(!formula.tools::is.one.sided(rho_model))stop("rho_model must be one sided formula")
  temp <- as.character(rho_model)
  if(length(temp)>1)stop("rho_model must be length 1")
  rho_model_vars <- all.vars(rho_model)
  rho_vars <- all.vars(rho_model)
  if("..tag" %in% rho_vars & dt_type %in% valid_dt_type()[c(1,3)])stop("rho model cannot include ..tag variable for dt type ", dt_type)

  # check the initial values vector
  if(!is.null(all_beta.start)){
     if(!is.numeric(all_beta.start))stop("all_beta.start must be a numeric vector")
     if(!is.vector(all_beta.start)) stop("all_beta.start must be a numeric vector")
  }

  # determine the number of p1_beta parameters as determined by the p_model
  # and the number of rho parameters as determined by the rho_model
  # Since we are only modelling p1, we can use the data frame directly
  n.p1_beta <- ncol(stats::model.matrix(p_model, data=data))

  # the rho model could have a tag argument
  # we need to temporatily augment the data frame with a ..tag variable with 2 values
  # in case the rho_model includes "..tag" as a term in the model.
  temp <- data[ rep(1:nrow(data), each=2),]
  temp$..tag <- as.character(rep(1:2, nrow(data))) # must be a factor
  #browser()
  n.rho_beta <- ncol(stats::model.matrix(rho_model, data=temp))

  # determine the initial values
  if(is.null(all_beta.start))all_beta.start <- 0
  all_beta.start <- rep(all_beta.start, length.out=n.p1_beta+n.rho_beta)
  names(all_beta.start) <- paste0("b",1:length(all_beta.start))

  bbmle::parnames(LP_TL_cond_lik)=names(all_beta.start)
  fit <- bbmle::mle2(LP_TL_cond_lik,
                     start=all_beta.start,
                     vecpar=TRUE,
                     data=list(data=data, dt_type=dt_type, p_model=p_model, rho_model=rho_model, what.return="negcll", trace=trace),
                     optimizer="nlminb")
  #browser()
  summary <- data.frame(
     p_model     = paste0(as.character(p_model),  collapse=""),
     rho_model   = paste0(as.character(rho_model), collapse=""),
     name_model  = paste0("p: ", toString(p_model),
                          ";  rho: ", toString(rho_model), collapse=""),
     cond.ll    = -fit@min,
     n.parms    = n.p1_beta+n.rho_beta,
     nobs       = sum(data$freq[substr(data$cap_hist,3,4) %in% c("10","01",11)]),
     method     = "TL cond ll"
  )

  res <- list(summary=summary,
              data=data,
              p_model=p_model,
              rho_model=rho_model,
              dt_type  =dt_type,
              name_model=paste0("p: ", toString(p_model),
                                ";  rho: ", toString(rho_model), collapse=""),
              fit=fit,
              datetime=Sys.time())
  class(res) <- make_class_LP_TL(dt_type)
  res

}
