#' Fit a Lincoln-Petersen Model using conditional likelihood
#'
#' This will take a data frame of capture histories, frequencies, and
#' additional covariates (e.g., strata and/or continuous covariates) and the model
#' for the capture probabilities and will use conditional likelihood (Huggins, 1989)
#' to fit the model. The population abundance is estimated using
#' a Horvitz-Thompson type estimator and the user can request abundance
#' estimates for sub-sets of the population.
#'

#' @template param.data
#' @template param.p_model
#' @param p_beta.start Initial values for call to optimization routine for the beta parameters (on the logit scale).
#' The values will be replicated to match
#' the number of initial beta parameters needed. Some care is needed here!
#' @param trace If trace flag is set in call when estimating functions

#' @details
#'
#' The frequency variable (\code{freq} in the \code{data} argument) is the number of animals with the corresponding capture history.
#'
#' Capture histories (\code{cap_hist} in the \code{data} argument) are character values of length 2.
#' \itemize{
#'   \item \strong{10}  Animals tagged but never seen again.
#'   \item \strong{11}  Animals tagged and recaptured and tag present at event 2.
#'   \item \strong{01}  Animals captured at event 2 that appear to be untagged.
#' }
#'
#'
#' @returns An list object of class *LP_fit* with abundance estimates and other information with the following elements
#' * **summary** A data frame with the model for the capture probabilities;
#' the conditional log-likelihood; the number of parameters; and method used to fit the model
#' * **data** A data frame with the raw data used in the fit
#' * **fit** Results of the fit from the optimizer
#' * **datetime** Date and time the fit was done
#'
#' After the fit is done, use the *LP_est()* function to get estimates of abundance.

#' @template author
#'
#' @importFrom formula.tools is.one.sided
#' @importFrom plyr is.formula
#' @importFrom stats as.formula model.matrix coef
#' @importFrom bbmle mle2

#' @example man-roxygen/ex.LP.R
#'
#' @export LP_fit
#'
#' @references
#' Huggins, R. M. 1989. On the Statistical Analysis of Capture Experiments.
#' Biometrika 76: 133--40.


#'

LP_fit <- function(data, p_model=~..time, p_beta.start=NULL, trace=FALSE){
  # Fit a model for the capture probability using conditional likelihood

  # check the data frame
  check.cap_hist.df(data)

  # check the model for p. Must be a formula and all the variables must be in data frame
  if(!plyr::is.formula(p_model))stop("p_model must be a formula")
  if(!formula.tools::is.one.sided(p_model))stop("p_model must be one sided formula")
  temp <- as.character(p_model)
  if(length(temp)>1)stop("p_model must be length 1")
  p_model_vars <- all.vars(p_model)
  temp <- p_model_vars %in% c("..time",names(data))
  if(!all(temp))stop("p_model refers to variables not in data :",
                     paste(p_model_vars[!temp],sep=", ", collapse=""))

  if(!is.null(p_beta.start)){
     if(!is.numeric(p_beta.start))stop("p_beta.start must be a numeric vector")
     if(!is.vector(p_beta.start)) stop("p_beta.start must be a numeric vector")
  }

  # determine the number of p_beta parameters as determined by the p_model
  # we need to temporatily augment the data frame with a time variable with 2 values
  # in case the p_model includes "..time" as a term in the model.
  temp <- data[ rep(1:nrow(data), each=2),]
  temp$..time <- as.character(rep(1:2, nrow(data))) # must be a factor
  #browser()
  n.beta <- ncol(stats::model.matrix(p_model, data=temp))

  # fit the model for the capture-probabilities
  # determine the initial values

  if(is.null(p_beta.start))p_beta.start <- 0
  p_beta.start <- rep(p_beta.start, length.out=n.beta)
  names(p_beta.start) <- paste0("b",1:length(p_beta.start))

  bbmle::parnames(LP_cond_lik)=names(p_beta.start)
  fit <- bbmle::mle2(LP_cond_lik,
                     start=p_beta.start,
                     vecpar=TRUE,
                     data=list(data=data,
                               p_model=p_model,
                               what.return="negcll",
                               trace=trace),
                     optimizer="nlminb")

  #browser()
  summary <- data.frame(
     p_model    = paste0(as.character(p_model), collapse=""),
     name_model = paste0("p: ", toString(p_model), collapse=""),
     cond.ll    = -fit@min,
     n.parms    = n.beta,
     nobs       = sum(data$freq),
     method     = "cond ll"
  )

  res <- list(summary=summary,
              data=data,
              p_model=p_model,
              name_model=paste("p: ", toString(p_model), collapse=""),
              fit=fit,
              datetime=Sys.time())
  class(res) <- "LP_fit"
  res

}
