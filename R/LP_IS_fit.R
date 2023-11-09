#' Fit a Lincoln-Petersen Model with incomplete stratification
#'
#' In some LP studies, stratification is only done on a random sample of unmarked
#' fish, e.g., only a sample of fish is sexed. Is is known as incomplete stratification.
#' This is a wrapper to the published code for the case of stratification by a
#' discrete covariate. At the moment, no other covariates are allowed, but see
#' the published code.
#'
#' @references Premarathna, W.A.L., Schwarz, C.J., Jones, T.S. (2018)
#' Partial stratification in two-sample capture–recapture experiments.
#' Environmetrics, 29:e2498. \doi{10.1002/env.2498}

#' @template param.data.IS
#' @template param.p_model
#' @param theta_model Model for theta (sampling fraction). Usually, this is set to be different
#' for the two sampling occasions, but you can constrain this to have equal sampling fractions at both occasions.
#' @param lambda_model Model for lambda category proportions. Usually this is set to different for the categories
#' but you can constrain this with a null matrix and the \code{logit_lamba_offset} parameter
#' @param logit_p_offset Used to fix capture probabilities at known values (seldom useful). Logit(p)=p_design %*% beta_p + logit_p_offset.
#' @param logit_theta_offset Used to fix sampling fractions at known values (seldom useful).
#' logit(theta) = theta_design %*% beta_theta + logit_theta_offset
#' @param logit_lambda_offset  Used to fix the sex ratio as a known value (e.g. .50) using
#'    logit(lambda) = lambda_design %*% beta_lambda + logit_lambda_offset.
#'    Set the design matrix to a matrix with all zeros. Notice that because the lambda proportions must sum to 1,
#'    only specify an offset matrix that is number of categories -1.
##' @param p_beta.start Initial values for call to optimization routine for the beta parameters (on the logit scale).
##' The values will be replicated to match
##' the number of initial beta parameters needed. Some care is needed here!
#' @param cat.unknown Value of character used to indicate the unknown stratum in the capture histories. Currently, this
#' is fixed to "U" regardless of what is specified.
#' @param trace If trace flag is set in call when estimating functions
#' @param control.optim Control values passed to optim() optimizer.

#' @details
#'
#' The frequency variable (\code{freq} in the \code{data} argument) is the number of animals with the corresponding capture history.
#'
#' Capture histories (\code{cap_hist} in the \code{data} argument) are character values of length 2.
#' The strata values are single character values with "U" typically representing a fish not measured
#' for the stratification variable. For example, consider the case where only a sample of unmarked fish
#' are examined for sex (M or F). Possible capture histories are:
#'
#' \itemize{
#'   \item \strong{M0}  Animals tagged and sexed as male but never seen again.
#'   \item \strong{MM}  Animals tagged and sexed as male and recaptured and tag present at event 2.
#'   \item \strong{0M}  Animals captured at event 2 that appears to be untagged and was sexed as male.
#'   \item \strong{F0}  Animals tagged and sexed as female but never seen again.
#'   \item \strong{FF}  Animals tagged and sexed as female and recaptured and tag present at event 2.
#'   \item \strong{0F}  Animals captured at event 2 that appears to be untagged and was sexed as female.
#'   \item \strong{U0}  Animals tagged and not sexed but never seen again.
#'   \item \strong{UU}  Animals tagged and not sexed and recaptured and tag present at event 2.
#'   \item \strong{0U}  Animals captured at event 2 that appears to be untagged and was not sexed.
#' }
#'
#' Capture histories such as \strong{UF} or \strong{UM} are not allowed since only UNTAGGED animals
#' are examined and sexed. Similarly, capture histories such as \strong{FM} or \strong{MF} are not allowed.
#'
#' @returns An list object of class *LP_IS_fit* with abundance estimates and other information with the following elements
#' * **summary** A data frame with the model for the capture probabilities, the sampling fractions at each capture occasion, and the category proportions;
#' the conditional log-likelihood; the number of parameters; the number of parameters, and method used to fit the model
#' * **data** A data frame with the raw data used in the fit
#' * **fit** Results of the fit including the estimates, SE, vcov, etc.
#' * **fit.call** Arguments used in the fit
#' * **datetime** Date and time the fit was done

#' @template author
#'
#' @importFrom formula.tools is.one.sided
#' @importFrom plyr is.formula
#' @importFrom stats as.formula model.matrix coef
#' @importFrom numDeriv hessian
#' @importFrom Matrix bdiag
#' @importFrom msm deltamethod


#' @examples
#'
#' data(data_wae_is_short)
#' fit <- Petersen::LP_IS_fit(data=data_wae_is_short, p_model=~..time)
#' fit$summary
#' est <- LP_IS_est(fit, N_hat=~1)
#' est$summary
#'
#' @export LP_IS_fit
#'


LP_IS_fit <- function(data, p_model,
                            theta_model=~-1+..time,
                            lambda_model=~-1+..cat,
                            logit_p_offset=0,
                            logit_theta_offset=0,
                            logit_lambda_offset=0,
                            cat.unknown="U",
                            p_beta.start=NULL, trace=FALSE,
                            control.optim=list(trace=0,maxit=1000)){
  # Fit a model for the capture probability with incomplete stratification
  # This is just a wrapper to the published code. At the moment, we do not allow
  # for other covariates in the p_model

  # check the data frame
  check.cap_hist_IS.df(data)

  # check the model for p. Must be a formula and all the variables must be in data frame
  if(!plyr::is.formula(p_model))stop("p_model must be a formula")
  if(!formula.tools::is.one.sided(p_model))stop("p_model must be one sided formula")
  temp <- as.character(p_model)
  if(length(temp)>1)stop("p_model must be length 1")
  p_model_vars <- all.vars(p_model)
  temp <- p_model_vars %in% c("..time","..cat", names(data))
  if(!all(temp))stop("p_model refers to variables not in data :",
                     paste(p_model_vars[!temp],sep=", ", collapse=""))
  # at the moment only allow ..time and ..cat in model for p
  temp <- p_model_vars %in% c("..time","..cat")
  if(!all(temp))stop("Sorry, current implementationd does not allow p_model to refers to other variables in data frame. :",
                     paste(p_model_vars[!temp],sep=", ", collapse=""))

  check.numeric(logit_p_offset,     min.value=0, max.value=0, req.length=1) # no real checking here
  check.numeric(logit_theta_offset, min.value=0, max.value=0, req.length=1)
  check.numeric(logit_lambda_offset,min.value=0, max.value=0, req.length=1)

  theta_model_vars <- all.vars(theta_model)
  temp <- theta_model_vars %in% c("..time")
  if(!all(temp))stop("theta_model refers to variables other than ..time :",
                     paste(theta_model_vars[!temp],sep=", ", collapse=""))

  lambda_model_vars <- all.vars(lambda_model)
  temp <- lambda_model_vars %in% c("..cat")
  if(!all(temp))stop("lambda_model refers to variables other than ..cat:",
                     paste(lambda_model_vars[!temp],sep=", ", collapse=""))

#  if(!is.null(p_beta.start)){
#     if(!is.numeric(p_beta.start))stop("p_beta.start must be a numeric vector")
#     if(!is.vector(p_beta.start)) stop("p_beta.start must be a numeric vector")
#  }

  if(!is.character(cat.unknown))stop("Code for unknown category must be single character of length 1")
  if(!length(cat.unknown)==1)   stop("Code for unknown category must be single character of length 1")
  if(!nchar(cat.unknown) ==1)   stop("Code for unknown category must be single character of length 1")
  if(cat.unknown != "U")        stop("Current implementation required unknown category to be 'U' - sorry")

  model.id = ""

  # give  the required design matrices for capture recapture probabilities,
  # theta(sampling(sexing) fractions ) and lambda(Category proportion)
  cat <- c(split_cap_hist(data$cap_hist, sep=""))
  strata <- sort(setdiff(unique(cat), c("0", cat.unknown))) # ignore the 0 code and U code

  dummy.data <- expand.grid( ..cat=as.factor(strata), ..time=as.factor(c(1,2)))
  captureDM <- model.matrix(p_model, dummy.data )
  #captureDM = create.DM(c(1,2,3,4))

  dummy.data <- expand.grid( ..time=as.factor(c(1,2)))
  thetaDM <- model.matrix(theta_model, dummy.data )
#  thetaDM   = create.DM(c(1,2))

  dummy.data <- expand.grid( ..cat=as.factor(strata))
  dummy.data <- dummy.data[-nrow(dummy.data),,drop=FALSE] # drop last category because must sum to 1
  lambdaDM   <- model.matrix(lambda_model, dummy.data )
  lambdaDM   <- lambdaDM[,-ncol(lambdaDM),drop=FALSE]
  #lambdaDM  = create.DM(c(1))

  #give the offset vectors(vectors of zero's should be given since no restriction)
  captureOFFSET = rep(logit_p_offset,     length.out=2*length(strata))  # for the two sampling occasions x number of categories
  thetaOFFSET   = rep(logit_theta_offset, length.out=2)  # for the two sampling occasions
  lambdaOFFSET  = rep(logit_lambda_offset,length.out=length(strata)-1)  # since category proportions must add to 1

  # set up the data structure
  Data = NULL
  Data$History   = data$cap_hist # capture histories of the raw data
  Data$counts    = data$freq     # counts of the corresponding capture histories
  Data$category  = strata        # categories in the data set

  model.id = paste0("p: ",      toString(p_model),
                         ";  theta: ",  toString(theta_model),
                         ";  lambda: ", toString(lambda_model),
                         ";  offsets ", paste0(c(captureOFFSET, thetaOFFSET, lambdaOFFSET),collapse=","),
                         collapse="")

  fit <- IS.fit.model(model.id=model.id,
                     Data     =Data,
                     captureDM=captureDM,
                     thetaDM  =thetaDM,
                     lambdaDM =lambdaDM,
                     captureOFFSET=logit_p_offset,
                     thetaOFFSET  =logit_theta_offset,
                     lambdaOFFSET =logit_theta_offset,
                     control.optim=control.optim)

  #browser()
  summary <- data.frame(
     p_model     = toString(p_model),
     theta_model = toString(theta_model),
     lambda_model= toString(lambda_model),
     name_model = model.id,
     cond.ll    = -fit$NLL,
     n.parms    = fit$np,
     nobs       = sum(data$freq),
     method     = "full ll"
  )

  res <- list(summary=summary,
              data=data,
              p_model=p_model,
              name_model=model.id,
              fit=fit,
              fit.call=list(
                          Data=Data,
                          strata=strata,
                          cat.unknown=cat.unknown,
                          p_model     =p_model,
                           logit_p_offset=logit_p_offset,
                          logit_theta_offset=logit_theta_offset,
                          logit_lambda_offset=logit_lambda_offset),
              datetime=Sys.time())
  class(res) <- "LP_IS_fit"
  res

}
