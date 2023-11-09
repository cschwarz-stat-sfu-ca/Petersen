#' Fit a Stratified-Petersen SPAS model.
#'
#' This function is a wrapper to fits a SPAS model(Schwarz, 2023; Schwarz and Taylor, 1998). Consult the SPAS package for more details.
#'
#' @template param.data
#' @param model.id Character string identifying the name of the model.
#' @template param.SPAS.autopool
#' @template param.SPAS.autopool.settings
#' @template param.SPAS.row.pool.in
#' @template param.quietly
#'
#' @importFrom SPAS SPAS.fit.model
#' @returns An list object of class *LP_SPAS_fit* with abundance estimates and other information with the following elements
#' * **summary** A data frame with the model for the capture probabilities;
#' the conditional log-likelihood; the number of parameters; the number of parameters, condition factor of the data matrix, and method used to fit the model
#' * **data** A data frame with the raw data used in the fit
#' * **fit** Results of the fit including the estimates, SE, vcov, etc.
#' * **row.pool.in**, **col.pool.in**, **autopool** Arguments used in the fit to indicate row, column, or automatic pooling used in the fit.
#' * **datetime** Date and time the fit was done
#'
#' After the fit is complete, use the *LP_SPAS_est()* function to extract the estimates, and the SPAS::SPAS.print.model() function to get a nicely
#' formatted report on the fit.
#'
#' @export LP_SPAS_fit
#' @examples
#' data(data_spas_harrison)
#'
#' fit <- Petersen::LP_SPAS_fit(data=data_spas_harrison,
#'                               model.id="Pooling rows 5/6",
#'                               row.pool.in=c(1,2,3,4,56,56),
#'                               col.pool.in=c(1,2,3,4,5,6),quietly=TRUE)
#' fit$summary
#' est <- Petersen::LP_SPAS_est(fit)
#' est$summary
#'
#' # make a nice report using the SPAS package functions
#' SPAS::SPAS.print.model(fit$fit)
#'
#' @references
#' Schwarz CJ (2023). _SPAS: Stratified-Petersen Analysis System_.
#' R package version 2023.3.31, <https://CRAN.R-project.org/package=SPAS>.
#'
#' Schwarz, C. J. and Taylor, C. G. (1998). The use of the
#' stratified-Petersen estimator in fisheries management: estimating the
#' number of pink salmon (Oncorhynchus gorbuscha) that spawn in the Fraser River.
#' Canadian Journal of Fisheries and Aquatic Sciences 55, 281-297.
#' https://doi.org/10.1139/f97-238


# this is just a wrapper to the SPAS function

LP_SPAS_fit <- function(data, model.id="Base model",
                              autopool=FALSE,
                              row.pool.in=NULL, col.pool.in=NULL,
                              min.released=100, min.inspected=50, min.recaps=50, min.rows=1, min.cols=1,
                              quietly=FALSE){


   #cat("check the data argument that it is ok\n")
   check.cap_hist_geographic.df(data)
   if(!is.character(model.id) || length(model.id) != 1)stop("invalid model id ", model.id)

   if(!is.logical(autopool) || length(autopool) !=1)stop("invalid values for autopool ", autopool)

   check.numeric(min.released,  min.value=1, max.value=Inf, req.length=1, check.whole=TRUE)
   check.numeric(min.inspected, min.value=1, max.value=Inf, req.length=1, check.whole=TRUE)
   check.numeric(min.recaps,    min.value=1, max.value=Inf, req.length=1, check.whole=TRUE)
   check.numeric(min.rows,      min.value=1, max.value=Inf, req.length=1, check.whole=TRUE)
   check.numeric(min.cols,      min.value=1, max.value=Inf, req.length=1, check.whole=TRUE)


   # make the SPAS data structure
   data.aug <- cbind(data, split_cap_hist(data$cap_hist, sep='..',  prefix='..gs', make.numeric=FALSE))

   # get the cross tabulation but the rows/columns are in the wrong order
   crosstab <- as.matrix(xtabs( freq ~ ..gs1 + ..gs2,
                      data=data.aug, exclude=NULL, na.action=na.pass))
   crosstab <- rbind(crosstab[-1,], crosstab[1,])
   crosstab <- cbind(crosstab[,-1], crosstab[,1])

   # now for the call to SPAS

   if(quietly){
     fit <- quiet.eval(SPAS::SPAS.fit.model(model.id=model.id,
                                            crosstab, autopool=autopool,
                                            row.pool.in=row.pool.in, col.pool.in=col.pool.in,
                                            row.physical.pool=FALSE,
                                            theta.pool=FALSE, CJSpool=FALSE,
                                            optMethod=c("nlminb"),
                                            optMethod.control=list(maxit = 50000),
                                            svd.cutoff=.0001, chisq.cutoff=.1,
                                            min.released  = min.released, # autopooling settings
                                            min.inspected = min.inspected,
                                            min.recaps    = min.recaps,
                                            min.rows=min.rows, min.cols=min.cols))
   }
   if(!quietly){
     fit <- SPAS::SPAS.fit.model(model.id=model.id,
                          crosstab, autopool=autopool,
                          row.pool.in=row.pool.in, col.pool.in=col.pool.in,
                           row.physical.pool=FALSE,
                          theta.pool=FALSE, CJSpool=FALSE,
                          optMethod=c("nlminb"),
                          optMethod.control=list(maxit = 50000),
                          svd.cutoff=.0001, chisq.cutoff=.1,
                          min.released  = min.released, # autopooling settings
                          min.inspected = min.inspected,
                          min.recaps    = min.recaps,
                          min.rows=min.rows, min.cols=min.cols)
   }

   summary <- data.frame(
       p_model    = paste0("Refer to row.pool/col.pool"),
       name_model = model.id,
       cond.ll    = fit$model.info$logL.cond,
       n.parms    = fit$model.info$np,
       nobs       = sum(data$freq),
       method     = "SPAS",
       cond.factor= fit$kappa.after.lp
   )

   res <- list(summary=summary,
               data=data,
               row.pool.in=row.pool.in, col.pool.in=col.pool.in, autopool=autopool,
               fit=fit,
               datetime=Sys.time())
   class(res) <- "LP_SPAS_fit"
   res
}

