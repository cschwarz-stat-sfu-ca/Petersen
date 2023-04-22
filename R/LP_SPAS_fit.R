#' Fit a Stratified-Petersen SPAS model.
#'
#' This function is a wrapper to fits a SPAS model. Consult the SPAS package for more details.
#'
#' @template param.data
#' @param model.id Character string identifying the name of the model.
#' @template param.SPAS.autopool
#' @template param.SPAS.autopool.settings
#' @template param.SPAS.row.pool.in
#'
#' @importFrom SPAS SPAS.fit.model
#' @return A list with many entries including the fit from SPAS and a summary line
#' in standard format. Refer to the vignettes for more details.
#' @export LP_SPAS_fit
#' @examples
#' data(data_spas_harrison)
#'
#' mod1 <- Petersen::LP_SPAS_fit(data=data_spas_harrison,
#'                       model.id="Pooling rows 5/6",
#'                       row.pool.in=c(1,2,3,4,56,56),
#'                       col.pool.in=c(1,2,3,4,5,6))
#' SPAS::SPAS.print.model(mod1$fit)
#' mod1$summary

# this is just a wrappr to the SPAS function

LP_SPAS_fit <- function(data, model.id="Base model",
                              autopool=FALSE,
                              row.pool.in=NULL, col.pool.in=NULL,
                              min.released=100, min.inspected=50, min.recaps=50, min.rows=1, min.cols=1){


   cat("check the data argument that it is ok")

   # make the SPAS data structure
   data.aug <- cbind(data, split_cap_hist(data$cap_hist, sep='..',  prefix='..gs', make.numeric=FALSE))

   # get the cross tabulation but the rows/columns are in the wrong order
   crosstab <- as.matrix(xtabs( freq ~ ..gs1 + ..gs2,
                      data=data.aug, exclude=NULL, na.action=na.pass))
   crosstab <- rbind(crosstab[-1,], crosstab[1,])
   crosstab <- cbind(crosstab[,-1], crosstab[,1])

   # now for the call to SPAS
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

