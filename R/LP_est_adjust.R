#' Estimate abundance after empirical adjustments for various factors.
#'
#' This will take a previous fit and return estimates of abundance after making various empirical adjustments
#'

#' @param N_hat Estimate of N that will be adjusted
#' @param N_hat_SE SE of the N_hat
#' @template param.conf_level
#' @param tag.retention.est Estimated tag retention probability
#' @param tag.retention.se  Estimated SE of tag retention probability
#' @param tag.reporting.est Estimated tag reporting probability
#' @param tag.reporting.se  Estimated SE of tag reporting probability
#' @param n.sim Number of simulation runs to make
#' @param trace If trace flag is set in call when estimating functions
#'
#' @return An list object with a summary data frame and a data frame with the adjustment factors.
#'
#' @details
#' The estimate and SE are converted to a beta distribution with equivalent
#' mean and SD as the estimate and se, and then simulated a large number of times
#' to get the mean and sd of all adjustments applied together.
#' Then the abundance is simulated (on the log scale), the product taken, and
#' the mean, sd, ci estimated directly.

#' @template author
#'
#' @importFrom stats rbeta rnorm sd

#' @examples
#'
#' data(data_rodli)
#' rodli.fit <- Petersen::LP_fit(data=data_rodli, p_model=~..time)
#' rodli.est <- Petersen::LP_est(rodli.fit)
#' Petersen::LP_est_adjust(rodli.est$summary$N_hat, rodli.est$summary$N_hat_SE,
#'           tag.retention.est=.90, tag.retention.se=.05)

#' @export LP_est_adjust
#'

LP_est_adjust <- function(N_hat, N_hat_SE, conf_level=0.95,
                          tag.retention.est=NULL, tag.retention.se=NULL,   # estimated tag retention and se
                          tag.reporting.est=NULL, tag.reporting.se=NULL,   # estimated tag reporting and SE
                          n.sim=10000,
                          trace=FALSE){
  # check the input estimates
  check.numeric(N_hat,    min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)
  check.numeric(N_hat_SE, min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)
  if(N_hat_SE > N_hat)warning("N_hat_SE > N_hat .. Did you accidently pass the variance of N_hat?")

  # check the confidence level
  check.conf_level(conf_level)

  # create the product of all adjustment factors
  if( is.null(tag.retention.est))adjust.tag.retention <-rep(1,n.sim)
  if(!is.null(tag.retention.est)){
    tag.retention.alpha.beta <- convert.estimate.se.to.alpha.beta (tag.retention.est, tag.retention.se)
    adjust.tag.retention <- rbeta(n.sim, shape1=tag.retention.alpha.beta[1], shape2=tag.retention.alpha.beta[2])
  }

  if( is.null(tag.reporting.est))adjust.tag.reporting <-rep(1,n.sim)
  if(!is.null(tag.reporting.est)){
    tag.reporting.alpha.beta <- convert.estimate.se.to.alpha.beta (tag.reporting.est, tag.reporting.se)
    adjust.tag.reporting <- rbeta(n.sim, shape1=tag.reporting.alpha.beta[1], shape2=tag.reporting.alpha.beta[2])
  }

  adjust <- adjust.tag.retention * adjust.tag.reporting
  adjust <- pmax(.001, adjust)  # bound away from 0 and 1

  # create the distribution of abundance on the log scale for each estimate
  estimates <- data.frame(N_hat_un    = as.vector(N_hat),
                          N_hat_un_SE = as.vector(N_hat_SE)
                          )
  estimates.adj <- plyr::adply(estimates,1, function(x){
     N_hat.dist.log <- rnorm(n.sim, mean=log(x$N_hat_un), sd=x$N_hat_un_SE/x$N_hat_un)

     # get the distribution after applying the adjustment factor
     N_hat.adjust <- exp(N_hat.dist.log + log(adjust))
     data.frame(N_hat_adj     = mean(N_hat.adjust),
                N_hat_adj_SE  = sd  (N_hat.adjust),
                N_hat_adj_LCL = quantile(N_hat.adjust, prob=(1-conf_level)/2),
                N_hat_adj_UCL = quantile(N_hat.adjust, prob=1-(1-conf_level)/2)
     )
  })

  #extract the N_hats if there are more than 1 and their respective SE
  summary <-estimates.adj

  #browser()
  adjustment <- data.frame(
                  factor=c("Tag retention"  , "Tag reporting"    ,"ALL"),
                  est   =c(ifelse(is.null(tag.retention.est),1,tag.retention.est),
                           ifelse(is.null(tag.reporting.est),1,tag.reporting.est),
                           mean(adjust)),
                  sd    =c(ifelse(is.null(tag.retention.se ),0,tag.retention.se),
                           ifelse(is.null(tag.reporting.se ),0,tag.reporting.se),
                           sd  (adjust))
  )

  res <- list(summary=summary,
              adjustment=adjustment,
              datetime=Sys.time()
              )

  res
}
