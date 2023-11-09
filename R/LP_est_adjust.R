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
#' @param n1.adjust.est     Adjustment to "n1". This should typically be a ratio of new n1 to old n1
#' @param n1.adjust.se      Adjustment to "n1" uncertainty
#' @param n2.adjust.est     Adjustment to "n2"  This should typically be a ratio of new n2 to old n2
#' @param n2.adjust.se      Adjustment to "n2" uncertainty
#' @param m2.adjust.est     Adjustment to "m2"  This should typically be a ratio of new m2 to old m2
#' @param m2.adjust.se      Adjustment to "m2" uncertainty
#'
#' @param n.sim Number of simulation runs to make
#' @param trace If trace flag is set in call when estimating functions
#'
#' @returns An list object with a summary data frame and a data frame with the adjustment factors with the following objects
#' **summary** A data frame with the adjusted abundance estimates, SE, and CI
#' **adjustment** a data frame showing the adjustment factors applied for tag retention, tag reporting, n1 n2 or m2.
#' **datetime** Date and time the adjustment was done

#'
#' @details
#' The estimate and SE are converted to a beta distribution for adjustment factors between 0 and 1 with equivalent
#' mean and SD as the estimate and se. The estimate and se are used in normal distribution for adjustment factors for n1, n2, and m2.
#' These adjustment factors are then simulated a large number of times and then multiplied together
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
#' res <- Petersen::LP_est_adjust(rodli.est$summary$N_hat, rodli.est$summary$N_hat_SE,
#'           tag.retention.est=.90, tag.retention.se=.05)
#' res$summary

#' @export LP_est_adjust
#'

LP_est_adjust <- function(N_hat, N_hat_SE, conf_level=0.95,
                          tag.retention.est=1, tag.retention.se=0,   # estimated tag retention and se
                          tag.reporting.est=1, tag.reporting.se=0,   # estimated tag reporting and SE
                          n1.adjust.est    =1, n1.adjust.se    =0,
                          n2.adjust.est    =1, n2.adjust.se    =0,
                          m2.adjust.est    =1, m2.adjust.se    =0,
                          n.sim=10000,
                          trace=FALSE){
  # check the input estimates
  check.numeric(N_hat,    min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)
  check.numeric(N_hat_SE, min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)
  if(N_hat_SE > N_hat)warning("N_hat_SE > N_hat .. Did you accidently pass the variance of N_hat?")

  check.numeric(tag.retention.est, min.value=0, max.value=1, req.length=1, check.whole=FALSE)
  check.numeric(tag.retention.se , min.value=0, max.value=1, req.length=1, check.whole=FALSE)

  check.numeric(tag.reporting.est, min.value=0, max.value=1, req.length=1, check.whole=FALSE)
  check.numeric(tag.reporting.se , min.value=0, max.value=1, req.length=1, check.whole=FALSE)

  check.numeric(n1.adjust.est, min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)
  check.numeric(n1.adjust.se , min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)

  check.numeric(n2.adjust.est, min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)
  check.numeric(n2.adjust.se , min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)

  check.numeric(m2.adjust.est, min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)
  check.numeric(m2.adjust.se , min.value=0, max.value=Inf, req.length=1, check.whole=FALSE)

  # check the confidence level
  check.conf_level(conf_level)


  # create the product of all adjustment factors
  if( tag.retention.se==0)adjust.tag.retention <-rep(tag.retention.est,n.sim)
  if(!tag.retention.se==0){
    tag.retention.alpha.beta <- convert.estimate.se.to.alpha.beta (tag.retention.est, tag.retention.se)
    adjust.tag.retention <- rbeta(n.sim, shape1=tag.retention.alpha.beta[1], shape2=tag.retention.alpha.beta[2])
  }

  if( tag.reporting.se==0)adjust.tag.reporting <-rep(tag.reporting.est,n.sim)
  if(!tag.reporting.se==0){
    tag.reporting.alpha.beta <- convert.estimate.se.to.alpha.beta (tag.reporting.est, tag.reporting.se)
    adjust.tag.reporting <- rbeta(n.sim, shape1=tag.reporting.alpha.beta[1], shape2=tag.reporting.alpha.beta[2])
  }

  if( n1.adjust.se==0)adjust.n1 <-rep(n1.adjust.est,n.sim)
  if(!n1.adjust.se==0){
    adjust.n1 <- rnorm(n.sim, mean=n1.adjust.est, sd=n1.adjust.se)
  }

  if( n2.adjust.se==0)adjust.n2 <-rep(n2.adjust.est,n.sim)
  if(!n2.adjust.se==0){
    adjust.n2 <- rnorm(n.sim, mean=n2.adjust.est, sd=n2.adjust.se)
  }

  if( m2.adjust.se==0)adjust.m2 <-rep(m2.adjust.est,n.sim)
  if(!m2.adjust.se==0){
    adjust.m2 <- rnorm(n.sim, mean=m2.adjust.est, sd=m2.adjust.se)
  }

  adjust <- adjust.tag.retention *
            adjust.tag.reporting *
            adjust.n1            *
            adjust.n2            /
            adjust.m2
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
                  factor=c("Tag retention"  ,
                           "Tag reporting"  ,
                           "n1",
                           "n2",
                           "m2",
                           "ALL"
                           ),
                  est   =c(mean(adjust.tag.retention),
                           mean(adjust.tag.reporting),
                           mean(adjust.n1),
                           mean(adjust.n2),
                           mean(adjust.m2),
                           mean(adjust)
                           ),
                  sd    = c(sd(adjust.tag.retention),
                           sd(adjust.tag.reporting),
                           sd(adjust.n1),
                           sd(adjust.n2),
                           sd(adjust.m2),
                           sd  (adjust)
                           )
  )

  res <- list(summary=summary,
              adjustment=adjustment,
              datetime=Sys.time()
              )

  res
}
