
#' Convert n1, n2, m2 to capture history data for use in estimating. Vectors are transformed to strata.
#'
#' @param n1 Number of animals marked and released. If a vector, there are the numbers in each stratum
#' @param n2 Number of animals examined for marks.
#' @param m2 Number of marked animals from n1 seen in n2.
#' @param strata Stratum labels. Added to the data frame if n1, n2, m2 are vectors.
#' @param stratum_var Name of variable to contain stratum labels.

#' @returns A data frame with the capture history, frequency.
#' \itemize{
#'   \item \strong{cap_hist} Capture history (see details below)
#'   \item \strong{freq} Number of times this capture history was observed
#' }
#' If the inputs are vectors, a stratum variable is also added.

#' @template author
#'
#' @importFrom checkmate assert_numeric assert_character

#' @examples
#'
#' # Rodli Tarn n1=109, n2=177, m2=57
#' cap_hist <- n1_n2_m2_to_cap_hist(n1=109, n2=177, m2=57)
#' cap_hist
#' fit <- Petersen::LP_fit(data=cap_hist, p_model=~..time)
#' fit$summary
#' # Now to get the estimated abundance
#' est <- Petersen::LP_est(fit, N_hat=~1)
#' est$summary
#'
#' #stratified example - Northern Pike stratified by sex
#' cap_hist <- n1_n2_m2_to_cap_hist(n1=c(4045,2777), n2=c(613,527), m2=c(89,68),
#'                 strata=c("F","M"),stratum_var="Sex")
#' cap_hist
#'

#' @export n1_n2_m2_to_cap_hist
#

n1_n2_m2_to_cap_hist <- function(n1, n2, m2,
                                 strata=paste0("Stratum_",1:length(m2)),
                                 stratum_var="Stratum"){

  checkmate::assert_numeric(n1, lower=1, finite=TRUE, any.missing=FALSE, min.len=1)
  checkmate::assert_numeric(n2, lower=1, finite=TRUE, any.missing=FALSE, min.len=1)
  checkmate::assert_numeric(m2, lower=0, finite=TRUE, any.missing=FALSE, min.len=1)
  checkmate::assert_character(strata, any.missing=FALSE, min.len=1)
  checkmate::assert_character(stratum_var, len=1, min.chars=1, any.missing=FALSE)

  # check that n1, n2, m2 all have same length
  if(sd(c(length(n1), length(n2), length(m2), length(strata)))>0){
      stop("n1,  n2,  m2 must be same length")
  }

  # create the history data frame
  df1 <- data.frame( cap_hist='10', freq=n1-m2)
  if(length(n1)>1)df1[,stratum_var] <- strata

  df2 <- data.frame( cap_hist='01', freq=n2-m2)
  if(length(n1)>1)df2[,stratum_var] <- strata

  df3 <- data.frame( cap_hist='11', freq=m2)
  if(length(n1)>1)df3[,stratum_var] <- strata

  df <- rbind(df1, df2, df3)
  df

}
