#' Compute summary statistics from the capture histories
#'
#' This function takes the capture histories and computes $n_1$, $n_2$ and $m_2$.
#'

#' @template param.data
#'
#' @return Summary statistics of *n1* (number observed at first occasion), *n2* (number observed at second occasion),
#' and *m2* number recaptured in second occasion.
#'

#' @examples
#'
#' data(data_rodli)
#' LP_summary_stats(data_rodli)
#'
#' @export LP_summary_stats
#'

LP_summary_stats <- function(data){
  # Compute summary statistics from capture-recapture data

  # check the data frame
  check.cap_hist.df(data)

  data <- cbind(data, split_cap_hist(data$cap_hist))
  res <- data.frame(
        n1 = sum(data$freq*(data$t1=="1")),
        n2 = sum(data$freq*(data$t2=='1')),
        m2 = sum(data$freq*(data$t1=='1' & data$t2=="1"))
        )
  # compute the recapture proportion and the marked fraction
  res$p.recap <- res$m2 / res$n1
  res$mf      <- res$m2 / res$n2

  res
}

#LP_summary_stats(rodli)
