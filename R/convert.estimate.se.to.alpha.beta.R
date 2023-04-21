#' Convert an estimate and SE to (alpha,beta) for a beta distribution.
#'
#' This computes the (alpha, beta) so that a beta distribution will
#' have the mean equal to the estimate and a SD equal to the SE.
#'
#' @param estimate A probability between 0 and 1.
#' @param SE Its standard rror

#' @return Vector of (alpha, beta)
#' @examples
#'
#' # Suppose that estimate was .90 with se of .05
#' ab <- convert.estimate.se.to.alpha.beta(.9, .05)
#' temp <- rbeta(1000, shape1=ab[1], shape2=ab[2])
#' mean(temp)
#' sd  (temp)
#'
#' @noRd

convert.estimate.se.to.alpha.beta <- function(estimate, SE){

  check.numeric(estimate, min.value=0, max.value=1,  check.whole=FALSE)
  check.numeric(SE      , min.value=0, max.value=.3, check.whole=FALSE)

  # use method of moments
  #   estimate = alpha/(alpha+beta)
  #   SE^2     = alpha*beta/(alpha+beta)^2/(alpha+beta+1)

  # see moment estimator at https://en.wikipedia.org/wiki/Beta_distribution
  alpha = (estimate *(1-estimate)/SE^2 -1) * estimate
  beta  = (estimate *(1-estimate)/SE^2 -1) * (1-estimate)
  c('alpha'=alpha, "beta"=beta)
}

