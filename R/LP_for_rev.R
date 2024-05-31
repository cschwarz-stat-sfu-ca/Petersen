#' Fit a combined FORWARD and REVERSE simple Lincoln-Petersen Model using pseudo-likelihood
#'
#' EXPERIMENTAL. This will take a data frame of capture histories, frequencies, and
#' additional covariates (e.g., strata and/or continuous covariates) for a simple forward Petersen
#' estimate plus estimates of escapement
#' and associated stock proportions with SE for backwards estimation.
#' DO NOT USE YET.
#'

#' @template param.data
#' @param E Escapement at one or more terminal areas. E, E.SE, G, G.SE must all have the same length
#' @param E.SE SE of the estimates of escapement
#' @param G Estimated proportion of the stock at the first capture location estimated using GSI and other methods
#' @param G.SE SE of the estimated stock proportion.
#' @param min.G Miniumum acceptable stock proportion during the bootstrap estimation of uncertainty
#' @param n.boot Number of bootstrap samples used to estimate the uncertainty
#' @param trace Should intermediate tracing be enabled (e.g. browser() stops)


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
#' A pseudo-likelihood is constructed consisting of the usual likelihood for a forward capture recapture and
#' marginal likelihoods for each of the escapement (E) and stock proportions (G) point estimates. I have not
#' integrated over the uncertainty in G and E.
#'
#' @returns An list object of class *LP_for_rev_est* with abundance estimates and measures of uncertainty
#' * **summary** A data frame with the pseudo-likelihood value, abundance estimates and SE
#' * **data** A data frame with the raw data used in the fit
#' * **fit** Results of the fit from the optimizer
#' * **datetime** Date and time the fit was done
#'
#' Unlike other routine, it is not necessary to use a *XX_est()* function to get estimates of abundance.

#' @template author
#'
#' @importFrom formula.tools is.one.sided
#' @importFrom plyr is.formula
#' @importFrom stats as.formula model.matrix coef
#' @importFrom bbmle mle2

#' @example man-roxygen/ex.for_rev_LP.R
#'
#' @export LP_for_rev_fit
#'

#'

LP_for_rev_fit <- function(data,
                           E,
                           E.SE,
                           G,
                           G.SE,
                           min.G=.01,
                           n.boot=100,
                           trace=FALSE){
  # Fit a model for the capture probability using conditional likelihood

  # check the data frame
  check.cap_hist.df(data)

  # check the E, E.SE, G, and G.SE parameters
  # Must all be the same length
  req.length=max(1, length(E), length(E.SE), length(G), length(G.SE))
  check.numeric(E,    min.value=1,     max.value=Inf,     req.length=req.length, check.whole=TRUE)
  check.numeric(E.SE, min.value=1,     max.value=Inf,     req.length=req.length, check.whole=FALSE)
  check.numeric(G,    min.value=min.G, max.value=1-min.G, req.length=req.length, check.whole=FALSE)
  check.numeric(G.SE, min.value=0    , max.value=Inf,     req.length=req.length, check.whole=FALSE)

  # define the negative of the log-likelihood functions
  # forward capture recaptured
  f.nloglike <- function(theta, f.data){
    # negative of forward log-likelihood given f.data$cap_hist and f.data$freq
    # theta = log(N), p1, p2
    N  <- exp(theta[1])
    p1 <- theta[2]
    p2 <- theta[3]

    f.data$p.hist <- 1
    f.data$p.hist <- f.data$p.hist * p1^(substr(f.data$cap_hist,1,1)=='1') * (1-p1)^(substr(f.data$cap_hist,1,1)=='0')
    f.data$p.hist <- f.data$p.hist * p2^(substr(f.data$cap_hist,2,2)=='1') * (1-p2)^(substr(f.data$cap_hist,2,2)=='0')
    p00 <- (1-p1)*(1-p2)

    n <- sum(f.data$freq)
    n00 <- N-n
    #browser()
    logL <- lgamma(N+1) - lgamma(N+1-n) +
      sum(f.data$freq * log(f.data$p.hist))+
      n00         * log(p00)
    -logL
  }

  # reverse capture-recapture negative log-likelihood is a simple binomial for each E, G pair.
  r.nloglike <- function(theta, r.data){
    # negative of reverse log-likelihood given r.data=E,  E.RSE G (escapement and GSI proportion for E)
    # Multiple values of E and G are allowed
    # theta <- log(N)
    N <- exp(theta[1])

    E <- r.data$E
    G <- r.data$G

    logL <- lgamma(N+1)- lgamma(N+1-E)+
      sum(E*log(G))+
      sum((N-E)*log(1-G))
    -logL
  }

  # both negative likelihoods combined
  b.nloglike <- function(theta, b.data){
    # both direction pseudo likelihood
    # theta = log(N) p1 p2
    # bdata has f.data and r.data terms

    f.data <- b.data$data
    r.data <- list(E=b.data$E, G=b.data$G)
    #browser()
    b.nLL <- f.nloglike(theta,    f.data) +
      r.nloglike(theta[1], r.data)
  }

  # parameters to the pseudo-likelihood are log(N), p1, p2
  n.beta <- 3

  # fit the model for the capture-probabilities using a simple LP to determine the initial values
  p_model=~..time # this should likely be passed as a parameter
  f.res <- Petersen::LP_fit(data, p_model=p_model)
  f.est <- Petersen::LP_est(f.res)


  beta.start <- c(log(f.est$summary$N_hat), f.est$detail$data.expand$p[1:2])

  fit <- optim(beta.start, b.nloglike, b.data=list(data=data, E=E, G=G),
               lower=.01, upper=c(Inf,.99, .99), method="L-BFGS-B")
  fit$par
  exp(fit$par[1])

  #browser()
  summary <- data.frame(
     p_model    = paste0(as.character(p_model), collapse=""),
     name_model = paste0("p: ", toString(p_model), " ForRev ", collapse=""),
     pseudo.ll  = -fit$value,
     n.parms    = n.beta,
     nobs       = sum(data$freq),
     method     = "pseudo ll",
     N_hat      = exp(fit$par[1])
  )

  res <- list(summary=summary,
              data=list(data=data, E=E, E.SE=E.SE, G=G, G.SE=G.SE),
              name_model="forward_reverse Petersen",
              fit=fit,
              datetime=Sys.time())
  class(res) <- "LP_for_rev_est"
  res

}

