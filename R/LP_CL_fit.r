
#' Fit the Chen-Lloyd model to estimate abundance using a non-parametric smoother for a covariates
#'
#' This will take a data frame of capture histories, frequencies, and a covariates
#' and will do a non-parametric smoother for the detection probabilities as a function
#' of the covariates and use this to estimate the population size.
#'

#' @template param.data
#' @param covariate Name of continuous covariate that influences capture probabilities at each event
#' @param centers   Centers of bins to group the covariates. We suggest no more than 30 bins in total
#' with fewer bins with smaller sample sizes. Of course with smaller sample sizes, a simple stratified
#' estimator may be easier to use.
#' @param h1,h2 Standard deviation of normal kernel for first sampling event. This should be between 1/2 and the 1.5x the
#' bin width. Larger values imply more smoothing. Smaller values imply less smoothing.
#' @template param.conf_level

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
#' @returns An list object of class *LP_CL_fit* with abundance estimates and other information with the following elements
#' * **summary** A data frame  with the estimates of abundance, SE, and CI
#' * **fit** Details on the Chen and Lloyd fit including the smoothed estimates of catchability, estimates abundance by category classes,
#' estimates of total abundance, plots of the estimated abundance curve and catchability curves, etc.
#' * **datetime** Date and time the fit was done
#'
#' @import ggplot2
#' @importFrom graphics hist
#' @importFrom stats quantile dnorm
#' @importFrom rlang .data

#' @examples
#'
#' library(Petersen)
#' data(data_NorthernPike)
#' res <- LP_CL_fit(data_NorthernPike, "length")
#' res$summary

#' @export LP_CL_fit
#' @references
#' SX Chen, CJ Lloyd (2000).
#' A nonparametric approach to the analysis of two-stage mark-recapture experiments.
#' Biometrika, 87, 633–649. \doi{10.1093/biomet/87.3.633}.



# These functions were supplied by Lloyd, C. and Chen, S. in regards to their
# paper in Biometrika. Received on 2005-12-03
#
# Documentation and integration into CL_fit added by C.Schwarz, 2023-02-22

LP_CL_fit <- function(data, covariate,
                   centers=hist(data[,covariate,drop=TRUE], breaks="Sturges", plot=FALSE)$mids,
                   h1=(centers[2]-centers[1])*.75,
                   h2=(centers[2]-centers[1])*.75,
                   conf_level=0.95){

  # check the data frame
  check.cap_hist.df(data)

  # check that covariate is in the data frame
  if(!is.character(covariate))   stop("Covariate name must be a character value")
  if(!length(covariate)==1)      stop("Only a single covariate is allowed")
  if(!covariate %in% names(data))stop("Covariate not in data frame: ", covariate)

  check.conf_level(conf_level)

# define functions
ChenLloydestimate <- function(data, h1, h2, w = FALSE, covariate="Covariate")
{
#       Estimate the effect of heterogeneity upon the petersen estimator
#       Data consists of a list "data" with elements
#           data$MM - number of marks placed out in the covariate class with respective center in data$centers
#           data$mm - number of recaptures from this covariate class with respective center in data$centers
#           data$uu - number of UNmarked   from this covariate class
#           data$center - center of this covariate class
#       An easy way to create this data from a list of capture histories and their covariates is
#       given by the "bin"  function
#
#       I'm not sure what the "w" parameter measures, but it appears to be some sort of smoothing using
#       a weighted average.
#
#       The smoothing parameters are h1 and h2
#           h1 = std dev of the smoothing kernel for the inital probability of capture rates
#           h2 = std dev of the smoothing kernel for the second probability of captur
  #browser()
  S <- function(x, h, xx = x, wgt = 0. * x + 1.){
#       Create the kernel density weights using a Normal distribution with std dev=h
	      k <- length(x)
	      M <- matrix(rep(xx, k), ncol = k) - matrix(rep(x, length(xx)), ncol = k, byrow = T)
	      M <- matrix(dnorm(M, sd = h), ncol = k) %*% diag(wgt)
	      sweep(M, 1., apply(M, 1., sum), "/")
  }

	MM <- as.vector(data$M)
	mm <- as.vector(data$m)
	uu <- as.vector(data$u)
	x <- data$center

#       Get the smoothing matrix based on a normal kernel with bandwidth = std dev = h1
	S1 <- S(x, h1)
#       Smooth the marks put out and the marks recaptured to estimate the recapture rates (smoothed)
	Mtl <- as.vector(S1 %*% MM)
	mtl <- as.vector(S1 %*% mm)
	S2 <- S(x, h2, wgt = w/mtl + (1. - w))
	ptl <- mtl/Mtl     # smoothed recapture rates
  data$ptl <- ptl

#       Estimate number in each population at second time period based on smoothed recapture rate
	N1 <- uu/ptl + MM
	A <- as.vector(S2 %*% (mtl * N1))
	B <- as.vector(S2 %*% mtl)
	N2 <- A/B
	dM <- diag(MM)
	SU <- as.vector(S2 %*% (N1 - MM))
	#       Variance matrix for N1
	Vm <- S1 %*% diag(MM * ptl * (1. - ptl)) %*% t(S1)
	VLA <- diag(SU * mtl * (Mtl - mtl)) + dM %*% Vm %*% dM
	VLA1 <- diag(1./A) %*% VLA %*% diag(1./A)
	VLB1 <- diag(1./B) %*% Vm %*% diag(1./B)
	CLAB1 <- diag(1./A) %*% dM %*% Vm %*% diag(1./B)
	VLN1 <- VLA1 + VLB1 - CLAB1 - t(CLAB1)
	VLA2 <- diag(1./A) %*% S2 %*% VLA %*% t(S2) %*% diag(1./A)
	VLB2 <- diag(1./B) %*% S2 %*% Vm %*% t(S2) %*% diag(1./B)
	CLAB2 <- diag(1./A) %*% S2 %*% dM %*% Vm %*% t(S2) %*% diag(
		1./B)
	VLN2 <- VLA2 + VLB2 - CLAB2 - t(CLAB2)
#
#       Now to assemble the output from the function to be returned to the user

	value <- NULL
	value$p <- ptl
	value$N1 <- N1
	value$N2 <- N2
	value$sd1j <- sqrt(diag(VLN1)) * N1
	value$sd2j <- sqrt(diag(VLN2)) * N2
	one <- matrix(rep(1., length(MM)), ncol = 1.)
	s1 <- sqrt(t(one) %*% diag(N1) %*% VLN1 %*% diag(N1) %*% one)
	s2 <- sqrt(t(one) %*% diag(N2) %*% VLN2 %*% diag(N2) %*% one)
#
#       Compute the simple petersen estimator assuming homogeneity of capture
#
	Nhom <- (sum(uu) * sum(MM))/sum(mm) + sum(MM)
        sdhom <-  sqrt(Nhom) * sqrt( sum(MM-mm)*sum(uu))/sum(mm)

#       The three estimates that are returned are the
#           simple N1 which is empirical area under curve; and
#           smoothed N2 which is area under smoothed curve; and
#           simple pooled Petersen (assuming homogeneity of capture,
	value$total.est <- c(sum(N1), sum(N2), Nhom )
	value$total.sd  <- c(s1,          s2,  sdhom)
	value$total.explain<- c("Area under empirical curve",
	                        "Area under smoothed curve",
	                        "Petersen")

	#  graphical ouput
	plotdata <- data.frame(
	                center=data$center,
	                M     =MM,
	                m     =mm,
	                N1    =value$N1,
	                N2    =value$N2,
	                sd1j  =value$sd1j,
	                sd2j  =value$sd2j,
	                ptl   =ptl        #smoothed recapture probabilities
	)

	plot1 <- ggplot(data=plotdata, aes(x=.data$center, y=.data$m/.data$M))+
	  ggtitle("Estimated recapture probability by classes",
	          subtitle="Smoothed probability of recapture is shown")+
	  geom_point()+
	  geom_line(aes(y=.data$ptl), linetype="dashed")+
	  xlab(covariate)+
	  ylab("Recapture probability")

	#browser()
	plot2<-ggplot(data=plotdata, aes(x=.data$center, y=.data$N1))+
	  geom_point()+
	  geom_line(aes(y=.data$sd1j),linetype="dashed")+
	  geom_line(aes(y=.data$sd2j),linetype="dotted")+
	  geom_line(aes(y=.data$N2))+
	  annotate("text", label=paste("Pop (SE): ",round(value$total.est[2])," ( ", round(value$total.sd[2])," )"),
	           x=Inf, y=Inf, hjust=1.5, vjust=1.5)+
	  xlab(covariate)+ylab("Estimated population size")+
  	ggtitle("Estimated frequency distribution - estimated population size ",
                subtitle="SE for smoothed and raw estimates are shown in dashed lines at bottom")
  value$plot1 <- plot1
  value$plot2 <- plot2
	value
}

"bin" <- function(x, pts){
	# Function bins a vector of numbers into bins with
	# center in the vector pts. It returns a vector of
	# frequencies of the same dimension as pts.
        # Modified by C.schwarz, 2005-12-02 to use -1000*min and 1000*max as the upper and lower bound of first and last class
	s.p <- sort(pts)
	brks <- (c(-10000.*min(s.p), s.p) + c(s.p, 10000.*max(s.p)))/2.
	hist(x, breaks = brks, plot = FALSE)$counts
}

"plugin.F" <-function(x, B = min(length(x), 49.), print = FALSE){
	#       Date: 5.5.98
	#       BANDWIDTH SELECTOR FOR CDF ESTIMATION
	#
	#       METHOD
	#       The sixth derivative is estimated assuming normality. From
	#       this, optimal bandwidth is used to estimate the fourth derivative.
	#       This is used to optimally estimate the second derivative and
	#       then the usual formula is used, see Wand and Jones, p72.
	#       If the selected bandwidth is negative then it is replaced by
	#       the direct normal plugin rule.
	#
	#       The derivatives are estimated using binned data to reduce the
	#       computation from order n^2 to order B^2, see Wand and Jones, p188.
	#
	#       VALUE:
	#          vector with 2 components. The first is the direct normal
	#          plugin, the second the 2stage normal plugin.
	#
	n <- length(x)
	grid <- function(x, n){
		min(x, na.rm = TRUE) + (-1.:(n + 1.))/n * (max(x, na.rm =TRUE) - min(x, na.rm = TRUE))
	}
	xx <- grid(x, B - 3.)
	X <- matrix(rep(xx, rep(B, B)), ncol = B)
	X <- X - t(X)
	counts <- hist(x, breaks = grid(x, B - 2.), plot = FALSE)$counts
	k <- matrix(rep(counts, rep(B, B)), ncol = B)
	k <- k * t(k)
	qs <- quantile(x, na.rm = TRUE)
	s <- (qs[4.] - qs[2.])/1.3400000000000001
	names(s) <- NULL
	h4 <- (1.2406999999999999 * s)/n^(1./7.)
	XX <- X/h4
	psi4 <- sum(k * (XX^4. - 6. * XX^2. + 3.) * dnorm(XX))/(n *	n * h4^5.)
	h2 <- 0.79790000000000005/(psi4 * n)^(1./5.)
	XX <- X/h2
	psi2 <- sum(k * (XX^2. - 1.) * dnorm(XX))/(n * n * h2^3.)
	if(print) {
		print(paste("Psi4= ", signif(psi4, 3.), "Psi2=", signif(psi2, 3.)))
		print(paste("Bandwidths used: h4=", signif(h4, 3.),	"h2=", signif(h2, 3.)))
	}
	h2stage <- (-0.56000000000000005/psi2/n)^(1./3.)
	h0stage <- (1.583 * s)/n^(1./3.)
	if(h2stage <= 0.) {
		h2stage <- h0stage
	}
	c(h0stage, h2stage)
}

#  extract the covariate vectors for fish captured at event 1, event 2, and both events
#  don't forget to expand by the frequencey
   covar.data <- NULL
   covar.data$center <- centers
   covar.data$all <- rep(data[,covariate, drop=TRUE], times=data$freq)
   select <- substr(data$cap_hist,1,1) == '1'
   covar.data$x1  <- rep(data[select,covariate,drop=TRUE], times=data$freq[select])
   select <- substr(data$cap_hist,2,2) == '1'
   covar.data$x2  <- rep(data[select,covariate,drop=TRUE], times=data$freq[select])
   select <- data$cap_hist == '11'
   covar.data$x12 <- rep(data[select,covariate,drop=TRUE], times=data$freq[select])

#  remove all  missing values as they cannot be used.

   covar.data$x1  <- covar.data$x1  [ !is.na(covar.data$x1)  ]     # select only those fish at least 18 inches in
   covar.data$x2  <- covar.data$x2  [ !is.na(covar.data$x2)  ]
   covar.data$x11 <- covar.data$x12 [ !is.na(covar.data$x12) ]

#  create the number of marked, unmarked, and recaptured into the bins centered at "centers"
   #browser()
   covar.data$M      <- bin(covar.data$x1 , covar.data$center)   # fish tagged
   covar.data$m      <- bin(covar.data$x12, covar.data$center)   # number of marks recaptured
   covar.data$u      <- bin(covar.data$x2 , covar.data$center) - covar.data$m # unmarked at time=2

#  estimate the approximate optimal band widths
   optimal.binwidth      <- NULL
   optimal.binwidth$t1   <- plugin.F( covar.data$x1)
   optimal.binwidth$t2   <- plugin.F( covar.data$x2)
   optimal.binwidth$t1t1 <- plugin.F( covar.data$x12)

#  This shows for this problem that optimal bandwidth is about 1 inch which happens to match our increments

#  Now to estimate the population size and get plots etc.
   res <- NULL
   res$fit <- ChenLloydestimate(covar.data,
                            h1= h1, #1, #optimal.binwidth$t1[1],
                            h2= h2, #1,  #optimal.binwidth$t2[2],
                            covariate=covariate)
   #browser()
   res$fit$optimal.binwidth <- optimal.binwidth
   res$fit$covar.data       <- covar.data

   # now to create the summary output
   summary <- data.frame(
                N_hat_f     = "~1",
                N_hat_rn    = "(Intercept)",
                N_hat       = res$fit$total.est[2],   # smoothed area under the curve
                N_hat_SE    = res$fit$total.sd [2],   # smoothed area under the curve
                N_hat_conf_level=conf_level,
                N_hat_conf_method="logN")
   summary$ N_hat_LCL   = exp(log(summary$N_hat)-qnorm(1-(1-conf_level)/2)*summary$N_hat_SE/summary$N_hat)
   summary$ N_hat_UCL   = exp(log(summary$N_hat)+qnorm(1-(1-conf_level)/2)*summary$N_hat_SE/summary$N_hat)
   summary$p_model = NA
   summary$cond.ll <- NA
   summary$n.parms <- NA
   summary$nobs    <- sum(data$freq)
   summary$method  <- "ChenLloyd"

   res$summary <- summary
   res$datetime <- Sys.time()

   class(res) <- "LP_CL_fit"  # chen lloyd fit
   res
}


