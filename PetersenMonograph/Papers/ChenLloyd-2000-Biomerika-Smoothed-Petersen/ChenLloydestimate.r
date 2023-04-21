# These functions were supplied by Lloyd, C. and Chen, S. in regards to their
# paper in Biometrika. Received on 2005-12-03
#
# Documentation added by C.Schwarz, 2005-12-03

"ChenLloydestimate" <- 
function(data, h1, h2, plt = T, w = F, pdf=F, covariatename="Covariate", title=" ")
{
#       Estimate the effect of heterogeneity upon the petersen estimator
#       Data consists of a list "data" with elements
#           data$MM - number of marks placed out in the covariate class with respective center in data$centers
#           data$mm - number of recaptures from this covariate class with respective center in data$centers
#           data$uu - number of UNmarked   from this covariate class
#           data$centre - centre of this covariate class
#       An easy way to create this data from a list of capture histories and their covariates is
#       given by the "bin"  function
#
#       I'm not sure what the "w" parameter measures, but it appears to be some sort of smoothing using
#       a weighted average.
#
#       The smoothing parameters are h1 and h2
#           h1 = std dev of the smoothing kernel for the recapture rates
#           h2 = std dev of the smoothing kernel for the 

  S <- function(x, h, xx = x, wgt = 0. * x + 1.){
#       Create the kernal density weights using a Normal distribution with std dev=h
	      k <- length(x)
	      M <- matrix(rep(xx, k), ncol = k) - matrix(rep(x, length(xx)), ncol = k, byrow = T)
	      M <- matrix(dnorm(M, sd = h), ncol = k) %*% diag(wgt)
	      sweep(M, 1., apply(M, 1., sum), "/")
  }

	MM <- as.vector(data$M)
	mm <- as.vector(data$m)
	uu <- as.vector(data$u)
	x <- data$centre

#       Get the smoothing matrix based on a normal kernal with bandwidth = std dev = h1
	S1 <- S(x, h1)
#       Smooth the marks put out and the marks recaptured to estimate the recapture rates (smoothed)
	Mtl <- as.vector(S1 %*% MM)
	mtl <- as.vector(S1 %*% mm)
	S2 <- S(x, h2, wgt = w/mtl + (1. - w))
	ptl <- mtl/Mtl     # smoothed recapture rates

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
	value$p <- round(1000. * ptl)/1000.
	value$N1 <- round(100. * N1)/100.
	value$N2 <- round(100. * N2)/100.
	value$sd1j <- round(sqrt(diag(VLN1)) * N1 * 100.)/100.
	value$sd2j <- round(sqrt(diag(VLN2)) * N2 * 100.)/100.
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
	value$total.est <- c(round(sum(N1)), round(sum(N2)), round(Nhom))
	value$total.sd <- round(10. * c(s1, s2, sdhom))/10.

	#    OPTIONAL GRAPHICAL OUTPUT
	if(plt == T) {
               if(pdf ==T){pdf("ChenLLoydestimateplot.pdf")}
#              if(pdf ==T){postscript("ChenLLoydestimateplot.ps")}

                par(mfcol=c(2,1))
		plot(x, mm/MM, xlab = covariatename, ylab = "weight")
		lines(x, ptl)
		title(paste(title,
                     " \n Estimated weight function - the prob of capture on the second occasion"))
		plot(x, N1, xlab = covariatename, ylab = "frequency")  # N1 is pop size at second sample
                lines(x, value$sd1j, lty=2.)   # standard error of estimate at each time point for ragged estimate
                lines(x, value$sd2j, lty=3.)   # standard error of smoothed estimate
		lines(x, N2)
                text(30,15000,labels=paste("Pop (SE) \n", 
                    value$total.est[1]," ( ", value$total.sd[1]," )"), pos=4)
		title("Estimated frequency distribution - estimated population size ")
                title(sub="SE for smoothed and raw estimates are shown in dashed lines at bottom")
                if(pdf==T){dev.off()}
	}
	value
}

"bin" <- function(x, pts)
{
	# Function bins a vector of numbers into bins with
	# centre in the vector pts. It returns a vector of
	# frequencies of the same dimension as pts.
        # Modified by C.schwarz, 2005-12-02 to use -1000*min and 1000*max as the upper and lower bound of first and last class
	s.p <- sort(pts)
	brks <- (c(-10000.*min(s.p), s.p) + c(s.p, 10000.*max(s.p)))/2.
	hist(x, breaks = brks, plot = F)$counts
}

"plugin.F" <- 
function(x, B = min(length(x), 49.), prin = F)
{
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
	grid <- function(x, n)
	{
		min(x, na.rm = T) + (-1.:(n + 1.))/n * (max(x, na.rm = 
			T) - min(x, na.rm = T))
	}
	xx <- grid(x, B - 3.)
	X <- matrix(rep(xx, rep(B, B)), ncol = B)
	X <- X - t(X)
	counts <- hist(x, breaks = grid(x, B - 2.), plot = F)$counts
	k <- matrix(rep(counts, rep(B, B)), ncol = B)
	k <- k * #
	t(k)
	qs <- quantile(x, na.rm = T)
	s <- (qs[4.] - qs[2.])/1.3400000000000001
	names(s) <- NULL
	h4 <- (1.2406999999999999 * s)/n^(1./7.)
	XX <- X/h4
	psi4 <- sum(k * (XX^4. - 6. * XX^2. + 3.) * dnorm(XX))/(n *
		n * h4^5.)
	h2 <- 0.79790000000000005/(psi4 * n)^(1./5.)
	XX <- X/h2
	psi2 <- sum(k * (XX^2. - 1.) * dnorm(XX))/(n * n * h2^3.)
	if(prin) {
		print(paste("Psi4= ", signif(psi4, 3.), "Psi2=", signif(
			psi2, 3.)))
		print(paste("Bandwidths used: h4=", signif(h4, 3.),
			"h2=", signif(h2, 3.)))
	}
	h2stage <- (-0.56000000000000005/psi2/n)^(1./3.)
	h0stage <- (1.583 * s)/n^(1./3.)
	if(h2stage <= 0.) {
		h2stage <- h0stage
	}
	c(h0stage, h2stage)
}

source("estimate.r")

library(Petersen)
data(NorthernPike)

range(NorthernPike$length)

#  extract the length vectors for fish captured in stratum1, stratum3, and both strata
   lengths <- NULL
   lengths$all <- NorthernPike$length
   lengths$x1  <- lengths$all[ substr(NorthernPike$cap_hist,1,1) == '1']
   lengths$x2  <- lengths$all[ substr(NorthernPike$cap_hist,2,2) == '1']
   lengths$x12 <- lengths$all[ NorthernPike$cap_hist == '11']

#  remove all  missing values as they cannot be used.
 
   data <- NULL
   data$x1 <- lengths$x1 [ !is.na(lengths$x1) & lengths$x1 > 17.99 ]     # select only those fish at least 18 inches in
   data$x2 <- lengths$x2 [ !is.na(lengths$x2) & lengths$x2 > 17.99 ]
   data$x11 <- lengths$x12 [ !is.na(lengths$x12) & lengths$x12 > 17.99]
   data
   print("Estimates for SP-SG phases")
   estimate(data, brief=F)

   
#  The SP-SG phases --------------------------------------------
#  create the number of marked, unmarked, and recaptured in 1 inch bins

   data$centres <- .5 + 16:45
   data$M      <- bin(lengths$x1 [ lengths$x1 > 17.99], data$centres)   # fish tagged
   data$m      <- bin(lengths$x12[ lengths$x12> 17.99], data$centres)   # number of marks recaptured
   data$u      <- bin(lengths$x2 [ lengths$x2 > 17.99], data$centres) - data$m # unmarked at time=2
   data
#  estimate the approximate optimal band widths
 
   plugin.F( lengths$x1)
   plugin.F( lengths$x2)
   plugin.F( lengths$x12)

#  This shows for this problem that optimal bandwidth is about 1 inch which happens to match our increments

#  Now to estimate the population size and get plots etc.
   print("Results from SP-SG phases")
   title<- "Lloyd and Chen (2000) estimator for population size SP-SG"
   ChenLloydestimate(data, 1, 1, covariatename="Length (in)", pdf=TRUE, title=title)


