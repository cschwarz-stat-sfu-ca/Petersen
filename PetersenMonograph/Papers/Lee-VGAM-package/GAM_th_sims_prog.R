############################################################################

## R-code for fitting GAM and M_th (capture-recapture) models from 

## Stoklosa and Huggins (2012).

## The code below is needed for Section 6.5 of the simulation comparison 

## in Yee, Stoklosa, and Huggins. Note that functions are not explained in 

## full detail.

############################################################################

############################################################################



## Function below constructs a required B-spline basis function.

bspline <- function(x, n.dx, deg = 3, deriv = 0) {
  
  xl <- min(x) - .001*diff(range(x))
  
  xr <- max(x) + .001*diff(range(x))
  
  dx <- (xr - xl)/n.dx
  
  knots <- seq(xl - deg*dx, xr + deg*dx, by = dx)
  
  splineDesign(knots, x, ord = deg + 1, derivs = rep(deriv, length(x)))
  
}



## Function needed for capture-recapture model setup.

covaz.all <- function(cz, hist, n) {
  
  cze <- matrix(rep(cz, tau*rep(1, NCOL(cz)*n)), ncol = NCOL(cz))
  
  ct <- diag(1, tau)[, -tau]
  
  cte <- matrix(rep(t(ct), n), ncol = (tau - 1), byrow = TRUE)
  
  bhist <- matrix(as.numeric(t(apply(hist, 1, cumsum)) > 0), ncol = tau)
  
  bhist <- cbind(rep(0, n), bhist[, 1:(tau - 1)])
  
  cbe <- as.vector(t(bhist))
  
  cbind(cze, cte, cbe)
  
}



## Function needed for capture-recapture model setup.

absorb <- function(m, k, n = nrow(m)/k) {
  
  m.dot <- matrix(0, n, ncol(m))
  
  m.a <- array(t(m), dim = c(ncol(m), k, n))
  
  for (j in 1:n) {
    
    if(is.vector(m.a[, , j])){m.dot[j, ] <- sum(m.a[, , j])
    
    }
    
    else {m.dot[j, ] <- apply(m.a[, , j], 1, sum)}
  }
  
  m.dot
  
}



## Function needed choosing smoothing parameter lambda.

GAM_cond_th <- function(X.d1, X.d2, X.b, Y, tau, n, n.dx, lambda1, lambda2, nGCV, p.ord) {
  
  n.lambda1 <- lambda.all2(lambda1, nGCV)
  
  GCV.output.fun_th <- function(j) {
    
    GCV.fun_th(X.d1, X.d2, X.b, Y, tau, n, n.dx, n.lambda1[j], lambda2, p.ord)
  }
  
  GCV.output_th <- as.matrix(1:nGCV)
  
  GCV_th1 <- apply(GCV.output_th, 1, GCV.output.fun_th)
  
  return(list(GCV_th = GCV_th1))
  
}



## Function needed choosing smoothing parameter lambda.

lambda.all2 <- function(lambda, nGCV) {
  
  n.lambda <- NULL
  
  for(i in 1:nGCV){
    
    n.lambda <- c(n.lambda, lambda)
    
    lambda <- lambda*10
    
  }
  
  n.lambda
}



## Function needed choosing smoothing parameter lambda (via GCV).

GCV.fun_th <- function(X.d1, X.d2, X.b, Y, tau, n, n.dx, lambda1, lambda2, p.ord) {
  
  GAM_th(X.d1, X.d2, X.b, Y, tau, n, n.dx, lambda1, lambda2, p.ord)$GCV.score_th
  
}



## GAM fitting functon for M_th using penalized MLE.

GAM_th <- function(X.d1, X.d2 = NULL, X.b = NULL, Y, tau, n, n.dx = 10, lambda1, lambda2 = NULL, p.ord, kappa = 10^-6) {
  
  B1 <- bspline(X.d1, n.dx, deriv = 0)
  
  K1 <- ncol(B1)      # No. of B-splines.
  
  D1 <- diff(diag(K1), diff = p.ord)  # Differencing matrix for the penalty of order 2.
  
  P.l1 <- lambda1*crossprod(D1)

  B2 <- bspline(X.d2, n.dx, deriv = 0)
  
  K2 <- ncol(B2)    # No. of B-splines.
  
  D2.2 <- diff(diag(K2), diff = p.ord)  # Differencing matrix for the penalty of order 2. 
  
  P.l2 <- lambda2*crossprod(D2.2)
  
  q <- ncol(X.b)
  
  R <- diag(c(kappa*rep(1, (K1 + K2)), rep(0, q)))
  
  globlik_3b <- function(para) {
    
    delta1 <- para[1:K1]
    
    delta2 <- para[(K1 + 1):(K1 + K2)]
    
    beta <- para[((K1 + K2) + 1):length(para)]
    
    eta <- drop(B1%*%delta1) + drop(B2%*%delta2) + drop(X.b%*%beta)
    
    pij <- c(1/(1 + exp(-eta)))
    
    pijstar <- matrix(1 - pij, nrow = n, byrow = TRUE)
    
    prodpij <- apply(pijstar, 1, prod)
    
    llik <- sum(Y*eta + log(1 - pij)) - sum(log(1 - prodpij)) - 
      
      (1/2)*t(delta1)%*%P.l1%*%delta1 - (1/2)*t(delta2)%*%P.l2%*%delta2
    
    -llik
    
  }
  
  para2 <- c(rep(0,K1), rep(0,K2), rep(0,q))
  
  para <- para2
  
  e1 <- optim(para, globlik_3b, method = "BFGS", hessian = TRUE, control = list(maxit = 10000))
  
  log.Likelihood_th <- 2*e1$value
  
  delta1 <- e1$par[1:K1]
  
  delta2 <- e1$par[(K1 + 1):(K1 + K2)]
  
  beta <- e1$par[((K1 + K2) + 1):length(para)]
  
  eta <- drop(B1%*%delta1) + drop(B2%*%delta2) + drop(X.b%*%beta)
  
  delta <- c(delta1, delta2)
  
  pr <- c(1/(1 + exp(-eta)))
  
  pi.mat <- matrix(1 - pr, nrow = n, byrow = TRUE)
  
  pi <- 1 - apply(pi.mat, 1, prod)
  
  pi <- rep(pi, tau)
  
  mu <- pr/pi
  
  pi <- 1 - apply(pi.mat, 1, prod)
  
  V <- mu*(1 - pr) - mu^2*rep(1 - pi, each = tau)
  
  W <- as.numeric(V)
  
  W1 <- diag(W)
    
  Bz <- cbind(B1, B2, X.b)
  
  P <- rbind(cbind(P.l1, matrix(0, nrow(P.l1), K2)), cbind(matrix(0, nrow(P.l2), K2), P.l2))
  
  P <- rbind(cbind(P, matrix(0, nrow(P), q)), matrix(0, q, (q + ncol(P))))
  
  trace.lambda <- sum(diag(Bz%*%(ginv((t(Bz)%*%W1%*%Bz + P + R), (t(Bz)%*%W1)))))
  
  varB1 <- (ginv(t(Bz)%*%W1%*%Bz + P + R)%*%(t(Bz)%*%W1%*%Bz)%*%ginv(t(Bz)%*%W1%*%Bz + P + R))
    
  pi <- 1- apply(pi.mat, 1, prod)
  
  N.est <- sum(1/pi)
  
  pi2 <- (1 - pi)/(pi)^2
  
  Xall.dot <- absorb(pr*Bz, tau)
  
  d2 <- apply(pi2*Xall.dot, 2, sum)
  
  varN.est <- sum((1 - pi)/(pi)^2) + d2%*%varB1%*%d2
  
  AIC_th <- log.Likelihood_th + 2*(trace.lambda)
  
  n.1 <- length(Y)
  
  rss <- sum((Y - mu)^2)
  
  GCV.score_th <- n.1*rss/(n.1 - 1.4*trace.lambda)^2
  
  list(pr.est = pr, beta.est_th = beta ,delta.est_th = delta, N.est_th = N.est,
       
       varN.est_th = varN.est, log.Likelihood_th = log.Likelihood_th,
       
       AIC_th = AIC_th, GCV.score_th = GCV.score_th, varB1 = varB1,
       
       EDF_th = trace.lambda)
  
} 

