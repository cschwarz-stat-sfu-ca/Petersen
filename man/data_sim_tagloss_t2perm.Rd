% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_sim_tagloss_t2perm}
\alias{data_sim_tagloss_t2perm}
\title{Simulated data for tag loss with second permanent tag.}
\format{
\subsection{\code{data_sim_tagloss_t2perm}}{

A data frame with many rows and 2 columns

\describe{
\item{\code{cap_hist}.}{Capture history (1000, 1010, 1P00, 1P0P, 1P1P, 0010).}
\item{\code{freq}.}{Number of fish with this history.}
}
}
}
\usage{
data(data_sim_tagloss_t2perm)
}
\description{
This is simulated data with the parameter values given in details.
}
\details{
\preformatted{data_sim_tagloss_t2perm <-LPTL_simulate(
      dt_type="t2perm",         # second permanent
      N=10000,
      cov1=function(N)         {rep(1,N)},
      cov2=function(cov1)      {rep(1,  length(cov1))},
      p1  =function(cov1, cov2){rep(.1, length(cov1))},
      pST =function(cov1, cov2){rep(.25,length(cov1))},
      rho1=function(cov1, cov2){rep(.70,length(cov1))},
      rho2=function(cov1, cov2){rep(1,  length(cov1))}, # permanent tag
      p2  =function(cov1, cov2){rep(.1, length(cov1))},
      seed=234523, trace=FALSE)}
}
\keyword{datasets}
