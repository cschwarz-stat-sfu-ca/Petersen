# sources all of the files needed for Partial stratification in two-sample capture-recapture experiments.
# path below should points to the Functions directory 
# setwd("U:\\Lasantha/Research/Rcode")
# setwd("C:\\Users/wpremara/Dropbox/Research/Rcode")
# setwd("C:\\Lasantha/SFU/Research/Rcode") # laptop
# load required libraries
library(numDeriv)# this is required to Calculate a numerical approximation to the Hessian matrix
library(plyr)
library(Matrix)
library(msm)
library(alabama)
library(knitr)
library(ggplot2)
library(Rsolnp)
library(directlabels)
library(gridExtra)
library(MASS)
library(matrixcalc)
library(MCMCpack) # this is needed for the dirichlet dstribution
library(coda) # to create a Markov Chain Monte Carlo object. 
library(R2WinBUGS) # this is needed for as.bugs.array

#Source all of the funtions in the Functions library 

files = list.files(path="./Functions", full.names=TRUE)
l_ply(files,source)

# source(file="AICc.R")
# source(file="get.data.R")
# source(file="create.DM.R")
# source(file="fit.model.R")
# source(file="expected.counts.R")
# source(file="helper.functions.R")
# source(file="initial.estimates.R")
# source(file="model.comparison.table.R")
# source(file="neg.log,likelihood.R")
# source(file="unpack.parm.R")
# source(file="pack.parm.R")
# source(file="print.output.R")
# source(file="unpack.parm.R")

