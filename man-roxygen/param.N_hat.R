#' @param  N_hat A formula requesting which abundance estimates should be formed. The formula are
#' expanded against the data frame to determine which records form part of the abundance estimate.
#' The formula is evaluated against the \code{data} frame used in the fit using the \code{model.matrix()} function,
#' and each column of the model
#' matrix is used to form an estimate.
#'
#' Some familiarity on how \code{model.matrix()} generates the model matrix of coefficients used in the expansion
#' is needed.
#' For example \code{N_hat=~1} creates a model matrix with 1 column (representing the intercept) and
#' so requests abundance over the entire population;
#' Specifying \code{N_hat=~-1+Sex} creates a model matrix with 2 columns (one for each sex) consisting of 0/1 depending
#' if that row of the data frame is M/F. Hence, two abundance estimates (one for each sex) is computed.
#' On the other hand, \code{N_hat=Sex} generates a model matrix where the first column is all 1's, and
#' a second column which is 0/1 depending if the row in the data frame is the "second" sex. Hence, this will
#' request the overall abundance (over both sexes) and the estimate of abundance for the second sex.
#'
#' In addition to the variables in the \code{data} frame, special variables include \code{..EF} to allow access to the expansion
#' factor so you can request a "truncated" Horvitz-Thompson estimator using \code{N_hat=~-1+I(as.numeric(..EF<1000))}
#' to only use those animals with expansion factors less than 1000 in forming the estimate.
