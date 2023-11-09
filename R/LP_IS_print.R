#' Print the results from a fit a Lincoln-Petersen Model with incomplete stratification
#'
#' @param IS.results Results from fitting an incomplete stratification model
#' @returns
#' A nicely formatted report showing the results of the fit.
#' * Model information (summary of arguments, model name, negative log-likelihood, number of parameters, AICc)
#' * Raw data used in the fit (history, frequency,categories)
#' * Initial values used in the optimization of the likelihood for the parameters
#' * Design matrix and offset values for the parameters
#' * Maximum likelihood estimates for the parameters and estimated abundance by category
#' * SE for the above
#' * Observed and expected counts for each capture history
#' * Residual plot constructed from the previous observed and expected counts
#'


#' @examples
#'
#' data(data_wae_is_short)
#' res <- Petersen::LP_IS_fit(data=data_wae_is_short, p_model=~-1 + ..cat:..time)
#' LP_IS_print(res)
#'
#' @export LP_IS_print
#'

LP_IS_print <- function(IS.results){
  #  Print the results from an imcomplete stratification model
  # This is just a wrapper to the published code.

  if(!inherits(IS.results, "LP_IS_fit"))stop("Only works with results from incomplete stratification model")

  IS.print.output(IS.results$fit)
  invisible()
}
