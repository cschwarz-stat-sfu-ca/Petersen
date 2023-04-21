#' Print the results from a fit a Lincoln-Petersen Model with incomplete stratification
#'
#' @param IS.results Results from fitting an incomplete stratification model


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
