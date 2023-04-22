#' @param  p_model Model for the captured probabilities. This can reference
#' other variables in the data frame, plus a special reserved term \code{..time} to indicate
#' a time dependence in the capture probabilities. For example, \code{p_model=~1} would indicate
#' that the capture probabilities are equal across the sampling events;
#' \code{p_model=~..time} would indicate that the capture probabilities vary by sampling events;
#' \code{p_model=~sex*..time} would indicate that the capture probabilities vary across
#' all combination of sampling events (\code{..time}) and a stratification variable (\code{sex}). The \code{sex} variable
#' also needs to be in the data frame.
#'
#' For some models (e.g., tag loss models), the \code{..time} variable cannot be used because
#' the models only have one capture probability (e.g., only for event 1).
#'

