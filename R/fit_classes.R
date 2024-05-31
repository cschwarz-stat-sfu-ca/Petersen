#' **LP_fit**, **LP_IS_fit**, **LP_SPAS_cit**, **CL_fit**, **LP_BTSPAS_fit_Diag**, **LP_BTSPAS_fit_NonDiag**, **LP_CL_fit** classes.
#'
#' We assign a "class" (one of the classes above) to the results from one of the fitting methods. This class
#' designation is **only used to ensure that estimation routines have the correct type of fit when finding estimates of abundance**,
#' and **model averaging only considers models of comparable class** when creating the
#' model averaging results of abundance.
#'
#' The structure of an object of the above classes is roughly comparable across classes. The object should be a
#' list with the following objects
#' * **summary** A data frame with the model for the parameters of the model;
#' the conditional log-likelihood; the number of parameters; the number of parameters, and method used to fit the model
#' * **data** A data frame with the raw data used in the fit
#' * **fit** Results of the fit including the estimates, SE, vcov, etc.
#' * **datetime** Date and time the fit was done
#'
#' Other objects may also be included in the list.
#'
#' After a fit performed, estimates of **ABUNDANCE** are extracted using the **xxx_est()** function corresponding to the **xxx_fit()** function
#' used to estimate parameters. This separation occurs because
#'
#' - abundance is a derived estimate; the fits are based on conditional likelihood (on
#' the observed data) and the abundance parameter does not appear in the conditional likelihood. Abundance is usually
#' estimated using a variation of a Horvitz-Thompson estimator.
#' - you can obtain estimates of overall abundance, subsets of the population (e.g., sex or other strata) from the same fit, and
#' so you don't need to do several fits to get several estimates of abundance.

#' @return NOTHING. This is not a function, but only documents how I use "classes".
#' @export fit_classes
#'

fit_classes <- function(){}

