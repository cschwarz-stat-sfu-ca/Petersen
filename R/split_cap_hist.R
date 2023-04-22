#' Split a vector of capture histories into a matrix with one column for each occasion
#'

#' @param cap_hist A vector of capture histories.
#' @param sep What separates the individual history values
#' @param n  Number of sampling events in each history
#' @param prefix Prefix for labeling columns of matrix
#' @param make.numeric Change the expanded columns to numeric from character?

#' @details @template data.cap_hist
#'
#' @return A matrix of capture histories with 1 column per sampling event
#'
#' @importFrom stringr str_split_fixed fixed

#' @examples
#'
#' # standard 2 character capture histor
#' data(data_rodli)
#' Petersen::split_cap_hist(data_rodli$cap_hist)
#'
#' # history vector with ".." separating the fields
#' test <- c("1..1","1..0")
#' split_cap_hist(test, sep=stringr::fixed(".."))

#'
#' @export split_cap_hist
#'

split_cap_hist <- function(cap_hist, sep="", n=2, prefix="t", make.numeric=FALSE){
  # split a vector of capture histories
  if(!is.character(cap_hist))stop("Capture histories must be character strings")
  if(!is.vector(cap_hist))stop("Capture histores must be a vector")
  if(!is.character(sep))stop("'sep' argument must be character of length 1")
  if(length(sep)>1)stop("'sep' argument must be character of length 1")
  if(!is.numeric(n))stop("'n' argument must be numeric and positive and at least 2")
  if(n<2)stop("'n' argument must be numeric and positive and at least 2")

  req.nchar= n+ (n-1)*nchar(sep)
  if(!all(nchar(cap_hist)>=req.nchar))stop("All capture histories must have at least ", req.nchar, ' characters')

  if(sep == "") res <- stringr::str_split_fixed(cap_hist, pattern=sep       , n=n)
  if(sep != "") res <- stringr::str_split_fixed(cap_hist, pattern=stringr::fixed(sep), n=n)
  if(make.numeric){
    res <- matrix(as.numeric(res),    # Convert to numeric matrix
                  ncol = ncol(res))
  }

  colnames(res)<- paste0(prefix, 1:n)
  res
}

#split_cap_hist(rodli$cap_hist)

#test <- c("1..1","1..0")
#split_cap_hist(test, sep=stringr::fixed(".."))
