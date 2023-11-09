
#' Convert capture history data to n, m and u for use in BTSPAS
#'
#'
#' @template param.data
#' @param sep Separator used between strata in cap_hit
#' @importFrom utils head tail


#' @details
#'
#' The frequency variable (\code{freq} in the \code{data} argument) is the number of animals with the corresponding capture history.
#'
#' Capture histories (\code{cap_hist} in the \code{data} argument) are character values of the format
#' \code{xx..yy} is a capture_history where \code{xx} and \code{yy} are the temporal stratum
#' (e.g., julian week) and \code{'..'} separates
#' the two temporal strata.
#'  If a fish is released in temporal stratum and never captured again, then \code{yy} is set to 0;
#'  if a fish is newly captured in temporal stratum \code{yy}, then \code{xx} is set to zero.
#'  For example, a capture history of \code{23..23} indicates animals released in temporal stratum
#'  23 and recaptured in temporal stratum 23; a capture history of \code{23..00}
#'   indicates animals released in temporal stratum
#'  23 and never seen again; a capture history of \code{00..23}
#'   indicates animals newly captured in temporal stratum
#'  23 at the second sampling event.
#'
#'. In the diagonal case, fish are only recovered in the same temporal stratum.
#'  In the non-diagonal case, fish are allowed to move among temporal strata.
#'
#'  It is not necessary to label the temporal strata starting at 1; BTSPAS will treat the smallest
#'  value of the temporal strata seen as the first stratum and will interpolate for temporal strata
#'  without any data. Temporal strata labels should be numeric, i.e., do NOT use A, B, C etc.

#'
#
#'
#' @returns A list with entries for the stratum index, n (number released), m matrix
#' of recoveries in the current, next, etc stratum, and u (number of unmarked fish)
#' captured in this recovery stratum.
#'
#' @examples
#'
#' data(data_btspas_diag1)
#' cap_hist_to_n_m_u(data_btspas_diag1)
#'
#' data(data_btspas_nondiag1)
#' cap_hist_to_n_m_u(data_btspas_nondiag1)

#'
#'
#' @export cap_hist_to_n_m_u
#'
#'

cap_hist_to_n_m_u <- function(data, sep=".."){

  # to avoid CRAN check errors
  freq <- NULL
  ..ts1 <- NULL

  # some basic data checking
  check.cap_hist_temporal.df(data)

  # Extract the two strata
  data.aug <- cbind(data, split_cap_hist(data$cap_hist, sep=sep, prefix="..ts", make.numeric=TRUE))

  # check non-diagonal case, i.e. ts1<=ts2 unless 0
  if( any(data.aug$..ts1 > 0 & data.aug$..ts2 > 0   &  (data.aug$..ts1 > data.aug$..ts2 )))
      stop("Release stratum must be <= recovery stratum ")

  # get the summary statistics
  n1 <- plyr::ddply(data.aug, "..ts1", plyr::summarize, n1=sum(freq))
  n1 <- n1[ n1$..ts1 >0,]
  # add missing values but warn the user
  time <- data.frame(..ts1=min(n1$..ts1):max(n1$..ts1))
  n1 <- merge(n1, time, all=TRUE)
  if(any(is.na(n1))){
    warning("*** Caution... Missing value for n1 set to 0 ***\n")
  }
  n1[is.na(n1)]<- 0

  #browser()
  u2 <- plyr::ddply(data.aug, "..ts2", plyr::summarize, u2=sum(freq[ ..ts1==0]))
  u2 <- u2[ u2$..ts2>0 , ]
  time <- data.frame(..ts2=min(n1$..ts1):max(u2$..ts2))
  u2 <- merge(u2, time, all=TRUE)
  # notice that we leave is.na(u2) as missing to distinguish from 0
  if(any(is.na(u2))){
    warning("*** Caution... Missing value for u2. These are not set to zero ***\n")
  }


  m2 <- matrix(0, nrow=max(n1$..ts1)-min(n1$..ts1)+2,
                  ncol=max(u2$..ts2)-min(n1$..ts1)+2)
  rownames(m2)<- as.character(c(0,min(n1$..ts1):max(n1$..ts1)))
  colnames(m2)<- as.character(c(0,min(n1$..ts1):max(u2$..ts2)))
  data.aug$..ts1c <- as.character(data.aug$..ts1)
  data.aug$..ts2c <- as.character(data.aug$..ts2)
  m2[ as.matrix(data.aug[,c("..ts1c","..ts2c")])] <- data.aug$freq
  m2 <- m2[ -1, -1]# get rid of 0 row and 0 column

  # rotate m2 from upper diagonal to left block
  # add a column at the end of m2 with the amount of shifting
  m2 <- cbind(m2, shift=-(0:(nrow(m2)-1)))

  # https://predictivehacks.com/?all-tips=how-to-circular-shift-vectors-in-r
  custom_roll <- function( x , n ){
    if( n == 0 | n %% length(x)==0) {
      return(x)
    }

    else if (abs(n)>length(x)) {
      new_n<- (abs(n)%%length(x))*sign(n)
      return(c( utils::tail(x,new_n) , utils::head(x,-new_n) ))
    }
    else {
      return(c( utils::tail(x,n) , utils::head(x,-n) ))
    }
  }

  m2 <- plyr::aaply(m2, 1, function (x){
     x <- custom_roll(x[ -length(x)], x["shift"])
     x
  })

  # remove all columns of m2 that are all zero at the end
  select <- !rev(cumprod(rev(apply(m2==0,2,all))))
  m2 <- m2[,select, drop=FALSE]
  rownames(m2) <- NULL
  colnames(m2) <- NULL

  list(..ts=as.vector(time$..ts2),
       n1=as.vector(n1$n1),
       m2=m2,
       u2=as.vector(u2$u2))
}
