
#' Estimating abundance of outgoing smolt - BTSPAS - diagonal case
#'
#' This is the first diagonal case dataset from BTSPAS.
#'
#' Consider an experiment to estimate the number of outgoing smolts on a small river. The
#' run of smolts extends over several weeks. As smolts migrate, they are captured and marked
#' with individually numbered tags and released at the first capture location using, for example, a
#' fishwheel. The migration continues, and a second fishwheel takes a second sample several
#' kilometers down stream. At the second fishwheel, the captures consist of a mixture of marked
#' (from the first fishwheel) and unmarked fish.
#'
#' The efficiency of the fishwheels varies over time in response to stream flow, run size passing
#' the wheel and other uncontrollable events. So it is unlikely that the capture probabilities are
#' equal over time at either location, i.e. are heterogeneous over time.
#'
#' We suppose that we can temporally stratify the data into, for example, weeks, where the
#' capture-probabilities are (mostly) homogeneous at each wheel in each week. Furthermore, suppose that
#' fish captured and marked in each week tend to migrate together so that they are
#' captured in a single subsequent stratum. For example,
#' suppose that in each julian week \eqn{j}{j}, \eqn{n1_j}{n1_j} fish are marked and released above the rotary screw trap.
#' Of these, \eqn{m2_j}{m2_j} are recaptured. All recaptures take place in the week of release,
#' i.e. the matrix of releases and recoveries is diagonal.
#' The \eqn{n1_j}{n1_j} and \eqn{m2_j}{m2_j} establish the capture efficiency of the second trap in julian week \eqn{j}{j}.
#'
#' At the same time, \eqn{u2_j}{u2_j} unmarked fish are captured at the screw trap.
#'
#' Capture-efficiency may be related to flow, so the log(flow) is also recorded.

#'
#' @format ## `data_btspas_diag1`
#' A data frame with many rows and 3 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history of the form `jweek..jweek' for fish that are recaptured
#' in the same julian week; '0..jweek' for unmarked fish newly captured in that julian week ; 'jweek..0' for fish
#' released in the julian week but never recaptured.}
#' \item{\code{freq}.}{Number of fish with this history.}
#' \item{\code{logflow}}{log(flow) for this julian week}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_btspas_diag1
#' @usage data(data_btspas_diag1)
"data_btspas_diag1"


#' Estimating abundance of salmon - BTSPAS - non-diagonal case
#'
#' This is the first non-diagonal case dataset from BTSPAS.
#'
#' Incoming sockeye salmon are captured on a first wheel, tagged with color tags
#' that vary by week, and recaptured on an upriver weir.
#' The upriver weir was not in operation for the first few weeks.
#
#'
#' @format ## `data_btspas_nondiag1`
#' A data frame with many rows and 3 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history of the form `week1..week2' for fish that are released on
#' week 1 and recaptured
#' on week 2 ; '0..week22' for unmarked fish newly captured in week 2; 'week1..0' for fish
#' released in week 1 but never recaptured.}
#' \item{\code{freq}.}{Number of fish with this history.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_btspas_nondiag1
#' @usage data(data_btspas_nondiag1)
"data_btspas_nondiag1"


#' Capture-recapture on Kokanee in Metolius River with tag loss
#'
#'  This is the data from  Hyun et al (2012).
#'  In August and September 2007, the period just before the spawning run,
#'  adult kokanee were collected by beach seining in the upper arm of the
#'  lake near the confluence with the Metolius River.
#'  Fish were tagged with nonpermanent, plastic T-bar anchor tags
#'   and then were released back into the lake.
#'    Randomly selected fish received single tags of one color,
#'    while the other fish received two tags of a second color (i.e., the double tags were identical in color).
#'  In late September through October, spawning ground surveys were conducted by 2–3 people walking
#'  abreast in a downstream direction (or floating, in sections where the water depth and flow were too great to allow walking).
#'  Instead of being physically recaptured, the fish were resighted as they swam freely in the clear,
#'  relatively shallow water within the spawning areas of the river.
#'  The total number of fish observed with or without a tag (or tags) was recorded for each section,
#'  and information on the number and color of tags for each marked fish was also noted.
#'
#'  Because fish were not handled, it is not possible to know which of the double tags were lost, and so only
#'  models with equal retention probabilities and non-distinguishable double tags should be fit. Note
#'  that the capture history for a lost of 1 indistinguishable tag is 111X rather than 1110 or 1101 (both of
#'  which are not allowed in the model with indistinguishable double tags)
#'
#' @format ## `data_kokanee_tagloss`
#' A data frame with many rows and 2 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history (1000, 1010, 1100, 1110, 1111).}
#' \item{\code{freq}.}{Number of fish with this history. Always 1}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_kokanee_tagloss
#' @usage data(data_kokanee_tagloss)
"data_kokanee_tagloss"


#' Capture-recapture experiment on Northern Pike in Mille Lacs, MN, in 2005.
#'
#'  Fish were tagged on the spawning grounds and recovered in the summer gillnet assessment.
#'  Fish were double tagged, and a tag loss analysis
#'  showed that tag loss was negligible. It will be ignored here.
#'  Length was measured a both times and didn't not change
#'  very much between the two sampling occasions. The value
#'  recorded below is the average of the two lengths if both
#'  lengths were present.
#' Fish that were not sexed or measured for length are ignored and not included
#'
#' @format ## `data_NorthernPike`
#' A data frame with many rows and 4 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history (10, 01, or 11). Note that
#' \eqn{n_{10}=n_1-m_2}{n10=n1-m2}; \eqn{n_{01}=n_2-m_2}{n01=n2-m2}; and \eqn{n_{11}=m2}{n11=m2}}
#' \item{\code{freq}.}{Number of fish with this history. Always 1}
#' \item{\code{Sex}.}{Sex of the fish. M=Male; F=Female}
#' \item{\code{Length}}{Length of the fish in inches.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_NorthernPike
#' @usage data(data_NorthernPike)
"data_NorthernPike"


#' Capture-recapture experiment on Northern Pike in Mille Lacs, MN, in 2005 with tagloss information.
#'
#'  Fish were tagged on the spawning grounds and recovered in the summer gillnet assessment.
#'  Fish were double tagged and the double tagging information is included here.
#'  Length was measured a both times and didn't not change
#'  very much between the two sampling occasions. The value
#'  recorded below is the average of the two lengths if both
#'  lengths were present.
#' Fish that were not sexed or measured for length are ignored and not included
#'
#' @format ## `data_NorthernPike_tagloss`
#' A data frame with many rows and 4 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history (10, 01, or 11). Note that
#' \eqn{n_{10}=n_1-m_2}{n10=n1-m2}; \eqn{n_{01}=n_2-m_2}{n01=n2-m2}; and \eqn{n_{11}=m2}{n11=m2}}
#' \item{\code{freq}.}{Number of fish with this history. Always 1}
#' \item{\code{Sex}.}{Sex of the fish. M=Male; F=Female}
#' \item{\code{Length}}{Length of the fish in inches.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_NorthernPike_tagloss
#' @usage data(data_NorthernPike_tagloss)
"data_NorthernPike_tagloss"


#' Capture-recapture experiment at Rodli Tarn.
#'
#' Ricker (1975) gives an example of work by Knut Dahl on estimating the number of brown trout
#' (\emph{Salmo truitta}) in some small Norwegian tarns. Between 100 and 200 trout were caught by
#' seining, marked by removing a fin (an example of a batch mark) and distributed
#' in a systematic fashion around the tarn to encourage mixing.
#' A total of \eqn{n_1}{n1}=109 fish were
#' captured, clipped and released,
#' \eqn{n_2}{n2}=177 fish were captured at the second occasion, and \eqn{m_2}{m2}=57
#' marked fish were recovered.
#'
#' @format ## `data_rodli`
#' A data frame with 3 rows and 2 columns:
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history (10, 01, or 11). Note that
#' \eqn{n_{10}=n_1-m_2}{n10=n1-m2}; \eqn{n_{01}=n_2-m_2}{n01=n2-m2}; and \eqn{n_{11}=m_2}{n11=m2}}
#' \item{\code{freq}.}{Number of fish with this history}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_rodli
#' @usage data(data_rodli)
"data_rodli"


#' Simulated data for tag loss with 2 distinguishable tags.
#'
#'  This is simulated data with the parameter values given in the description.
#'
#'  \preformatted{data_sim_tagloss_twoD <-LPTL_simulate(
#'       dt_type="twoD",         # two distinguishable tags
#'       N=10000,
#'       cov1=function(N)         {rep(1,N)},
#'       cov2=function(cov1)      {rep(1,  length(cov1))},
#'       p1  =function(cov1, cov2){rep(.1, length(cov1))},
#'       pST =function(cov1, cov2){rep(.25,length(cov1))},
#'       rho1=function(cov1, cov2){rep(.70,length(cov1))},
#'       rho2=function(cov1, cov2){rep(.80,length(cov1))},
#'       p2  =function(cov1, cov2){rep(.1, length(cov1))},
#'       seed=234523, trace=FALSE)}

#'
#' @format ## `data_sim_tagloss_twoD`
#' A data frame with many rows and 2 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history (1000, 1010, 1100, 1110, 1101, 1111).}
#' \item{\code{freq}.}{Number of fish with this history.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_sim_tagloss_twoD
#' @usage data(data_sim_tagloss_twoD)
"data_sim_tagloss_twoD"


#' Simulated data for tag loss with second permanent tag.
#'
#'  This is simulated data with the parameter values given in details.
#'
#'  \preformatted{data_sim_tagloss_t2perm <-LPTL_simulate(
#'       dt_type="t2perm",         # second permanent
#'       N=10000,
#'       cov1=function(N)         {rep(1,N)},
#'       cov2=function(cov1)      {rep(1,  length(cov1))},
#'       p1  =function(cov1, cov2){rep(.1, length(cov1))},
#'       pST =function(cov1, cov2){rep(.25,length(cov1))},
#'       rho1=function(cov1, cov2){rep(.70,length(cov1))},
#'       rho2=function(cov1, cov2){rep(1,  length(cov1))}, # permanent tag
#'       p2  =function(cov1, cov2){rep(.1, length(cov1))},
#'       seed=234523, trace=FALSE)}

#'
#' @format ## `data_sim_tagloss_t2perm`
#' A data frame with many rows and 2 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history (1000, 1010, 1P00, 1P0P, 1P1P, 0010).}
#' \item{\code{freq}.}{Number of fish with this history.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_sim_tagloss_t2perm
#' @usage data(data_sim_tagloss_t2perm)
"data_sim_tagloss_t2perm"


#' Simulated data for reward tags used to estimate reporting rate
#'
#'  This is simulated data with the parameter values given in details.
#'
#'  \preformatted{
#'  data_sim_reward <-LP_TL_simulate(
#'       dt_type=dt_type,  #  permanent tag
#'       N=10000,
#'       cov1=function(N)         {rep(1,N)},
#'       cov2=function(cov1)      {rep(1,  length(cov1))},
#'       p1  =function(cov1, cov2){rep(.1, length(cov1))},
#'       pST =function(cov1, cov2){rep(.75,length(cov1))},
#'       rho1=function(cov1, cov2){rep(.70,length(cov1))},
#'       rho2=function(cov1, cov2){rep(1,  length(cov1))},  # permanent second tag
#'       p2  =function(cov1, cov2){rep(.1, length(cov1))},
#'       seed=45985, trace=FALSE)
#' # we don't have fish with both tags
#' data_sim_reward$cap_hist <- gsub("1P", "0P", data_sim_reward$cap_hist)
#' }
#'
#' @format ## `data_sim_reward`
#' A data frame with many rows and 2 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history (1000, 1010, 0P00, 0P0P, 0010).}
#' \item{\code{freq}.}{Number of fish with this history.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_sim_reward
#' @usage data(data_sim_reward)
"data_sim_reward"





#' Estimating abundance of salmon - SPAS - Harrison River
#'
#' Incoming sockeye salmon are captured on a first wheel, tagged with color tags
#' that vary by week, and recaptured on several spawning areas.
#'
#' @format ## `data_spas harrison`
#' A data frame with many rows and 2 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history of the form `week..area' for fish that are released on
#' week and recaptured
#' area ; '0..area' for unmarked fish newly captured in area; 'week..0' for fish
#' released in week  but never recaptured.}
#' \item{\code{freq}.}{Number of fish with this history.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_spas_harrison
#' @usage data(data_spas_harrison)
"data_spas_harrison"


#'  Walleye data with incomplete stratification with length covariate
#'
#'  Data used in
#'     Premarathna, W.A.L., Schwarz, C.J., Jones, T.S. (2018)
#'     Partial stratification in two-sample capture–recapture experiments.
#'     Environmetrics, 29:e2498. https://doi.org/10.1002/env.2498
#'
#'  Fish were tagged on the spawning grounds and
#'  recovered in the summer gillnet assessment.
#'
#'  Length was measured a both times and didn't not change
#'  very much between the two sampling occasions. The value
#'  recorded below is the average of the two lengths if both
#'  lengths were present.
#'
#'  Rather than sexing all of the fish, only a sub-sample of unmarked fish
#'  is sexed at each sampling occasion. Possible capture histories are then
#'    M0, F0, MM, FF, U0, UU, 0M, 0F

#' @format `data_wae_is_long`
#' A data frame with many rows and 3 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history with possible histories as noted below}
#' \item{\code{freq}.}{Number of fish with this history.}
#' \item{\code{length}}{Length of fish (inches)}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_wae_is_long
#' @usage data(data_wae_is_long)
"data_wae_is_long"


#'  Walleye data with incomplete stratification with no covariates and condensed
#'
#'  Data used in
#'     Premarathna, W.A.L., Schwarz, C.J., Jones, T.S. (2018)
#'     Partial stratification in two-sample capture–recapture experiments.
#'     Environmetrics, 29:e2498. https://doi.org/10.1002/env.2498
#'
#'  Data is slightly different from that in paper above because some fish did not
#'  have length measured and so were drop from data_wae_is_long and this is the
#'  condensed version of data_wae_is_long.
#'
#'  Fish were tagged on the spawning grounds and
#'  recovered in the summer gillnet assessment.
#'
#'  Rather than sexing all of the fish, only a sub-sample of unmarked fish
#'  is sexed at each sampling occasion. Possible capture histories are then
#'    M0, F0, MM, FF, U0, UU, 0M, 0F

#' @format `data_wae_is_short`
#' A data frame with many rows and 2 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history with possible histories as noted below}
#' \item{\code{freq}.}{Number of fish with this history.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_wae_is_short
#' @usage data(data_wae_is_short)
"data_wae_is_short"


#'  Yukon River data used for Reverse Capture-Recapture example.
#'
#'  Data from
#'     Hamazaki, T. and DeCovich, N. (2014).
#'     Application of the Genetic Mark–Recapture Technique for Run Size Estimation of Yukon River Chinook Salmon.
#'     North American Journal of Fisheries Management, 34, 276-286.
#'     DOI: 10.1080/02755947.2013.869283

#'  This is the data from the 2011 data in Table 2 of the above paper.

#' Estimated that total escapement to Canada (plus  harvest)
#' was 66,225 (SE 1574)
#' Estimated that proportion of stock that was Canadian was .34644 (SE .030)
#' We converted this into a "sample size" and number of fish with Cdn genetics
#' that gave the same SE.

#' @format `data_yukon_reverse`
#' A data frame with many rows and 4 columns
#'
#' \describe{
#' \item{\code{cap_hist}.}{Capture history with possible histories as noted below}
#' \item{\code{freq}.}{Number of fish with this history.}
#' \item{\code{SE}.}{SE of the number of fish with this history}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_yukon_reverse
#' @usage data(data_yukon_reverse)
"data_yukon_reverse"
