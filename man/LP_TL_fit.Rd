% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LP_TL_fit.R
\name{LP_TL_fit}
\alias{LP_TL_fit}
\title{Fit a Lincoln-Petersen Model with Tag Loss using conditional likelihood}
\usage{
LP_TL_fit(
  data,
  dt_type = NULL,
  p_model,
  rho_model,
  all_beta.start = NULL,
  trace = FALSE
)
}
\arguments{
\item{data}{Data frame containing the variables:
\itemize{
\item \strong{cap_hist} Capture history (see details below)
\item \strong{freq} Number of times this capture history was observed
}
plus any other covariates (e.g. discrete strata and/or continuous covariates) to be used
in the model fitting.}

\item{dt_type}{Double Tag type. Valid values are
\code{notD}, \code{twoD}, and \code{t2perm}
for two indistinguishable tags; two distinguishable tags, when the second tag is a permanent tag and cannot be lost,
respectively.}

\item{p_model}{Model for the captured probabilities. This can reference
other variables in the data frame, plus a special reserved term \code{..time} to indicate
a time dependence in the capture probabilities. For example, \code{p_model=~1} would indicate
that the capture probabilities are equal across the sampling events;
\code{p_model=~..time} would indicate that the capture probabilities vary by sampling events;
\code{p_model=~sex*..time} would indicate that the capture probabilities vary across
all combination of sampling events (\code{..time}) and a stratification variable (\code{sex}). The \code{sex} variable
also needs to be in the data frame.

For some models (e.g., tag loss models), the \code{..time} variable cannot be used because
the models only have one capture probability (e.g., only for event 1).}

\item{rho_model}{Model for retention probabilities}

\item{all_beta.start}{Initial values for call to optimization routine for the beta parameters (on the logit scale).
The values will be replicated to match
the number of initial beta parameters needed. Some care is needed here since the parameter order are for the p1 probabilities
and then for the rho probabilities}

\item{trace}{If trace flag is set in call when estimating functions}
}
\value{
An list object with the fit information
}
\description{
This will take a data frame of capture histories, frequencies, and
additional covariates (e.g., strata and/or continuous covariates) and the model
for p1 and the tag retention probabilities
and will use conditional likelihood (conditional on capture at time 2)
to fit the model. The population abundance is estimated using
a Horvitz-Thompson type estimator and the user can request abundance
estimates for sub-sets of the population.
}
\details{
The frequency variable (\code{freq} in the \code{data} argument) is the number of animals with the corresponding capture history.

Capture histories (\code{cap_hist} in the \code{data} argument) are character values of length 4.

If the tag loss model is two indistinguishable tags (\code{dt_type="notD"}), then valid capture histories are:
\itemize{
\item \strong{1100}  Animals double tagged but never seen again.
\item \strong{111X}  Animals double tagged, but only 1 tag was present when animal recaptured at second event.
\item \strong{1111}  Animals double tagged and both tags present when animal recaptured at second event.
\item \strong{1000}  Animals single tagged and never seen again.
\item \strong{0100}  Animals single tagged and never seen again.
\item \strong{1010}  Animals single tagged and recaptured with the single tag.
\item \strong{0101}  Animals single tagged and recaptured with the single tag.
\item \strong{0010}  Animals APPARENTLY captured for the first time at event 2. This includes animals that are
newly captured, plus fish that were tagged and lost all their tags, and were captured again
}

If the tag loss model is two distinguishable tags (\code{dt_type="twoD"}), then valid capture histories are the same
as above except the history \code{111X} is replaced by:
\itemize{
\item \strong{1110} Animals double tagged, but only the first of the double tags applied  was present when animal recaptured at event 2,
\item \strong{1101} Animals double tagged, but only the second of the double tags applied was present when animal recaptured at event 2.
}

If the second tag is a permanent batch mark (\code{dt_type="t2perm"}), then valid capture histories are:
\itemize{
\item \strong{1P00}  Animals double tagged but never seen again.
\item \strong{1P0P}  Animals double tagged,but non-permanent tag missing when animal recaptured at second event.
\item \strong{1P1P}  Animals double tagged and both tags present when animal recaptured at second event.
\item \strong{1000}  Animals single tagged and never seen again.
\item \strong{0P00}  Animals single tagged with a permanent batch mark only and never seen again.
\item \strong{1010}  Animals single tagged and recaptured with the single tag.
\item \strong{0P0P}  Animals single tagged with the permanent batch mark and recaptured with the permanent tag.
\item \strong{0010}  Animals APPARENTLY captured for the first time at event 2. This includes animals that are
newly captured, plus fish that were tagged and lost all their tags, and were captured again
}
}
\examples{

data(data_kokanee_tagloss)
Petersen::LP_TL_fit(data=data_kokanee_tagloss, p_model=~1, rho_model=~1, dt_type="notD")

}
\author{
Schwarz, C. J. \email{cschwarz.stat.sfu.ca@gmail.com}.
}
