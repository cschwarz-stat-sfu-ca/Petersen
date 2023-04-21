# Petersen

Use of the Petersen Capture-Recapture estimates in fisheries management.

## Versions and installation

  * **CRAN**  Download the **Petersen** package

  * **Github** To install the latest development version from Github, 
    install the newest version of the **devtools** package; then run
```
devtools::install_github("cschwarz-stat-sfu-ca/Petersen", dependencies = TRUE,
                        build_vignettes = TRUE)
```

## Features

The Petersen-method is the simplest of more general capture-recapture
methods which are extensively reviewed in Williams et al. (2002).
Despite the Petersen method's simplicity, many of the properties of the estimator, 
and the effects of violations of assumptions are similar to these more complex capture-recapture
studies. Consequently, a firm understanding of the basic principles
learned from studying this method are extremely useful to develop an
intuitive understanding of the larger class of studies.

The purpose of this R package is to bring together a wide body of older
and newer literature on the design and analysis of the "simple"
two-sample capture-recapture study. This monograph builds upon the
comprehensive summaries found in Ricker (1975), Seber (1982), and
William et AL (2002), and incorporates newer works that have not yet
summarized. While the primary emphasis is on the application to
fisheries management, the methods are directly applicable to many other studies.

The core of the package is the use of conditional likelihood estimation that 
allows for covariates which are not observed on animals not handled.

The packages includes functions for the analysis of

- simple studies with no covariates or strata
- simple stratified-Petersen studies or fixed continuous covariates
- incompletely stratified studies where only a sub-sample is stratified to save costs
- geographically stratified studies (wrapper to *SPAS*)
- temporally-stratified studies (wrapper to *BTSPAS*)
- double tagging studies including reward tagging studies
- multiple-Petersen studies (call *RMark/MARK* and use mark-resight methods therein)
- forward and reverse-Petersen studies and their combination

## References

Petersen, C. G. J. (1896). The Yearly Immigration of Young Plaice into
the Limfjord from the German Sea, Etc. Report Danish Biological Station
6, 1--48.

Seber, G. A. F. (1982). The Estimation of Animal Abundance and Related
Parameters. 2nd ed. London: Griffin.


Williams, B. K., J. D. Nichols, and M. J. Conroy. (2002). Analysis and
Management of Animal Populations. New York: Academic Press.

