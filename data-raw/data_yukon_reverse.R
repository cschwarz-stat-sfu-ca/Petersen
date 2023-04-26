## code to prepare `data_yukon_reverse` dataset

## This is the data from the Hamazaki (2014) paper on estimating
## abundance using reverse mark-recapture data

## This is the data from the 2011 data in Table 2 of the above paper.

## Estimated that total escapement to Canada (plus  harvest)
## was 66,225 (SE 1574)
## Estimated that proportion of stock that was Canadian was .34644 (SE .030)
## We need to convert the latter into a sample size with similar SE

library(usethis)

phat <- .34644
phat.se <- .030

n_genetic <- trunc(phat*(1-phat)/phat.se^2)
m2  <- round(n_genetic * phat)
cat(n_genetic, m2, m2/n_genetic, "\n")

data_yukon_reverse <- data.frame(
  cap_hist=c("10",        "11",     "01"),
  freq    =c( 66225-m2,    m2,      n_genetic-m2),
  se      =c( 1574,       sqrt(n_genetic*phat*(1-phat)), 0),
  comment =c("Escapement + harvest - genetic samples recaptured cdn fish",
             "Genetic sample that are cdn fish",
             "Genetic sample that are non-cdn fish")
)

data_yukon_reverse

usethis::use_data(data_yukon_reverse, overwrite = TRUE)
