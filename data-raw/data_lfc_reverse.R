## code to prepare `data_lfc_reverse` dataset

## This is the data provided by Kaitlyn Dionne, DFO.

## Arbeider et al (2020) proposed to estimate the run size of Lower Fraser River Coho (LFC) using a geographically stratified reverse-capture
## method. Briefly, a LFC coho swim upstream, they are sampled near New Westminister, BC, which is downstream from several major rivers upwhich
## are large spawning populations. These sample fish are assigned to the spawning population using genetic and other methods. These
## spawning populations are identified as the Chilliwack Hatchery (denoted *C*), the Lilloet River natural spawning population (deonted *L*),
## the Nicomen Slough population (denote *N*) and all other population (denoted as *0*). Notice that the sample fish at New Westminister
## are NOT physically tagged, and population assignment is through genetic and other measures.

## The upstream migration extends over two months (September to October) and is divided into 3 temporal strata corresponding
## to *Early* (denoted *1E*), *Peak* (denoted *2P*) and *Late* (denoted *3L*). The digits 1, 2, 3 in front of the codes ensures
## that the temporal strata are sorted temporally, but this is merely a convenience and does not affect the results.

## The spawning populations at *C*, *L*, and *N* are estimated by a variety of methods (see Arbeider, et al. 2020).
## Each of the population estimates also has an estimated (which will be ignored for now).

library(usethis)

data_lfc_reverse.csv <- textConnection("
cap_hist, freq, SE
C..1E,  23, NA
C..2P,  44, NA
C..3L,   6, NA
C..0, 66127, 13039
L..1E,  12, NA
L..2P,   8, NA
L..3L,   2, NA
L..0, 14300, 403
N..1E,   1, NA
N..2P,  13, NA
N..3L,   3, NA
N..0, 15274, 461
0..1E, 37,  NA
0..2P, 146, NA
0..3L,  46, NA")

data_lfc_reverse <- read.csv(data_lfc_reverse.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)

data_lfc_reverse

usethis::use_data(data_lfc_reverse, overwrite = TRUE)
