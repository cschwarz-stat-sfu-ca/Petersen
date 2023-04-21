## code to prepare `kokanee_tagloss` dataset goes here
## This is the data in Hyan et al (2012)

data.csv <- textConnection("
 cap_hist, freq, comment
 1000, 2589, singe tag;  not seen again
 1010,  218, A single tag; recovered
 1111,   35, BB double tagged; recovered with double tags
 111X,   24, Double tagged; recovered with a single tag but cannot tell which tag was lost
 1100,  432, double tagged; not seen again
 0010, 11167, U - apparently captured for first time at t2 but this includes fish that lost all tags"
 )

data_kokanee_tagloss <- read.csv(data.csv, header=TRUE, as.is=TRUE, strip.white=TRUE, )
data_kokanee_tagloss$cap_hist <- formatC(data_kokanee_tagloss$cap_hist, width=4, digits=0, format="f", flag="0")

usethis::use_data(data_kokanee_tagloss, overwrite = TRUE)
