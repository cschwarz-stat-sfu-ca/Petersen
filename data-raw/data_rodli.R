## code to prepare `rodli` dataset goes here

data_rodli <- data.frame(
  cap_hist=c("11", "10", "01"),
  freq    =c( 57,     52,      120)
)
usethis::use_data(data_rodli, overwrite = TRUE)
