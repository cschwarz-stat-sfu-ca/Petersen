## code to prepare btspas nondiagonal example1 dataset goes here

# THis is based on the Conne River 2009 data, except truncated to a small number of releases


data.csv <- textConnection(
"Date,tagged,   0, 1, 2, 3,   Tot-recoveries,Untagged,Tot-caught,WaterTemp
2009-04-29, 25, 1, 2, 1, 0,     4,   0,   0,   7.8
2009-04-30, 75, 3, 8, 0, 0,    19, 133, 133,   7.0
2009-05-01, 97, 4,16, 0, 0,    20, 158, 161,   7.5
2009-05-02,127, 6, 0, 0, 3,    9,  128, 130,   7.7
2009-05-03,  0, 0, 0, 0, 0,     0,   0,   0,   9.3
2009-05-04,  0, 0, 0, 0, 0,     0,   0,   0,  11.7
2009-05-05,216, 8,12, 4, 1,    25,  859, 893,  11.6
2009-05-06,215, 1,26, 2, 0,    29,  427, 445,  11.7
2009-05-07,205,11,13,11, 0,    35,  849, 892,  11.1
2009-05-08,202, 3,23, 2, 0,    28,  488, 507,   9.8
2009-05-09  ,0, 0, 0, 0, 0,     0,   22,  22,  10.9
2009-05-10,  0, 0, 0, 0, 0,     0,   22,  22,  12.2
2009-05-11,  0, 0, 0, 0, 0,     0,    1,   1, 14.3")

data <- read.csv(data.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
data$Date <- as.Date(data$Date)

# create capture histories for this data
# Notice that we use .. as a separator between the temporal strata

# Fish captured in same temporal stratum as released
data$relday <- lubridate::yday(data$Date)

# create capture-histories for releases and recoveries

# Fish released but not recaptured (temporal stratum 0)
data1 <- plyr::adply(data, 1, function(x){
   cap_hist <- paste0(formatC(x$relday, digits=0, width=3, format="f", flag="0"),
                      "..",
                      formatC(c(x$relday+c(0,1,2,3),0), digits=0, width=3, format="f", flag="0"))
   freq= c(x$X0, x$X1, x$X2, x$X3, x$tagged - x$Tot.recoveries)
   data.frame(cap_hist=cap_hist,freq=freq)
   })

# Newly untagged fish
data2 <- plyr::adply(data, 1, function(x){
   cap_hist <- paste0(formatC(0,           digits=0, width=3, format="f", flag="0"),
                      "..",
                      formatC(c(x$relday), digits=0, width=3, format="f", flag="0"))
   freq= c(x$Untagged)
   data.frame(cap_hist=cap_hist,freq=freq)
   })

data_btspas_nondiag1 <- plyr::rbind.fill(data1, data2)
data_btspas_nondiag1 <- data_btspas_nondiag1[ order(data_btspas_nondiag1$relday, data_btspas_nondiag1$cap_hist),]
data_btspas_nondiag1 <- data_btspas_nondiag1[,c("cap_hist","freq")]
rownames(data_btspas_nondiag1)<- NULL

#remove 0 counts
data_btspas_nondiag1 <- data_btspas_nondiag1[data_btspas_nondiag1$freq > 0,]

usethis::use_data(data_btspas_nondiag1, overwrite = TRUE)



