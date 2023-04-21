## code to prepare data set for SPAS analysis

# This is the harrison river data. Fish were captured at a trap downstream
# over several weeks, and then 6 different spawning areas are surveyed.
# It is a combination of temporand and geographic stratitication.

data.csv <- textConnection("
  4   ,      2   ,      1   ,     1   ,     0   ,     0   ,   130
 12   ,      7   ,     14   ,     1   ,     3   ,     0   ,   330
  7   ,     11   ,     41   ,     9   ,     1   ,     1   ,   790
  1   ,     13   ,     40   ,    12   ,     9   ,     1   ,   667
  0   ,      1   ,      8   ,     8   ,     3   ,     0   ,   309
  0   ,      0   ,      0   ,     0   ,     0   ,     1   ,    65
744   ,   1187   ,   2136   ,   951   ,   608   ,   127   ,     0")

har.data <- as.matrix(read.csv(data.csv, header=FALSE))
har.data

# create capture histories for this data
# Notice that we use .. as a separator between the temporal strata

# Fish captured in the temporal stratum release but not recaptured

cap_hist <- paste0(formatC(1:(nrow(har.data)-1), width=2, digits=0, format="f", flag="0"), "..",
                         formatC(0, width=2, digits=0, format="f",flag="0"))
freq     <- as.vector(har.data[,ncol(har.data),drop=TRUE][-nrow(har.data)])
data1 <- data.frame(cap_hist=cap_hist, freq=freq)

# Fish released and recaptured


data2 <- plyr::adply(har.data[-nrow(har.data),-ncol(har.data)],1,function(x){
   cap_hist <- paste0("  ..",
                      letters[1:length(x)])
   freq     <- as.vector(x)
   data.frame(cap_hist=cap_hist, freq=freq)
})
substr(data2$cap_hist,1,2) <- rep( formatC(1:(nrow(har.data)-1), width=2, digits=0, format="f", flag="0"),
                                   each=ncol(har.data)-1)

# Fish newly captured in spawning areas
data3 <- data.frame(
    cap_hist = as.vector(paste0(formatC(0          ,width=2, digits=0, format="f",flag="0"), "..",
                         letters[1:(ncol(har.data)-1)])),
    freq     = as.vector(har.data[ncol(har.data),-ncol(har.data)])
)


data_spas_harrison <- plyr::rbind.fill(data1, data2, data3)
data_spas_harrison <- cbind(data_spas_harrison, Petersen::split_cap_hist(data_spas_harrison$cap_hist,
                                                               sep="..", prefix="..st"))

data_spas_harrison <- data_spas_harrison[ order(data_spas_harrison$..st2, data_spas_harrison$cap_hist), ]
data_spas_harrison <- data_spas_harrison[,c("cap_hist","freq")]
rownames(data_spas_harrison) <- NULL
data_spas_harrison

usethis::use_data(data_spas_harrison, overwrite = TRUE)
