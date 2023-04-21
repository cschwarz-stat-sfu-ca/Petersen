## code to prepare btspas nondiagonal example1 dataset goes here

# THis is based on the Nass River sockeye in 2014, with stratification by week


data.csv <- textConnection("
   xxx   ,  n1   , W23 , W24 , W25 , W26 , W27 , W28 , W29   ,   W30   ,   W31   ,   W32   ,   W33   ,   W34   ,   W35   ,   W36   ,   W37   ,   W38   ,   W39   ,   W40   ,   W41
   W23   ,  34   ,  0  ,  0  ,  0  ,  0  ,  0  ,  11 ,  2  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0
   W24   ,  280  ,  0  ,  0  ,  0  ,  0  ,  3  , 118 ,  18  ,  6  ,  3  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0
   W25   ,  607  ,  0  ,  0  ,  0  ,  0  ,  1  , 149 ,  76  ,  33  ,  16  ,  3  ,  1  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0
   W26   ,  603  ,  0  ,  0  ,  0  ,  0  ,  0  ,  42 ,  163  ,  69  ,  34  ,  1  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0
   W27   ,  388  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1  ,  39  ,  84  ,  35  ,  1  ,  3  ,  1  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0
   W28   ,  212  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  1  ,  29  ,  60  ,  14  ,  5  ,  2  ,  0  ,  1  ,  0  ,  0  ,  0  ,  0  ,  0
   W29   ,  468  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  3  ,  82  ,  84  ,  19  ,  11  ,  2  ,  1  ,  0  ,  0  ,  0  ,  0  ,  0
   W30   ,  586  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  20  ,  154  ,  56  ,  53  ,  9  ,  1  ,  1  ,  1  ,  0  ,  0  ,  0
   W31   ,  512  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  16  ,  44  ,  146  ,  83  ,  25  ,  3  ,  3  ,  0  ,  1  ,  0
   W32   ,  458  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  2  ,  52  ,  125  ,  82  ,  13  ,  14  ,  3  ,  0  ,  0
   W33   ,  479  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  83  ,  133  ,  43  ,  12  ,  5  ,  2  ,  0
   W34   ,  329  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  12  ,  96  ,  60  ,  29  ,  8  ,  0  ,  1
   W35   ,  248  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  19  ,  72  ,  37  ,  17  ,  7  ,  0
   W36   ,  201  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  7  ,  55  ,  23  ,  6  ,  1
   W37   ,  104  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  7  ,  3  ,  3  ,  0
   W00   ,  0   ,  0  ,  0  ,  0  ,  0  ,1420 ,27886  ,  9284  ,  4357  ,  11871  ,  14896  ,  9910  ,  16526  ,  17443  ,  16485  ,  6776  ,  4644  ,  2190  ,  1066  ,  166
")

data <- read.csv(data.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)

data$relday <- as.numeric(substring(data$xxx,2))

# create capture histories for this data
# Notice that we use .. as a separator between the temporal strata

# create capture-histories for releases and recoveries

# Fish released but not recaptured (temporal stratum 0)
data1 <- plyr::adply(data[-nrow(data),], 1, function(x){
  cap_hist <- paste0(formatC(x$relday, digits=0, width=3, format="f", flag="0"),
                     "..000")
  freq= x$n1 - sum(x[-c(1:2,length(x))])
  #browser()
  data.frame(cap_hist=cap_hist,freq=freq)
})


# recaptures
data2 <- plyr::adply(data[-nrow(data),], 1, function(x){
   recsel <- grepl("W", colnames(x))
   recweek = as.numeric(substring(colnames(x)[recsel],2))
   cap_hist <- paste0(formatC(x$relday, digits=0, width=3, format="f", flag="0"),
                      "..",
                      formatC(recweek, digits=0, width=3, format="f", flag="0"))
   freq= as.vector(unlist(x[recsel]))
   #browser()
   data.frame(cap_hist=cap_hist,freq=freq)
   })


# Newly untagged fish
data3 <- plyr::adply(data[nrow(data),,drop=FALSE], 1, function(x){
  recsel <- grepl("W", colnames(x))
  recweek = as.numeric(substring(colnames(x)[recsel],2))
  cap_hist <- paste0(formatC(0,           digits=0, width=3, format="f", flag="0"),
                      "..",
                      formatC(recweek, digits=0, width=3, format="f", flag="0"))
  freq= as.vector(unlist(x[recsel]))
  data.frame(cap_hist=cap_hist,freq=freq)
   })

data_btspas_nondiag1 <- plyr::rbind.fill(data1, data2, data3)
data_btspas_nondiag1 <- cbind(data_btspas_nondiag1, split_cap_hist(data_btspas_nondiag1$cap_hist, sep="..", prefix="..ts", make.numeric=TRUE))
data_btspas_nondiag1 <- data_btspas_nondiag1[ order(data_btspas_nondiag1$..ts2, data_btspas_nondiag1$..ts1),]

data_btspas_nondiag1 <- data_btspas_nondiag1[,c("cap_hist","freq")]
rownames(data_btspas_nondiag1)<- NULL

#remove 0 counts excepts for the newly capture animals
select <- data_btspas_nondiag1$freq > 0 |
          data_btspas_nondiag1$freq ==0 & substr(data_btspas_nondiag1$cap_hist,1,3) == '000'
data_btspas_nondiag1 <- data_btspas_nondiag1[select,]
rownames(data_btspas_nondiag1) <- NULL

temp <- cbind(data_btspas_nondiag1, split_cap_hist(data_btspas_nondiag1$cap_hist, sep="..", make.numeric=TRUE,
                                                   prefix="..ts"))


xtabs(freq ~ ..ts1+..ts2, data=temp, exclude=NULL, na.action=na.pass)
usethis::use_data(data_btspas_nondiag1, overwrite = TRUE)



