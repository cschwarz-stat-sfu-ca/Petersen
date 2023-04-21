## code to prepare btspas diagonal example1 dataset goes here

data.csv <- textConnection(
'jweek, n1, m2,    u2, logflow
 9,    0,    0,  4135, 6.617212
10, 1465,   51, 10452, 6.51217
11, 1106,  121,  2199, 7.193686
12,  229,   25,   655, 6.960754
13,   20,    0,   308, 7.008376
14,  177,   17,   719, 6.761573
15,  702,   74,   973, 6.905753
16,  633,   94,   972, 7.062314
17, 1370,   62,  2386, 7.600188
18,  283,   10,   469, 8.246509
19,  647,   32,   897, 8.110298
20,  276,   11,   426, 8.035001
21,  277,   13,   407, 7.859965
22,  333,   15,   526, 7.774255
23, 3981,  242, 39969, 7.709116
24, 3988,   55, 17580, 7.653766
25, 2889,  115,  7928, 7.622105
26, 3119,  198,  6918, 7.593734
27, 2478,   80,  3578, 7.585063
28, 1292,   71,  1713, 7.291072
29, 2326,  153,  4212, 6.55556
30, 2528,  156,  5037, 6.227665
31, 2338,  275,  3315, 6.278789
32, 1012,  101,  1300, 6.273685
33,  729,   66,   989, 6.241111
34,  333,   44,   444, 6.687999
35,  269,   33,   339, 7.222566
36,   77,    7,   107, 7.097194
37,   62,    9,    79, 6.949993
38,   26,    3,    41, 6.168714
39,   20,    1,    23, 6.113682')

data <- read.csv(data.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)

# to reduce running times, we will collapse every 2 rows together
data$jweek <- trunc((data$jweek-1)/2)

data <- plyr::ddply(data, "jweek", plyr::summarize,
              n1 = sum(n1),
              m2 = sum(m2),
              u2 = sum(u2),
              logflow=mean(logflow))

# create capture histories for this data
# Notice that we use .. as a separator between the temporal strata

# Fish captured in same temporal stratum as released
data1 <- data
data1$cap_hist <- paste0(formatC(data1$jweek, width=2, digits=0, format="f", flag="0"), "..",
                         formatC(data1$jweek, width=2, digits=0, format="f",flag="0"))
data1$freq     <- data1$m2

# Fish released but not recaptured (temporal stratum 0)
data2 <- data
data2$cap_hist <- paste0(formatC(data1$jweek,width=2, digits=0, format="f",flag="0"), "..",
                         formatC(0          ,width=2, digits=0, format="f",flag="0"))
data2$freq     <- data2$n1 - data2$m2

# Fish newly captured in temporal stratum
data3 <- data
data3$cap_hist <- paste0(formatC(0          ,width=2, digits=0, format="f",flag="0"), "..",
                         formatC(data1$jweek,width=2, digits=0, format="f",flag="0"))
data3$freq     <- data3$u2


data_btspas_diag1 <- plyr::rbind.fill(data1, data2, data3)
data_btspas_diag1$n1 <- NULL
data_btspas_diag1$m2 <- NULL
data_btspas_diag1$u2 <- NULL

data_btspas_diag1 <- data_btspas_diag1[ order(data_btspas_diag1$jweek, data_btspas_diag1$cap_hist), ]
data_btspas_diag1 <- data_btspas_diag1[,c("cap_hist","freq","logflow")]
rownames(data_btspas_diag1) <- NULL
data_btspas_diag1

usethis::use_data(data_btspas_diag1, overwrite = TRUE)
