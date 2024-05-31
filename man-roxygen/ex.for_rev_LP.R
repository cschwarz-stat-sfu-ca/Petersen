# Example of combined forward and reverse MR

# get some data
n1 <-  500
n2 <- 1500
m2 <-  150

f.data <- data.frame(cap_hist=c("10","11","01"), freq=c(n1 - m2,  m2,  n2 - m2))
f.data

E = 1500
E.SE = 150
G = .2
G.SE = .05

res <- LP_for_rev_fit(data=f.data,
                      E=E,
                      E.SE=E.SE,
                      G=G,
                      G.SE=G.SE)
