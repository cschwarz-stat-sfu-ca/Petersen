############################################################################
############################################################################
##
## The R-code given below contains all the R-code that was used for the 
## examples/results/simulations presented in the submitted JSS manuscript: 
## "The VGAM Package for Capture-Recapture Data Using the Conditional 
## Likelihood" by Yee, Stoklosa, and Huggins. We give the section numbers 
## where the R-code was used.
##
############################################################################
############################################################################

## [1.] Check that the latest version of "VGAM" has been installed.
onthismachine <- installed.packages()
dim.otm <- dimnames(onthismachine)
if (any(dim.otm[[1]] == "VGAM")) {
## VGAM installed, but need to check it is recent enough.
  if (packageVersion("VGAM") < "0.9.4") {
    warning("'VGAM' 0.9-4 or later required. Updating it now.")
    install.packages("VGAM")
  }
} else {
## VGAM is not installed so install it.
  install.packages("VGAM")
}

## [2.] Check that the latest version of "mra" has been installed.
onthismachine <- installed.packages()
dim.otm <- dimnames(onthismachine)
if (any(dim.otm[[1]] == "mra")) {
## Do nothing since mra has been already installed.
} else {
## mra  is not installed so install it.
  install.packages("mra")
}

## [3.] Check that the latest version of "RMark" has been installed.
## Note: installation details are at
##       http://www.phidot.org/software/mark/downloads/
## Also note that RMark requires external software which is not 
## currently available for Linux.
onthismachine <- installed.packages()
dim.otm <- dimnames(onthismachine)
if (any(dim.otm[[1]] == "RMark")) {
## Do nothing since RMark has been already installed.
} else {
## RMark is not installed so install it.
  install.packages("RMark")
}

## [4.] Check that the latest version of "MASS" has been installed.
onthismachine <- installed.packages()
dim.otm <- dimnames(onthismachine)
if (any(dim.otm[[1]] == "MASS")) {
## Do nothing since MASS has been already installed.
} else {
## MASS is not installed so install it.
  install.packages("MASS")
}

## [5.] Check that the latest version of "splines" has been installed.
onthismachine <- installed.packages()
dim.otm <- dimnames(onthismachine)
if (any(dim.otm[[1]] == "splines")) {
## Do nothing since splines has been already installed.
} else {
## splines is not installed so install it.
  install.packages("splines")
}

## Now load in the packages and run the code!
library("VGAM")
library("mra")
library("RMark")
library("MASS")
library("splines")
rm(list = ls())

options(prompt = "R> ", continue = "+  ")

############################################################################

## Section 5.

args(posbinomial)

args(posbernoulli.t)

args(posbernoulli.b)

args(posbernoulli.tb)

############################################################################

## We use two example data sets (replication material), both of which are 
## contained in the VGAM package. Code used for the simulations are
## also given.

## Section 6.1. (Example using the deer mice data).

head(deermice, 4)

deermice$age <- 2 - as.numeric(deermice$age)  
deermice$sex <- 1 - as.numeric(deermice$sex)  

M.0 <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1,
            posbernoulli.t(parallel = TRUE ~ 1), data = deermice)

M.b <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1, 
            posbernoulli.b, data = deermice)

M.t <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1, 
            posbernoulli.t, data = deermice)

M.h <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ weight + sex + age,
            posbernoulli.t(parallel = TRUE ~ weight + sex + age),
            data = deermice)

M.th <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ weight + sex + age,
             posbernoulli.t, data = deermice)

M.tb <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1, 
             posbernoulli.tb, data = deermice)

M.bh <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ weight + sex + age,
             posbernoulli.b, data = deermice)

M.tbh <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ weight + sex + age,
             posbernoulli.tb, data = deermice)

c(M.bh@extra$N.hat, M.bh@extra$SE.N.hat)
c(logLik(M.bh), AIC(M.bh))

Table <-
  rbind(c(round(M.tbh@extra$N.hat, 2), round(M.bh@extra$N.hat, 2), 
          round(M.tb@extra$N.hat, 2),  round(M.th@extra$N.hat, 2),
          round(M.h@extra$N.hat, 2),   round(M.b@extra$N.hat, 2),
          round(M.t@extra$N.hat, 2),   round(M.0@extra$N.hat, 2)), 
               
c(round(M.tbh@extra$SE.N.hat, 2), round(M.bh@extra$SE.N.hat, 2), 
  round(M.tb@extra$SE.N.hat, 2),  round(M.th@extra$SE.N.hat, 2),
  round(M.h@extra$SE.N.hat, 2),   round(M.b@extra$SE.N.hat, 2),
  round(M.t@extra$SE.N.hat, 2),   round(M.0@extra$SE.N.hat, 2)), 

-2 * c(round(logLik(M.tbh), 2), round(logLik(M.bh), 2),
       round(logLik(M.tb), 2),  round(logLik(M.th), 2),
       round(logLik(M.h), 2),   round(logLik(M.b), 2),
       round(logLik(M.t), 2),   round(logLik(M.0), 2)), 

c(round(AIC(M.tbh), 2), round(AIC(M.bh), 2), round(AIC(M.tb), 2), 
  round(AIC(M.th), 2), round(AIC(M.h), 2), round(AIC(M.b), 2), 
  round(AIC(M.t), 2), round(AIC(M.0), 2)))

colnames(Table) <- c("M.tbh", "M.bh", "M.tb", "M.th", "M.h",
                     "M.b", "M.t", "M.0")
rownames(Table) <- c("N.hat", "S.E.","-2ln(L)", "AIC")

Table

round(coef(M.bh), 2)
round(sqrt(diag(vcov(M.bh))), 2)

fit.bh <- vgam(cbind(y1, y2, y3, y4, y5, y6) ~
               s(weight, df = 3) + sex + age,
               posbernoulli.b, data = deermice)
plot(fit.bh, se = TRUE, las = 1, lcol = "blue", scol = "orange",
     rcol = "purple", scale = 5)

summary(fit.bh)  # The p-value tests for linearity

par(mfrow = c(2, 2))
par(las = 1, cex = 1.1, mar = c(3.8, 4, 0.5, 0.2) + 0.1)
par(mgp = c(2.3, 1, 0)) 

plot(fit.bh, se = TRUE, las = 1, lcol = "blue", scol = "orange",
     rcol = "purple", scale = 5, mgp = c(2.0, 1, 0))

############################################################################

## Section 6.2. (Example using the prinia bird data).

data("prinia", package = "VGAM")

head(prinia, 4)[, 1:4]

M.h.GAM <- vgam(cbind(cap, noncap) ~ s(length, df = 3) + fat, 
                posbinomial(omit.constant = TRUE,
                            parallel = TRUE ~ s(length, df = 3) + fat),
                data = prinia)

M.h.GAM@extra$N.hat     
M.h.GAM@extra$SE.N.hat

par(mfrow = c(1, 1))

plot.info <- plot(M.h.GAM, se = TRUE, las = 1, lcol = "blue",
                  scol = "orange", rcol = "purple", scale = 5)

tau <- 19
D <- nrow(prinia)
info.fit2 <- plot.info@preplot[[1]]
fat.effect <- coef(M.h.GAM)["fat"]  # Fat index effect.
intercept <- coef(M.h.GAM)["(Intercept)"]  # Overall intercept.

ooo <- order(info.fit2$x)
centering.const <- mean(prinia$length) - coef(M.h.GAM)["s(length, df = 3)"]

plotframe <- data.frame(lin.pred.b = intercept + fat.effect * 1 +
                                     centering.const + info.fit2$y[ooo],
                        lin.pred.0 = intercept + fat.effect * 0 +
                                     centering.const + info.fit2$y[ooo],
                        x2 = info.fit2$x[ooo])

plotframe <- transform(plotframe,
                       up.lin.pred.b = lin.pred.b + 2*info.fit2$se.y[ooo],
                       lo.lin.pred.b = lin.pred.b - 2*info.fit2$se.y[ooo],
                       up.lin.pred.0 = lin.pred.0 + 2*info.fit2$se.y[ooo],
                       lo.lin.pred.0 = lin.pred.0 - 2*info.fit2$se.y[ooo])

plotframe <- transform(plotframe,
                       fv.b    = logit(lin.pred.b,    inverse = TRUE),
                       up.fv.b = logit(up.lin.pred.b, inverse = TRUE),
                       lo.fv.b = logit(lo.lin.pred.b, inverse = TRUE),
                       fv.0    = logit(lin.pred.0,    inverse = TRUE),
                       up.fv.0 = logit(up.lin.pred.0, inverse = TRUE),
                       lo.fv.0 = logit(lo.lin.pred.0, inverse = TRUE))

with(plotframe,
     matplot(x2, cbind(up.fv.b, fv.b, lo.fv.b), type = "l", col = "blue",
             lty = c(2, 1, 2), las = 1, cex.lab = 1.5, lwd = 2,
             main = "", ylab = "", xlab = "Wing length (standardized)"))
mtext( ~ hat(p), side = 2, cex = 1.4, line = 4, adj = 0.5, las = 1)
with(plotframe, matlines(x2, cbind(up.fv.0, fv.0, lo.fv.0),
                         col = "darkorange", lty = c(2, 1, 2)), lwd = 2)
legend("topleft", legend = c("Fat present", "Fat not present"), bty = "n",
       lwd = 2, col = c("blue", "darkorange"), merge = TRUE, cex = 1.5)



############################################################################

## Section 6.3. (A time-varying covariate example).
## Reference: Huggins, R. H., Biometrika, 1989, 76(1), 133--140.

head(Huggins89table1, 4)

Hdata <- transform(Huggins89table1, x3.tij = t01,
                   T02 = t02, T03 = t03, T04 = t04, T05 = t05, T06 = t06,
                   T07 = t07, T08 = t08, T09 = t09, T10 = t10)
Hdata <- subset(Hdata,
                y01 + y02 + y03 + y04 + y05 + y06 + y07 + y08 + y09 + y10 > 0)


fit.th <-
   vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08, y09, y10) ~  x2 + x3.tij,
        xij = list(x3.tij ~ t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 +
                            t09 + t10 - 1),
        posbernoulli.t(parallel.t = TRUE ~ x2 + x3.tij), 
        data = Hdata, trace = FALSE, 
        form2 = ~ x2 + x3.tij + t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 +
                                t09 + t10)
constraints(fit.th, matrix = TRUE)



fit.tbh <-
  vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08, y09, y10) ~  x2 + x3.tij,
       xij = list(x3.tij ~ t01 + t02 + t03 + t04 + t05 + t06 +
                                 t07 + t08 + t09 + t10 +
                                 T02 + T03 + T04 + T05 + T06 +
                                 T07 + T08 + T09 + T10 - 1),
       posbernoulli.tb(parallel.t = TRUE ~ x2 + x3.tij),
       data = Hdata, trace = FALSE,
       form2 = ~  x2 + x3.tij +
                  t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 + t09 + t10 +
                        T02 + T03 + T04 + T05 + T06 + T07 + T08 + T09 + T10)

c(logLik(fit.th), AIC(fit.th))
c(logLik(fit.tbh), AIC(fit.tbh))


head(constraints(fit.tbh, matrix = TRUE), 4)
tail(constraints(fit.tbh, matrix = TRUE), 4)


coef(fit.tbh)
sqrt(diag(vcov(fit.tbh))) 


fit.tbh@extra$N.hat
fit.tbh@extra$SE.N.hat





## Alternative method to fit the above model.
## Uses Select() to select certain variables from pdata.

Hdata <- subset(Huggins89table1, rowSums(Select(Huggins89table1, "y")) > 0)
Hdata.T <- Select(Hdata, "t")  # A 10-column submatrix copy
colnames(Hdata.T) <- gsub("t", "T", colnames(Hdata.T))  # Rename colnames
Hdata <- data.frame(Hdata, Hdata.T)  # A copy made of the "t" variables
Hdata <- transform(Hdata, x3.tij = y01)
Form2 <- Select(Hdata, prefix = TRUE, as.formula = TRUE)  # Bloated
Xij   <- Select(Hdata, c("t", "T"), as.formula = TRUE,
                sort = FALSE, rhs = "0", lhs = "x3.tij", exclude = "T01")
fit.tbh <- vglm(Select(Hdata, "y") ~ x2 + x3.tij,
                form2 = Form2,  xij = list(Xij),
                posbernoulli.tb(parallel.t = TRUE ~ x2 + x3.tij),
                data = Hdata, trace = FALSE)
coef(fit.tbh)






############################################################################

## Section 6.4. (Ephemeral and enduring memory effects)
  
## Reference: Yang H. C. and Chao A., Biometrics, 2005, 61(4), 1010--1017.


## This is a lag-1 model (ephemeral memory lasts one time period):

deermice <- transform(deermice, Lag1 = y1)
M.tbh.lag1 <-
  vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight + Lag1,
       posbernoulli.tb(parallel.t = FALSE ~ 0,
                       parallel.b = FALSE ~ 0,
                       drop.b = FALSE ~ 1),
       xij = list(Lag1 ~ fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                         fill(y5) + fill(y6) +
                         y1 + y2 + y3 + y4 + y5),
       form2 = ~ sex + weight + Lag1 +
                 fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                 fill(y5) + fill(y6) +
                 y1 + y2 + y3 + y4 + y5 + y6,
       data = deermice)

coef(M.tbh.lag1)
coef(M.tbh.lag1)["(Intercept):1"]  # The enduring effect
coef(M.tbh.lag1)["Lag1"]           # The ephemeral effect




## This is an alternative method to fit the same model as M.tbh.lag1:

deermice <- transform(deermice, Lag1 = y1)
deermice <- transform(deermice, f1 = y1, f2 = y1, f3 = y1, f4 = y1,
                                f5 = y1, f6 = y1)
tau <- 6
H2 <- H3 <- cbind(rep(1, 2*tau-1))
H4 <- cbind(c(rep(0, tau), rep(1, tau-1)))
M.tbh.lag1.method2 <-
  vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight + Lag1,
       posbernoulli.tb(parallel.b = TRUE ~ 0, parallel.t = TRUE ~ 0),
       constraints = list("(Intercept)" = cbind(H4, 1), sex = H2, weight= H3, 
                          Lag1 = H4),
       xij = list(Lag1 ~ f1 + f2 + f3 + f4 + f5 + f6 +
                         y1 + y2 + y3 + y4 + y5),
       form2 = Select(deermice, prefix = TRUE, as.formula = TRUE),
       data = deermice)
coef(M.tbh.lag1.method2)



## This is an lag-2 model: ephemeral memory lasts for two time periods.

deermice <- transform(deermice, Lag2 = y1)
M.bh.lag2 <-
  vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight + Lag2,
       posbernoulli.tb(parallel.t = FALSE ~ 0,
                       parallel.b = FALSE ~ 0,
                       drop.b = FALSE ~ 1),
       xij = list(Lag2 ~ fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                         fill(y5) + fill(y6) +
                         y1 + pmax(y1, y2) + pmax(y2, y3) + pmax(y3, y4) + 
                         pmax(y4, y5)),
       form2 = ~ sex + weight + Lag2 +
                 fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                 fill(y5) + fill(y6) +
                 y1 + pmax(y1, y2) + pmax(y2, y3) + pmax(y3, y4) + 
                 pmax(y4, y5) + y6,
       data = deermice)
coef(M.bh.lag2)



  

############################################################################

## Section 6.5. (Some timing and reliability tests)
  
## Time comparisons for VGAM and mra R-packages. 

## Using the VGAM R-Package.
  
start_0 <- Sys.time()
M_0 <- vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08,
                  y09, y10, y11, y12, y13, y14,
                  y15, y16, y17, y18, y19) ~ 1,
            posbernoulli.t(parallel = TRUE ~ 1), data = prinia) 
end_0 <- Sys.time()

start_t <- Sys.time()
M_t <- vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08,
                  y09, y10, y11, y12, y13, y14, y15,
                  y16, y17, y18, y19) ~ 1,
            posbernoulli.t, data = prinia) 
end_t <- Sys.time()

start_h <- Sys.time()
M_h <- vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08,
                  y09, y10, y11, y12, y13, y14, y15,
                  y16, y17, y18, y19) ~ fat,
            posbernoulli.t(parallel = TRUE ~ fat), data = prinia)
end_h <- Sys.time()

start_b <- Sys.time()
M_b <- vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08,
                  y09, y10, y11, y12, y13, y14, y15,
                  y16, y17, y18, y19) ~ 1, 
            posbernoulli.b, data = prinia)
end_b <- Sys.time()

start_tb <- Sys.time()
M_tb <- vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08,
                   y09, y10, y11, y12, y13, y14, y15,
                   y16, y17, y18, y19) ~ 1, 
             posbernoulli.tb, data = prinia)
end_tb <- Sys.time()

start_th <- Sys.time()
M_th <- vglm(cbind(y01, y02, y03, y04, y05, y06, y07,
                   y08, y09, y10, y11, y12, y13,
                   y14, y15, y16, y17, y18, y19) ~ fat, 
             posbernoulli.t, data = prinia)
end_th <- Sys.time()

start_bh <- Sys.time()
M_bh <- vglm(cbind(y01, y02, y03, y04, y05, y06, y07,
                    y08, y09, y10, y11, y12, y13,
                    y14, y15, y16, y17, y18, y19) ~ fat, 
              posbernoulli.b, data = prinia)
end_bh <- Sys.time()

start_tbh <- Sys.time()
M_tbh <- vglm(cbind(y01, y02, y03, y04, y05, y06, y07,
                   y08, y09, y10, y11, y12, y13,
                   y14, y15, y16, y17, y18, y19) ~ fat, 
             posbernoulli.tb, data = prinia)
end_tbh <- Sys.time()

VGAM_times1 <- c(end_0 - start_0, end_t - start_t, end_h - start_h, 
                 end_b - start_b, end_tb - start_tb, end_th - start_th,
                 end_bh - start_bh, end_tbh - start_tbh)
VGAM_times <- c(VGAM_times1, sum(VGAM_times1))

## Using the mra R-package.

cap.hist <- matrix(c(as.numeric(unlist(prinia[, 5:23]))), ncol = 19)

ct <- as.factor(1:ncol(cap.hist))
attr(ct, "nan") <- nrow(cap.hist)
x.obs2 <- as.numeric(prinia$fat)
attr(x.obs2, "ns") <- ncol(cap.hist)

start_0 <- Sys.time()
hug.0 <- F.huggins.estim( ~ 1, NULL, cap.hist)
end_0 <- Sys.time()
start_t <- Sys.time()
hug.t <- F.huggins.estim( ~ tvar(ct), NULL, cap.hist)
end_t <- Sys.time()
start_h <- Sys.time()
hug.h <- F.huggins.estim( ~ ivar(x.obs2), NULL, cap.hist)
end_h <- Sys.time()


# Note: Models _b and _tb omitted here because they may crash,
# e.g., on Linux machines.


start_th <- Sys.time() 
hug.th <- F.huggins.estim( ~ ivar(x.obs2) + tvar(ct), NULL, cap.hist) 
end_th <- Sys.time()

mra_times1 <- c(end_0 - start_0, end_t - start_t, end_h - start_h, 
                end_b - start_b, end_tb - start_tb, end_th - start_th,
                NA, NA)
mra_times <- c(mra_times1, sum(mra_times1, na.rm = TRUE))

## Using the RMark R-package.

prinia$ch <- apply(prinia[, -(1:4)], 1, paste, collapse = "")
prinia.proc <- process.data(prinia, model = "Huggins")
prinia.ddl <- make.design.data(prinia.proc)

start_0 <- Sys.time()
M0 <- try(mark(prinia.proc, prinia.ddl, 
               model.parameters = list(p = list(formula = ~ 1,
                                           share = TRUE))))
end_0 <- Sys.time()

start_t <- Sys.time()
Mt <- try(mark(prinia.proc, prinia.ddl, 
               model.parameters = list(p = list(formula = ~ time,
                                           share = TRUE))))
end_t <- Sys.time()

start_h <- Sys.time()
Mh <- try(mark(prinia.proc, prinia.ddl, 
               model.parameters = list(p = list(formula = ~ fat,
                                           share = TRUE))))
end_h <- Sys.time()

start_b <- Sys.time()
Mb <- try(mark(prinia.proc, prinia.ddl, 
               model.parameters = list(p = list(formula = ~ 1 + c,
                                           share = TRUE))))
end_b <- Sys.time()

start_tb <- Sys.time()
Mtb <- try(mark(prinia.proc, prinia.ddl, 
                model.parameters = list(p = list(formula = ~ time + c,
                                            share = TRUE))))
end_tb <- Sys.time()

start_th <- Sys.time()
Mth <- try(mark(prinia.proc, prinia.ddl, 
                model.parameters = list(p = list(formula = ~ time + fat,
                                            share = TRUE))))
end_th <- Sys.time()

start_bh <- Sys.time()
Mbh <- try(mark(prinia.proc, prinia.ddl, 
                model.parameters = list(p = list(formula = ~ fat + c,
                                            share = TRUE))))
end_bh <- Sys.time()

start_tbh <- Sys.time()
Mtbh <- try(mark(prinia.proc, prinia.ddl, 
                 model.parameters = list(p = list(formula = ~ time + fat + c,
                                             share = TRUE))))
end_tbh <- Sys.time()

RMark_times1 <- c(end_0 - start_0, end_t - start_t, end_h - start_h, 
                 end_b - start_b, end_tb - start_tb, end_th - start_th,
                 end_bh - start_bh, end_tbh - start_tbh)
RMark_times <- c(RMark_times1, sum(RMark_times1))

## Combine all the results in a table.

table_time <- cbind(signif(VGAM_times,  digits = 2),
                    signif(mra_times,   digits = 2),
                    signif(RMark_times, digits = 4))
rownames(table_time) <- c('M_0', 'M_t', 'M_h', 'M_b',
                          'M_tb', 'M_th', 'M_bh', 'M_tbh', 'Total')
colnames(table_time) <- c('VGAM comp. times', 'mra comp. times', 
                          'RMark comp. times')
table_time

## Simulations. 
## The simulations below are for GAM M_th models and follow 
## a similar simulation set-up as in Stoklosa and Huggins (2012).
## Reference: Stoklosa, J and Huggins, R. M., Computational 
## Statistics & Data Analysis, 2012, 56(2), 408--417.

## Simulations are compared with the GAM methods of Stoklosa and 
## Huggins (2012). We therefore require their code which is 
## sourced below from a separate R-file called GAM_th_sims_prog.R. 
## Also, see Stoklosa and Huggins (2012) for further details on 
## the fitting code.

source("GAM_th_sims_prog.R")
N.sim <- 100  # True population size.
tau <- 8      # Number of capture occasions.
sim <- 500    # Simulations will take about 6 hours, reduce sim if needed.
## Linear predictor with unknown smooths with inverse logit function.

P.s1 <- function(x1, x2, ti, z) {
  1/(1 + exp( - ( -1.2 + 2*sin(2*x1 - 1) + cos(2*x2) - 0.3*z + ti)))
}

## Simulated capture history matrix and covariates.

sim.data1 <- function(N.sim, tau, i) {      
  set.seed(666+i);  # Set seed to reproduce the results here.
  x.s1 <- runif(N.sim, -3, 0)
  x.s2 <- runif(N.sim, -3, 0)
  x1 <- rep(x.s1, each = tau)
  x2 <- rep(x.s2, each = tau)
  x1 <- matrix(x1, N.sim, tau, byrow = TRUE)
  x2 <- matrix(x2, N.sim, tau, byrow = TRUE)
  z.s <- rbinom(N.sim, 1, 0.6)
  z <- rep(z.s, each = tau)
  z <- matrix(z, N.sim, tau, byrow = TRUE)
  ti <- matrix(rep(c(0, log(2), 0, log(2), 0, log(2), 0, log(2)), N.sim),
               N.sim, tau, byrow = TRUE)     
  P.sim <- P.s1(x1, x2, ti, z)
  YRU <- matrix(runif(N.sim*tau), N.sim, tau)
  Yij <- (YRU < P.sim) + 0    
  Y.s <- apply(Yij, 1, sum)
  nocap <- (Y.s == 0)  
  Yij <- Yij[!nocap, ]
  x.s1 <- x.s1[!nocap]
  x.s2 <- x.s2[!nocap]
  z.s <- z.s[!nocap]
  n <- length(z.s)
  list(Yij = Yij, x1 = x.s1, x2 = x.s2, z = z.s, n = n)
}

## 95% coverage probability function.

cov.pc <- function(true, est, esd, per = 0.95) {
  aa <- (abs(est - true)< -qnorm((1 - per)/2)*esd)
  apply(aa, 2, mean)
}

start <- Sys.time()
B <- NULL
for (i in 1:sim) {
  cap.sim <- sim.data1(N.sim, tau, i)
  yij <- c(t(cap.sim$Yij))
  x.obs1 <- cap.sim$x1
  x.obs2 <- cap.sim$x2
  z.obs <- cap.sim$z
  hist <- cap.sim$Yij
  n <- nrow(hist)
  
  data.vgam <- data.frame(cbind(hist, x.obs1, x.obs2, z.obs))
  colnames(data.vgam) <- c("y1", "y2", "y3", "y4", "y5",
                           "y6", "y7", "y8", "x1", "x2", "z1")
  
## Linear Model.
  
  ests1_th <- vglm(cbind(y1, y2, y3, y4, y5, y6, y7, y8) ~  x1 + x2 + z1,
                   posbernoulli.t, data = data.vgam)
  
  VGLM.N.est_th <- ests1_th@extra$N.hat
  VGLM.varN.est_th <- (ests1_th@extra$SE.N.hat)^2

## GAM fit using methods/code from Stoklosa and Huggins (2012).

  cz <- cbind(rep(1, n), x.obs1, x.obs2, z.obs)
  czall <- covaz.all(cz, hist, n)
  czall <- czall[, -ncol(czall)]
  Y <- yij
  X <- czall
  X.d1 <- as.matrix(X[, 2])
  X.d2 <- as.matrix(X[, 3])
  X.b <- as.matrix(X[, c(1, 4:ncol(X))])
  lambda1 <- 0.001       
  lambda2 = lambda2.int <- 0.001  # Starting lambda values
  n.dx <- 8  # Number of knots.
  p.ord <- 2 # Polynomial order.
  nGCV <- 5  # Number of GCV steps.
  GCV.list <- NULL
  for(j in 1:nGCV)  {
    est1 <- GAM_cond_th(X.d1, X.d2, X.b, Y, tau, n, n.dx, lambda1, lambda2, nGCV, p.ord)
    GCV.list <- cbind(GCV.list, est1$GCV_th)
    lambda2 <- lambda2*10
    }
  lambda1_opt <- lambda1*10^(which(GCV.list == min(GCV.list), arr.ind = TRUE)[1]-1)
  lambda2_opt <- lambda2.int*10^(which(GCV.list == min(GCV.list), arr.ind = TRUE)[2]-1)
  ests2_th <- GAM_th(X.d1, X.d2, X.b, Y, tau, n, n.dx, lambda1_opt, lambda2_opt, p.ord)
  GAM.N.est_th <- ests2_th$N.est_th
  GAM.varN.est_th <- ests2_th$varN.est_th 

## VGAM.
  
  ests3_th <- vgam(cbind(y1, y2, y3, y4, y5, y6, y7, y8) ~
                   s(x1, df = 3) + s(x2, df = 3) + z1,
                   posbernoulli.t, data = data.vgam)
  
  VGAM.N.est_th <- ests3_th@extra$N.hat
  VGAM.varN.est_th <- (ests3_th@extra$SE.N.hat)^2
  
  B <- rbind(B, c(cap.sim$n,
                  VGLM.N.est_th, GAM.N.est_th, VGAM.N.est_th,
                  sqrt(VGLM.varN.est_th),
                  sqrt(GAM.varN.est_th),
                  sqrt(VGAM.varN.est_th)
                  ))
  }

means <- apply(B, 2, mean)
sd.means <- apply(B, 2, sd)
medians <- apply(B, 2, median)
mads <- apply(B, 2, mad)
MSE <- apply((B[, c(2:4)] - N.sim)^2, 2, mean)
MAE <- apply(abs(B[, c(2:4)] - N.sim), 2, median)
coverage <- cov.pc(N.sim, B[, c(2:4)], B[, c(5:7)])

results<-rbind(means[1:4],medians[1:4],sd.means[1:4],mads[1:4],
               c(NA,means[5:7]),c(NA,medians[5:7]),
               c(NA,MSE),c(NA,MAE),
               c(NA,coverage))
results

end <- Sys.time()
end - start

