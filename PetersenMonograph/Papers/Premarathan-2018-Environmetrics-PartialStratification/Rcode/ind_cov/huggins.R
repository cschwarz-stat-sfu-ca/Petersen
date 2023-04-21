
###############################################################################
############ Capture Closed models  and Huggins models using RMark ############
###############################################################################

setwd("U:\\Lasantha/Research/Rcode")
library(RMark)
options(width=350)

## "convert.inp" - Convert MARK input file to RMark dataframe
##    covariates - names to be assigned to the covariates defined in the inp file
##    first two colums in the .inp file are names as "ch" and "freq"
data <- convert.inp("test", covariates=c("sex", "length"))
head(data)
data$length.sq= (data$length)^2 # add a column of (length)^2 to the data frame
data$sex  <- as.factor(ifelse(data$sex==1, "Male", "Female"))
data
summary(data)

# Process encounter history dataframe for MARK analysis
#Prior to analyzing the data, this function initializes several variables 
# (e.g., number of capture occasions, time intervals) that are often specific
# to the capture-recapture model being fitted to the data. It also is used to
# 1) define groups in the data that represent different levels of one or more factor
#    covariates (e.g., sex), 
# 2) define time intervals between capture occasions (if not 1), and 
# 3) create an age structure for the data, if any.

# Process data
data.processed=process.data(data,model="Huggins",groups=c("sex"))
#
# Create default design data
#data.ddl=make.design.data(data.processed)

#  Define parameter models
p.dot = list(formula=~1,share=TRUE)
p.time = list(formula=~time,share=TRUE)
p.sex = list(formula=~sex,share=TRUE)
p.time.sex =  list(formula=~time+sex,share=TRUE)

p.length.sq =list(formula=~length+length.sq,share=TRUE)
p.length.time =list(formula=~length+length.sq + time,share=TRUE)
p.length.sex =list(formula=~length+length.sq + sex,share=TRUE)
p.length.time.sex =list(formula=~length+length.sq + sex+time,share=TRUE)



############  Capture Closed models ###########################################
# Not vary by sex or time , p=c 
closed.tdot=mark(data.processed,model="Closed",model.parameters=list(p=p.dot))

#  constant p and constant c but different, Not vary by sex or time 
closed.tdot.2=mark(data.processed,model="Closed")

# Not vary by sex but time varying,   p=c ,
closed.time=mark(data.processed,model="Closed",model.parameters=list(p=p.time))

# Closed heterogeneity : (vary by sex but Not time), p=c
closed.sex=mark(data.processed,model="Closed",
                model.parameters=list(p=p.sex),adjust=TRUE)

# Closed heterogeneity : (vary by time and sex), p=c
closed.sex=mark(data.processed,model="Closed",
                model.parameters=list(p=p.time.sex),adjust=TRUE)



######################### Huggins models #####################################

# Not vary by sex or time , p=c 
huggins.tdot.sdot=mark(data.processed,model="Huggins", 
                       model.parameters=list(p=p.length.sq))

# Not vary by sex but time varying,   p=c ,
huggins.time=mark(data.processed,model="Huggins",
                  model.parameters=list(p=p.length.time))

# Huggins heterogeneity (vary by sex but not time), p=c
huggins.sex=mark(data.processed,model="Huggins",
                 model.parameters=list(p=p.length.sex),adjust=TRUE)

# Huggins heterogeneity (vary by time and sex), p=c
huggins.time.sex=mark(data.processed,model="Huggins",
                      model.parameters=list(p=p.length.time.sex),adjust=TRUE)


##################################################################################

huggins.sex$results$beta

# create values of length to use for predictions
min.length=min(data$length)
max.length=max(data$length)
length.values=seq(from=min.length,to=max.length,length=25)

# Create design dataframes for MARK model specification
data.ddl=make.design.data(data.processed)

#Compute estimates of real parameters for multiple covariate values
# Computes real estimates for a dataframe of covariate values and the var-cov matrix of the real estimates.

# make predictions for rows of 'data.ddl' associated with
# females (par.index=1) and males (par.index=3)
pred.model <- covariate.predictions(huggins.time.sex, 
                                  data=data.frame(length=length.values,length.sq=length.values^2),
                                  indices=c(1,3))

# store values of sex in pred.model
female.rows=which(pred.model$estimates$par.index==1)
male.rows=which(pred.model$estimates$par.index==3)
pred.model$estimates$sex <- NA
pred.model$estimates$sex[female.rows] <- "female"
pred.model$estimates$sex[male.rows] <- "male"
head(pred.model$estimates)


pred=pred.model$estimates

library(ggplot2)
pred.ggplot <- ggplot(pred, aes(x=length, y=estimate, group=sex)) +
                  geom_line(size=1.5, aes(colour = sex)) +
                  theme(legend.position=c(0.15,0.8), legend.justification=c(1,0)) +
                  xlab("Length") + 
                  ylab("Estimated cap prob")

pred.ggplot

aa=1
xval=20
xsq = xval^2
exp(42.6273520 -4.7076212*xval +  0.1213267*xsq + 0.1755777*aa ) / (1+ exp(42.6273880 -4.7076212*xval +  0.1213267*xsq + 0.1755777*aa ))