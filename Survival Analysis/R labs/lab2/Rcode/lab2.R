#############################################
### Lab 2: Kaplan-Meier survival estimate ###
#############################################

# We use the same data as with the lab 1
# However, we are going to take censoring into account

# (a)
# Import data
nhl1 = read.csv("C:/Applied_Survival_Analysis_Jan2016/lab2/data/nhl1.csv")
nhl1

# KM estimates
library(survival)
nhl.fit = survfit(Surv(nhltime,fail) ~ 1,data = nhl1)

# Add censored = T to include censored event times
summary(nhl.fit,censored = T)

# (a) OPTIONAL
km = data.frame(t = nhl.fit$time, dj = nhl.fit$n.event,rj = nhl.fit$n.risk)
km

km$lambda.hat = km$dj/km$rj
km$surv = cumprod(1-km$lambda.hat)
km

# (b) 95%CI using the "log-log" approach
# Be careful here: the default in R is not the "log-log" approach
nhl.fit = survfit(Surv(nhltime,fail) ~ 1,data = nhl1,conf.type = "log-log")
summary(nhl.fit,censored = T)

setwd("C:/Applied_Survival_Analysis_Jan2016/lab2/graphs")

pdf('lab2KmSurv.pdf',width = 6,height = 6)
plot(nhl.fit,main = "Bone Marrow Transplant for Non-Hodgkin's Lymphoma",
     ylab = "Survival Probability",xlab = "Time to Death or Relapse (months)")
dev.off()

# 95%CI using our code
A = qnorm(0.975)*sqrt( cumsum( km$dj/((km$rj-km$dj)*km$rj) )/log(km$surv)^2 )
km$lower = km$surv^exp(A)
km$upper = km$surv^exp(-A)
round(km,4)

# (c): Median survival time, lower and upper quartile
quantile(nhl.fit,probs = c(0.25,0.5,0.75),conf.int = F)

# (d): Nelson-Aalen estimator of the cumulative hazard function
nhl.na = basehaz(coxph(Surv(nhltime,fail) ~ 1,data = nhl1,ties = "breslow"))
nhl.na

# nhl.na is a data frame with the Nelson-Aalen estimator and the survival times
# Be aware of the following syntax of the plot function
pdf('lab2NA.pdf',width = 6,height = 6)
plot(hazard ~ time,data = nhl.na, type = "b",
     main = "Bone Marrow Transplant for Non-Hodgkin's Lymphoma",
     xlab = "Time to Death or Relapse (months)",xlim = c(0,40),ylim = c(0,1.5),
     ylab = "Nelson-Aalen Cumulative Hazard")
dev.off()

#################################################
# (e): Life table when the data are not grouped #
#################################################

# Import the nurshome data set
# see lecture 1 for a description of the data
nurshome = read.csv("C:/Applied_Survival_Analysis_Jan2016/lab2/data/nurshome.csv")

# Keep only treated 
nurshome = nurshome[nurshome$rx == 1,]
head(nurshome)

# Group length of stay into 100-day intervals
los100 = floor(nurshome$los/100)

# Number of subjects who experienced the event during each interval
died = tapply(nurshome$fail,los100,sum)

# Why not type table(los100[nurshome$fail==1]) ???
# They're supposed to be the same! See
died
table(los100[nurshome$fail==1]) # There is no category for los100 = 10

# Frequency table of length of stay grouped into 100-day intervals
total = table(los100)

# Number of censorings at each interval
censor = total - died

LT.treated = data.frame(LOS100 = unique(los100), died = c(died), censored = c(censor))

#los100 must have one more length than everyone else
los100 = sort(LT.treated$LOS100)
lt = length(los100)
los100[lt+1] = NA

# Life table for the angina example
library(KMsurv)
nursLT = lifetab(los100, sum(total), LT.treated$censored, LT.treated$died)
nursLT

setwd("C:/Applied_Survival_Analysis_Jan2016/lab2/graphs")

#Plot of the survival
pdf("nhact_surv.pdf",height = 5.5,width = 5.5)
plot(los100[-(lt+1)],nursLT[,5], type="s", 
     ylab="Survival", xlab="LOS (100-day intervals)")
dev.off()

#Plot of the hazard
pdf("nhact_haz.pdf",height = 5.5, width = 5.5)
plot(los100[-(lt+1)],nursLT[,7], type="s", 
     ylab="Hazard", xlab="LOS (100-day intervals)")
dev.off()

# (f): R function to compute Kaplan-Meier survival estimates
KM = function(surv,delta)
{
  # surv: vector of survival times
  # delta: vector of failure indicator
  data = data.frame(time = surv,fail = delta)
  
  # Sort the dataset with respect to time
  # Pay special attention to this
  data = data[order(data$time),]
  
  time = data$time
  fail = data$fail
  
  # Distinct failure times
  dist.times = unique(time[fail==1])
  K = length(dist.times)
  
  # Number of events for each failure time
  dj = table(time[fail==1])
  
  # Number at risk for each failure time
  # We present a simple but not computationally efficient way to find it
  rj = rep(NA,K)
  
  for (i in 1:K)
  {
    rj[i] = sum(time >= dist.times[i])
  }
  
  # Hazard estimates
  lambda.hat = dj/rj
  
  # Km estimates
  surv.hat = cumprod(1 - lambda.hat)
  names(surv.hat) = NULL
  
  # Output
  out = data.frame(time = dist.times, survival = surv.hat)
  
  return(round(out,4))
}

# We use seed for reproducibility
set.seed(2)

n = 30
# Simulate the (uncensored) survival times
surv = exp(rnorm(n, mean = - 1, sd = 0.5))

# Simulate the censoring timess
cens = rexp(n, rate = 2)

# Failure indicator
# 1, if survival time < censoring time
# 0, if survival time >= censoring time
delta = 1*(surv < cens)

# Observed survival time
# time = surv if surv < cens
# time = cens if surv >= cens
time = pmin(surv,cens)

# Using survfit
data = data.frame(time = time, fail = delta)
fit = survfit( Surv(time,fail) ~ 1, data = data)
summary(fit)

# Our function
KM(time,delta)