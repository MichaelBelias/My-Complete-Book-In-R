##########################
### Lecture 2 - R Code ###
##########################

# leukemia data
leukem = read.csv("C:/Users/christofer/Dropbox/Applied_Survival_Analysis_Jan2016/Data/leukem.csv")
head(leukem)

# Emprical survival estimate of control arm, assuming no censoring
library(survival)
empCtl = survfit(Surv(weeks) ~ 1,subset = trt == 0,data = leukem,conf.type = "log-log")
plot(empCtl,ylab = "Estimated survival probability",conf.int = T,
     xlab = "Time from Remission to Relapse (weeks)",main = "Empirical Survivor Function (Control arm)")

# KM plot for treated
KMtreat = survfit(Surv(weeks,remiss) ~ 1, subset = trt == 1,data = leukem,conf.type = "log-log")
plot(KMtreat,mark.time = F,main = "Kaplan-Meier Survival estimate (6MP arm)",
     ylab = "Estimated survival probability",conf.int = T,
     xlab = "Time from Remission to Relapse (weeks)")

# To get R output
summary(KMtreat,censored = T) # To match with the command sts li in STATA

# To obtain quantiles for treated
quantile(KMtreat,prob = c(0.25,0.5,0.75),conf.int = F)

###################################################
### Nelson-Aalen estimator of cumulative hazard ###
### and Fleming-Harrington of survival          ###
###################################################

# We consider the nurshome dataset
nurshome <- read.csv("C:/Users/christofer/Dropbox/Applied_Survival_Analysis_Jan2016/Data/nurshome.csv")

# A subset of the data
nurshome = nurshome[nurshome$rx==0 & nurshome$gen==0 
                    & nurshome$health==2 & nurshome$mar==1,]

# We can get the Nelson-Aalen estimator of cumulative hazard in R
# (i): using a Cox model, 
# since the breslow estimator of baseline cumulative hazard reduces to the Nelson-Aalen estimator
fit = basehaz(coxph(Surv(los,fail) ~ 1,data = nurshome,ties = "breslow"))

# Comparing Fleming-Harrington estimate with KM estimate
fit$skm = survfit(Surv(los,fail) ~ 1,data = nurshome)$surv
fit$sfh = exp(-fit$haz)
fit

