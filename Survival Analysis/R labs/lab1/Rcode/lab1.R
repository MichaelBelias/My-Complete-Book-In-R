#########################################################################
### Time to death or relapse (in months) after bone marrow transplant ###
### for 12 patients with non-Hodgkin's lymphoma                       ### 
#########################################################################

# (a)
# Time to death or relapse
time = c(1,2,2,2,3,5,6,7,8,16,17,34)

# Create a data frame with the survival time, 
# the number of deaths at each time point (assuming no censoring)
# and the empirical survival estimate
surv.table = data.frame(time = unique(time), failures = c(table(time)) )

# Number of patients being alive
surv.table$alive = 12 - cumsum(surv.table$failures)

# Survival probability
surv.table$surv = surv.table$alive/12
surv.table

# (b)
# Don't forget to import library survival
library(survival)
nhl.data = data.frame(time = time, delta = 1)

# Use the Surv function to declare survival data
nhl.fit = survfit(Surv(time,delta) ~ 1,data = nhl.data)
summary(nhl.fit)

# (c)
# Set working directory
setwd("C:/Applied_Survival_Analysis_Jan2016/lab1/graphs")

pdf('lab1EmpSurv.pdf',width = 5,height = 5)
plot(nhl.fit,conf.int = F,mark.time = F,xlab = "Time to death or relapse",
     ylab = "Survival probability")
dev.off()

# (d)
# Calculation of standard errors
# Just the SE of a binomial proportion, p(1-p)/n
surv.table$se = sqrt(surv.table$surv*(1-surv.table$surv)/12)
surv.table[,c("time","se")]

# (e)
# Median survival
quantile(nhl.fit,probs = c(0.25,0.50,0.75),conf.int = F)



