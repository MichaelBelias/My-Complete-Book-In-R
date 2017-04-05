#############################################
### LAB 4: Cox Proportional Hazards Model ###
#############################################

# Fitting Cox Model and handling of ties: Nursing Home Data
# Import data
nurshome = read.csv("C:/Applied_Survival_Analysis_Jan2016/lab4/data/nurshome.csv")

# See few lines of the data
nurshome[1:5,]

# (a): The Efron method is the default!
library(survival)
fit = coxph( Surv(los,fail) ~ married,data = nurshome)
summary(fit)

# (b): Methods for tie handling
# Breslow method
fit.breslow = coxph( Surv(los,fail) ~ married,data = nurshome, ties = "breslow")
summary(fit.breslow)

# Efron method
fit.efron = coxph( Surv(los,fail) ~ married,data = nurshome, ties = "efron")
summary(fit.efron)

# Exact (discrete)
fit.exact = coxph( Surv(los,fail) ~ married,data = nurshome, ties = "exact")
summary(fit.exact)

# (c): Logrank and wilcoxon test
# Logrank
survdiff(Surv(los,fail) ~ married,data = nurshome)
chi2.logrank = survdiff(Surv(los,fail) ~ married,data = nurshome)$chisq
chi2.logrank

# Peto & Peto modification of the Gehan-Wilcoxon test
survdiff(Surv(los,fail) ~ married,data = nurshome, rho = 1)

# The logrank test is equivalent to 
# the cox model using the "exact" option
c(fit.exact$score,chi2.logrank)

fit.exact$score - chi2.logrank

# Simulated data from a clinical trial
# Treament group: T ~ exp(rate = 2)
# Control group: T ~ exp(rate = 4)
set.seed(5)
n = 5000

# These are uncensored survival times
timeTrt = rexp(n,rate = 2)
timeCont = rexp(n,rate = 4)
surv = c(timeTrt,timeCont)

# Censoring time
# Mean = 1/rate for exponentially distributed times
cens = rexp(2*n,rate = 1/5)

# Failure indicator
delta = 1*(surv < cens) 

# Observed survival time
time = pmin(surv,cens)

# Save data in a data frame
data = data.frame(time = time,fail = delta, trt = c(rep(1,n),rep(0,n) ))

# Fit a cox model
fit.cox = coxph(Surv(time,fail) ~ trt,data = data)
summary(fit.cox)


