#########################
### Interactions in R ###
#########################

mac <- read.csv("C:/Applied_Survival_Analysis_Jan2016/lab6/data/mac.csv")
head(mac)

# Trt as a factor
trt = rep(NA,nrow(mac))
trt[mac$rif==1] = "rif"
trt[mac$clari==1] = "clari"
trt[mac$rif==0 & mac$clari==0] = "both"
trt = factor(trt,levels = c("both","rif","clari"))
mac$trt = trt

library(survival)
fitmin = coxph(Surv(dthtime,dthstat) ~ trt + cd4,data = mac)
summary(fitmin)

fitmax = coxph(Surv(dthtime,dthstat) ~ trt*cd4,data = mac)
summary(fitmax)

library(MASS)
stepAIC(fitmin,
        scope = list(lower = formula(fitmin),upper = formula(fitmax)),
        k = 3,trace = 10,direction = "forward")
