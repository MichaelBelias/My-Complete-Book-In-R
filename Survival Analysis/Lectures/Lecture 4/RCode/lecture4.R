##########################################################
### Code for lecture 4: Cox Proportional Hazards Model ###
##########################################################
library(survival)
# Fitting Cox Model and handling of ties: Nursing Home Data
# Import data
leukemia = read.csv("h:/teaching/athens/BIO223 (Survival-Yiannoutsos)/data/leukem.csv")
leukemia.fit = coxph( Surv(weeks,remiss) ~ trt,data = leukemia)
summary(leukemia.fit)

# Fecundability example
#import data
fecund<- read.csv("h:/teaching/athens/BIO223 (Survival-Yiannoutsos)/data/fecundability.csv")
fecund.fit = coxph( Surv(cycle,censor) ~ smoker, data = fecund)
summary(fecund.fit)


