#########################################################
### Logrank test as a stratified Mantel-Haenszel test ###
#########################################################

leukem = read.csv("C:/Applied_Survival_Analysis_Jan2016/lab3/data/leukem.csv")
head(leukem)

# Unique death times
times.uni = sort(unique(leukem$weeks[leukem$remiss==1]))
K = length(times.uni)

# We're going to get the number of deaths and the number at risk by group
# using the options of summary function
library(survival)
fit = survfit(Surv(weeks,remiss) ~ trt, data = leukem)

# Please, use plausible times to get a meaningful result
sumFit = summary(fit,times = times.uni,extend = T) 
sumFit

# Number of deaths and number at risk in the control group
dj.grp0 = sumFit$n.event[sumFit$strata=="trt=0"]
rj.grp0 = sumFit$n.risk[sumFit$strata=="trt=0"]

# Number of deaths and number at risk in the treatment group
dj.grp1 = sumFit$n.event[sumFit$strata=="trt=1"]
rj.grp1 = sumFit$n.risk[sumFit$strata=="trt=1"]

# To perform a Mantel-Haenszel analysis 
# we need to save our data into an a 2*2*K matrix
# i.e., a series of 2*2 group*failure contingency tables for each death time

# see the results in a table
cbind(DeathContr = dj.grp0,AliveContr = rj.grp0-dj.grp0,
      DeathTrt = dj.grp1,AliveContr = rj.grp1-dj.grp1)

MHarray = array(t(cbind(dj.grp0,dj.grp1,rj.grp0-dj.grp0,rj.grp1-dj.grp1)),dim = c(2,2,K),
                dimnames = list(Group = c("Control","6-MP"),Death = c("Yes","No"),
                                Time = times.uni))

# and in an array (for example see the first 5 contingency tables)
MHarray[,,1:5]

fitMH = mantelhaen.test(MHarray,correct = F)
fitMH

# Compare
library(survival)
fitlogR = survdiff(Surv(weeks,remiss) ~ trt,data = leukem)
fitlogR

fitMH$statistic - fitlogR$chisq