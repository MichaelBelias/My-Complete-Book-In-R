#####################################################################
### R notes for Lecture 5: More on Cox Proportional Hazards Model ###
#####################################################################

# MAC dataset, here we're interested in time to MAC disease
# NOT time to death
mac = read.csv("h:/teaching/athens/BIO223 (Survival-Yiannoutsos)/data/mac.csv")

# See some lines of the data
mac[1:10,c("patid","macstat","mactime","karnof","rif","clari","cd4")]

# Number of events
table(mac$macstat)
summary(mac$mactime)

# Delete 26 patients with zero time
table(mac$mactime==0)
mac = mac[mac$mactime!=0,]

# (a)-(i): Time to MAC disease in relation
# to karnofsky score and treatment group
library(survival)
fit.mac1 = coxph( Surv(mactime,macstat) ~ karnof + rif + clari,data = mac)
summary(fit.mac1)

# a-(ii): Contruncting 95% CIs for the hazard ratios
# First calculate CI's for log(HR) and 
# then back-transform to the HR scale
exp(confint(fit.mac1))

# or,
exp( cbind( coef(fit.mac1)-qnorm(0.975)*sqrt(diag(vcov(fit.mac1))),
       coef(fit.mac1),
       coef(fit.mac1)+qnorm(0.975)*sqrt(diag(vcov(fit.mac1))) ))

# a-(iii): Wald test for the karnofsky score
z.tests = coef(fit.mac1)/sqrt(diag(vcov(fit.mac1)))
chi2 = z.tests^2
chi2[1]

# a-(iv): Add CD4 cell count in the model
fit.mac2 = coxph( Surv(mactime,macstat) ~ karnof + rif + clari + cd4,data = mac)
summary(fit.mac2)

# LR test comparing the model with and without CD4
anova(fit.mac1,fit.mac2)

# or,
-2*(fit.mac1$loglik - fit.mac2$loglik)[2]

# a-(v): Test for an overall treatment effect using a (multivariate) Wald test
# after taking into account the Karnofsky score and CD4 count
fit.mac2 = coxph( Surv(mactime,macstat) ~ karnof + rif + clari + cd4,data = mac)
coef(fit.mac2)

# Please, have a look at the GLM notes
# We are contructing a chi-square test with 2 degrees of freedon
chi2 = t(coef(fit.mac2)[2:3]) %*% solve(vcov(fit.mac2)[2:3,2:3]) %*% coef(fit.mac2)[2:3]
chi2

#Or, alternatively:
library(aod)
chi2<-wald.test(b=coef(fit.mac2),Sigma=vcov(fit.mac2), Terms=2:3)
chi2
# p-value
1-pchisq(chi2,df=2)

# a-(vi): To test whether there is a DIFFERENCE between the 
# rif and clari treatment arms, 
# see the coefficient of rif 
fit.mac3 = coxph( Surv(mactime,macstat) ~ karnof + rif + I(rif + clari) + cd4,data = mac)
summary(fit.mac3)

# Wald test of the difference betweeen the treatment arms
L<-matrix(c(0, -1, 1, 0), nrow=1, ncol=4)
wald.test(b=coef(fit.mac2), Sigma=vcov(fit.mac3), L=L)
##################################################################
### Survival Function and Predicted Medians: Nursing Home Data ###
##################################################################
nurshome = read.csv("h:/teaching/athens/BIO223 (Survival-Yiannoutsos)/data/nurshome.csv")
nurshome[1:4,]

# Two-way frequency table
tab = table(nurshome$married,nurshome$health)
tab

# To see the percentages (by row and by col)
prop.table(tab,1)
prop.table(tab,2)

# b-(i) : Cox PH model
library(survival)
fit.cox = coxph( Surv(los,fail) ~ married + health, data = nurshome)
summary(fit.cox)

# b-(ii): Calculation of median by a KM approach
# First, create the grouping factor
nurshome$group = c(1*(nurshome$mar==0 & nurshome$health==2) +
                   2*(nurshome$mar==0 & nurshome$health==5) +
                   3*(nurshome$mar==1 & nurshome$health==2) +
                   4*(nurshome$mar==1 & nurshome$health==5))

# We're not interested in the category of group = 0
table(nurshome$group)

# We can create a new data frame excluding group = 0
# Of course, several other options are available to deal with issue ...
nurshome2 = nurshome[nurshome$group!=0,]

# Specify the order of the values and also label the values
nurshome2$group = factor(nurshome2$group,levels = 1:4,labels = c("Single, healthy",
                                                                 "Single, unhealthy",
                                                                 "Married, healthy",
                                                                 "Married, unhealthy"))

fit.KM = survfit( Surv(los,fail) ~ group,data = nurshome2)
fit.KM

# b-(iii) We calculate the medians assuming proportional hazards
newdata = data.frame(married = c(0,0,1,1),health = c(2,5,2,5))
surv.cox = survfit(fit.cox, newdata = newdata)
surv.cox

# b-(iv)  Calculate survival for groups which do not exist in the data
newdata2 = data.frame(married = c(0,0,1,1),health = c(0, 2, 0 ,2))
surv.cox2 = survfit(fit.cox, newdata = newdata2)
summary(surv.cox2)

setwd("h:/teaching/athens/BIO223 (Survival-Yiannoutsos)/lectures/R/lect5/")

pdf("nurshSurvEst2.pdf",height = 6,width = 6)
plot(surv.cox2 ,mark.time = F,xlab = "Length of Stay of Resident (days)",
     ylab = "Survival probability",lty = 1:4,
     col = c("blue","red","green","orange"))
legend("topright",lty = 1:4,col = c("blue","red","green","orange"),bty = "n",
       legend = c("Single, healthy","Single, unhealthy","Married, healthy",
                  "Married, unhealthy"))
dev.off()

# b-(iv) Breslow estimate of baseline cumulative hazard
# Design matrix of covariates (Don't include an intercept!)
x = cbind(nurshome$married,nurshome$health)
beta.hat = coef(fit.cox)

brescHaz = function(surv,fail,x,beta)
{
  # Sort data with respect to time
  data = cbind(surv,fail,x)
  data = data[order(surv),]
  
  time = data[,1]
  delta = data[,2]
  x = data[,3:ncol(data)]
  
  # Distinct event times
  times.uni = unique(time[delta==1])
  K = length(times.uni)
  n = length(time)
  
  # Number of failures for each event time
  dj = table(time[delta==1])
  
  # Find where each risk set starts
  ind = match(times.uni,time)
  
  # Linear predictor,
  # %*% denotes matrix multiplication
  Xbeta = c(x %*% beta)
  eXbeta = exp(Xbeta)
  
  # Breslow estimate of baseline cumulative hazard
  BrcumHaz = rep(NA,K)
  
  for (i in 1:K)
  {
    BrcumHaz[i] = dj[i]/sum(eXbeta[ind[i]:n])
  }
  BrcumHaz = cumsum(BrcumHaz)
  
  out = data.frame(time = times.uni, surv0 = exp(-BrcumHaz))
  return(out)
}

# Our estimation
fit = brescHaz(nurshome$los,nurshome$fail,x,beta.hat)

# Using survfit
fit2 = survfit(fit.cox,newdata = data.frame(married = 0,health = 0),type = "aalen")

# Compare
fit[1:10,]
cbind(summary(fit2)$time,summary(fit2)$surv)[1:10,]

################################ P-year survival ###########################
newdata3 <- data.frame(married = c(0,0,1,1),
                       health = c(2,5,2,5), 
                       fail= c(1,1,1,1), los=365)

predict(fit.cox, newdata=newdata3, type="expected")