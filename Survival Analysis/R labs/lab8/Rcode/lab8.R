###########################################
### LAB 8: Parametric Survival Analysis ###
###########################################

# Import the data into R
nurshome = read.csv("C:/Applied_Survival_Analysis_Jan2016/lab7/data/nurshome.csv")
nurshome$losyr = nurshome$los/365
head(nurshome)

# (i): Fit an exponential model
library(survival)
fitExp = survreg( Surv(losyr,fail) ~ gender,data = nurshome,dist = "exp")
summary(fitExp)

# corresponding cox model
fitCox = coxph( Surv(losyr,fail) ~ gender,data = nurshome)
summary(fitCox)

# Tranform coefficients into a proportional hazards parameterization
- coef(fitExp)

# Get the KM estimates
fitKM = survfit( Surv(losyr,fail) ~ gender,data = nurshome)

# (ii): Comparing exponential model with KM
setwd("C:/Applied_Survival_Analysis_Jan2016/lab8/graphs")

pdf("ExpKM.pdf",height = 5.5, width = 5.5)
plot(fitKM,mark.time = F,main = "Predicted survival for Exponential model VS KM",
     xlab = "Length of Stay (years)",ylab = "Survival probability",lty = 1:2,
     col = c("blue","red"))
# Add the curves fitted by the exponential model
pct = seq(0,1,by = 0.001)
lines(predict(fitExp,newdata = data.frame(gender = 0),type = "quantile",
              p = pct),1-pct,lty = 3,col = "green")
lines(predict(fitExp,newdata = data.frame(gender = 1),type = "quantile",
              p = pct),1-pct,lty = 4,col = "orange")
legend("topright",bty = "n",lty = 1:4,col = c("blue","red","green","orange"),
       legend = c("KM: Females","KM: Males","EXP: Females","EXP: Males"),ncol = 2)
dev.off()

# (iii): Fit a Weibull model
fitWei = survreg( Surv(losyr,fail) ~ gender,data = nurshome,dist = "weibull")
summary(fitWei)

# Tranform coefficients into a proportional hazards parameterization
- coef(fitWei)/fitWei$scale

# Comparing Weibull model with KM
pdf("WeiKM.pdf",height = 5.5, width = 5.5)
plot(fitKM,mark.time = F,main = "Predicted survival for Weibull model VS KM",
     xlab = "Length of Stay (years)",ylab = "Survival probability",lty = 1:2,
     col = c("blue","red"))
# Add curves fitted by the Weibull model
pct = seq(0,1,by = 0.001)
lines(predict(fitWei,newdata = data.frame(gender = 0),type = "quantile",
              p = pct),1-pct,lty = 3,col = "green")
lines(predict(fitWei,newdata = data.frame(gender = 1),type = "quantile",
              p = pct),1-pct,lty = 4,col = "orange")
legend("topright",bty = "n",lty = 1:4,col = c("blue","red","green","orange"),
       legend = c("KM: Females","KM: Males","WEI: Females","WEI: Males"),ncol = 2)
dev.off()

# (iv): log-cumulative hazard plot
pdf("logChaz.pdf",height = 5.5, width = 5.5)
plot(fitKM,mark.time = F,main = "Estimated log cumulative hazard vs log(time)", 
     xlab = "Length of Stay (years on log-scale)",ylab = "Log cumulative hazard",lty = 1:2,
     col = c("blue","red"), fun = "cloglog")
dev.off()

# (v): Mean and median times
ExpModel = data.frame(Mean = rep(NA,2),Median = rep(NA,2)) 
rownames(ExpModel) = c("Females","Males")

# Coefficients in PH form
beta = - coef(fitExp)

ExpModel[,1] = 1/exp(cumsum(beta))
ExpModel[,2] = - log(0.5)/exp(cumsum(beta))
ExpModel

# Weibull model
WeiModel = data.frame(Mean = rep(NA,2),Median = rep(NA,2)) 
rownames(WeiModel) = c("Females","Males")

# Coefficients in PH form
beta = - coef(fitWei)/fitWei$scale
kappa = 1/fitWei$scale

WeiModel[,1] = gamma(1/kappa + 1)*exp(cumsum(beta))^(-1/kappa)
WeiModel[,2] = (- log(0.5)/exp(cumsum(beta)) )^(1/kappa)
WeiModel

# (vi): Computing the loglikelihood for the Weibull model
f = function(param)
{
  # x: design matrix including the intercept
  p = ncol(x)
  
  # Param: The parameter vector
  beta = param[1:p]
  kappa = param[p+1]
  
  lambda = exp(x%*%beta)
  
  logl = sum(delta*(log(kappa) + log(lambda) + (kappa-1)*log(time)) - lambda*time^kappa)
  
  return(-logl)
}

# Fit the weibull model
fitWei = survreg( Surv(losyr,fail) ~ gender,data = nurshome,dist = "weibull")

# Tranform coefficients into a proportional hazards parameterization
beta = - coef(fitWei)/fitWei$scale
kappa = 1/fitWei$scale
x = cbind(1,nurshome$gender)
time = nurshome$losyr
delta = nurshome$fail

# Compare
f(c(beta,kappa))
fitWei$loglik[2]

f(c(beta,kappa)) + fitWei$loglik[2]

# Use the nlm function to maximize the loglikelihood
# Starting values
starting = c(rep(0,2),1)

nlm(f,starting,gradtol = 1e-10)
beta;kappa