###################################################
### LAB 6: Model Selection in Survival Analysis ###
###################################################

# (a): Collet’s Approach for Model Selection:
mac <- read.csv("C:/Applied_Survival_Analysis_Jan2016/lab6/data/mac.csv")

# The variables we're interested in are:
vars = c("agecat","sex","cd4","karnof","ivdrug","antiret","rif","clari")

# Note that in this lab we're focusing on time to death
mac[1:5,c("dthtime","dthstat",vars)]

#####################################################################
### Step 1: Fit univariate models to choose candidate predictors. ###
### Use criterion of p <= 0.15 to identify predictors.            ###
#####################################################################

time.expr = "Surv(dthtime,dthstat) ~ "

# See what's happening here
paste(time.expr,vars[1])
formula(paste(time.expr,vars[1]))

# Constructing the Table1
table1 = data.frame(Estimate = rep(NA,8),SE = NA,Pvalue = NA)
rownames(table1) = vars

library(survival)
for (i in 1:(length(vars)-2))
{
  # Fit univariate models
  fit = coxph(formula(paste(time.expr,vars[i])),data = mac)
  
  # Save results
  res = summary(fit)
  
  table1[i,] = c(coef(fit),sqrt(vcov(fit)),res$logtest[3])
}

# We need both rif and clari to evaluate the treatment effect!
fit = coxph(Surv(dthtime,dthstat) ~ rif + clari,data = mac)
res = summary(fit)

table1[7:8,"Estimate"] = coef(fit)
table1[7:8,"SE"] = sqrt(diag(vcov(fit)))
table1[7:8,"Pvalue"] = res$coef[,5]
table1$HR = exp(table1$Estimate)
round(table1,3)

# Wald test for the combined effect of treatment
res$waldtest

######################################################
### Step 2 (i): Fit a multivariate model with all  ###
### significant predictors (p <= 0.15) from Step 1 ###
######################################################

vars[table1$Pvalue<=0.15]

fitStep2i = coxph(Surv(dthtime,dthstat) ~ agecat + sex + cd4 + karnof 
                  + antiret,data = mac)
summary(fitStep2i)

# (ii): Use backward selection to eliminate non-significant predictors 
# in a multivariate  framework using the AIC criterion
# Type install.packages("MASS") if needed
library(MASS)
fitStep2ii = stepAIC(fitStep2i, direction = "backward", k = 3,trace = 10)

##########################################################################
### Step 3: Use forward selection to add any variables not significant ###
### at Step 1 to the multivariate model obtained at the end of Step 2  ###                     
##########################################################################

# Force the variables that are significant at the end of Step 2 into the model,
# by using the scope argument i.e.,
# tell R which is the minimum and maximum possible model to be examined

# To evaluate the signifance of treatment jointly, 
# you must create a factor with the treatment arms
trt = rep(NA,nrow(mac))
trt[mac$rif==1] = "rif"
trt[mac$clari==1] = "clari"
trt[mac$rif==0 & mac$clari==0] = "both"
trt = factor(trt,levels = c("both","rif","clari"))
mac$trt = trt

max.model = formula(Surv(dthtime,dthstat) ~ agecat + sex + cd4 + karnof + 
                      ivdrug + antiret + trt)

fitStep3 = stepAIC(fitStep2i,
                 scope=list(lower = formula(fitStep2i),upper = max.model),
                 k = 3,trace = 10,direction = "forward")

########################################################################
### Step 4: (i) Do final pruning of main-effects model using forward ###        
### forward stepwise regression                                      ###                           
########################################################################

# Starting model (NULL model)
fit = coxph(Surv(dthtime,dthstat) ~ 1,data = mac)

fitStep4i = stepAIC(fit,
                 scope=list(lower = formula(fit),upper = formula(fitStep2i)),
                 k = 3,trace = 10,direction = "both")

# (ii) Adding interaction terms
mac$agsex = mac$agecat*mac$sex 
mac$agcd4 = mac$agecat*mac$cd4 
mac$agkar = mac$agecat*mac$karnof 
mac$aganti = mac$agecat*mac$antiret 
mac$sexcd4 = mac$sex*mac$cd4 
mac$sexkar = mac$sex*mac$karnof 
mac$sexanti = mac$sex*mac$antiret 
mac$cd4kar = mac$cd4*mac$karnof 
mac$cd4anti = mac$cd4*mac$antiret 
mac$karanti = mac$karnof*mac$antiret

max.model = formula(Surv(dthtime,dthstat) ~ agecat + sex + cd4 + karnof 
                    + antiret + agsex 
                    + agcd4 + agkar + aganti + sexcd4 + sexkar + sexanti 
                    + cd4kar + cd4anti + karanti)

# Fit the maximum model
fit = coxph(max.model,data = mac)

fitStep4ii = stepAIC(fit,
                     scope = list(lower = formula(fitStep2i),upper = max.model),
                   k = 3,trace = 0,direction = "both")
fitStep4ii$anova

##############################################
### Step 5: alternate coding of covariates ###
##############################################

# See the formula of the last model
formula(fitStep4ii)

# (i): cd4cat instead of cd4
fitStep5i = coxph(Surv(dthtime, dthstat) ~ agecat + sex + cd4cat + 
                    karnof + antiret + sexanti + karanti,data = mac)

# (ii): age instead of agecat
fitStep5ii = coxph(Surv(dthtime, dthstat) ~ age + sex + cd4 + karnof 
                   + antiret + sexanti + karanti,data = mac)

# Summarize the results in a table
models = data.frame(Model = c("Step 2 (i)","Step 2 (ii)","Step 3","Step 4 (i)",
                              "Step 4 (ii)","Step 5 (i)","Step 5 (ii)"),
                    Covariates = NA,twlogl = NA,q = NA,AIC = NA)

md = c("fitStep2i","fitStep2ii","fitStep3","fitStep4i",
       "fitStep4ii","fitStep5i","fitStep5ii")

for (i in 1:7)
{
  obj = get(md[i])
  models[i,2] = paste(names(coef(obj)),collapse = " ")
  models[i,-(1:2)] = c(-2*obj$logl[2],extractAIC(obj,k=3))
}

models

#################################################
### Assessing the overall fit of the model by ###
### checking the residuals                    ###
#################################################

fit = coxph(Surv(dthtime, dthstat) ~ age + sex + cd4 + karnof + antiret,data = mac)
#summary(fit)

# (b)-(i): Generalized (Cox-Snell) Residuals
mac$csres = mac$dthstat - residuals(fit,type = "martingale")

setwd("C:/Applied_Survival_Analysis_Jan2016/lab6/graphs")

pdf("Surv_coxsnell.pdf",height = 5,width = 5)
fitcs = survfit(Surv(csres,dthstat)~1,data = mac)
plot(fitcs,fun = "cloglog",conf.int = F,mark.time = F, 
     xlab = "Cox-Snell residuals (log scale)",ylab = "Log(-log Scs)")
abline(h = 0,col = "red")
abline(v = 1,col = "red") # Don't forget the log scale!
dev.off()

# b-(ii): Martingale Residuals
mac$mg = residuals(fit,type = "martingale")
mac$betaz = predict(fit,type = "lp")

pdf("mgale.pdf",height = 5,width = 5)
plot(mg ~ betaz,data = mac,xlab = "Linear prediction",
     ylab = "Martingale residuals")
dev.off()

# b-(iii): Deviance residuals
mac$devres = residuals(fit,type = "deviance")

# Deviance residuals vs linear prediction and other covariates
pdf("devres.pdf",height = 5.5,width = 5.5)
par(mfrow = c(2,2))
plot(devres ~ betaz,data = mac,xlab = "Linear prediction",
     ylab = "Deviance residuals")
abline(h = 0,col = "red")

plot(devres ~ age,data = mac,xlab = "Age (years)",
     ylab = "Deviance residuals")
abline(h = 0,col = "red")

plot(devres ~ cd4,data = mac,xlab = "CD4+ cell count",
     ylab = "Deviance residuals")
abline(h = 0,col = "red")

plot(devres ~ karnof,data = mac,xlab = "Karnofsky score status",
     ylab = "Deviance residuals")
abline(h = 0,col = "red")
dev.off()

# (b)-(iv): Weighted Schoenfeld Residuals
temp = cox.zph(fit,transform = "identity")

# Age
pdf("sclSchAge.pdf",height = 5.5,width = 5.5)
plot(temp[1],main = "scaled Schoenfeld residuals of Age")
dev.off()

# CD4
pdf("sclSchCd4.pdf",height = 5.5,width = 5.5)
plot(temp[3],main ="scaled Schoenfeld residuals of CD4")
dev.off()

################################################
# (b)-(v): Checking the functional form of CD4 #
################################################

# Fit a model without CD4
coxNoCD4 = coxph(Surv(dthtime, dthstat) ~ age + sex + karnof + antiret,data = mac)
# summary(coxNoCD4)
mac$mgNoCD4 = residuals(coxNoCD4,type = "martingale")

# Regress CD4 on the other covariates
lmCD4 = lm(cd4 ~ age + sex + karnof + antiret,data = mac)
# summary(lmCD4)
mac$CD4res = lmCD4$res

pdf("funCD4.pdf",height = 5.0,width = 5.0)
plot(mgNoCD4 ~ CD4res,data = mac, ylab = "CD4-free martingale residuals",
     xlab = "CD4 residuals")
lines(smooth.spline(mac$CD4res, mac$mgNoCD4, df = 4), col = "red", lwd = 2)
dev.off()

# Try sqrt(CD4), very popular transformation in longitudinal data analysis!
mac$sqCD4res = lm(I(sqrt(cd4)) ~ age + sex + karnof + antiret,data = mac)$res

pdf("funsqCD4.pdf",height = 5.0,width = 5.0)
plot(mgNoCD4 ~ sqCD4res,data = mac, ylab = "sqrt(CD4)-free martingale residuals",
     xlab = "sqrt(CD4) residuals")
lines(smooth.spline(mac$sqCD4res, mac$mgNoCD4, df = 4), col = "red", lwd = 2)
dev.off()
