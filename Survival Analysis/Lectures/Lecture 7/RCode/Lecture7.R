##########################################
### LAB 7: Assessing the PH Assumption ###
##########################################

# Import the data
nurshome = read.csv("C:/Applied_Survival_Analysis_Jan2016/lab7/data/nurshome.csv")
head(nurshome)

# Set the working directory
setwd("H:/teaching/Athens/BIO223 (Survival-Yiannoutsos)/lectures/R/lect7")

# (i): Assessing the PH Assumption for gender and marital status separately
# KM estimates by gender  (1=male, 0=female)
library(survival)
pdf("KMgender.pdf",height = 5.5,width = 5.5)
fitKMgen = survfit( Surv(los,fail) ~ gender,data = nurshome)
plot(fitKMgen,mark.time = F,fun = "cloglog",main = "Evaluation of PH Assumption",
     xlab = "Length of Stay (days on a log-scale)",lty = 1:2,col = c("blue","red"),
     ylab = "Ln[-Ln(Survival Probabilities)]")
legend("topleft",lty = 1:2,col = c("blue","red"),bty="n",legend = c("Female","Male"))
dev.off()

# KM estimates by marital status (1=married, 0=not married)
fitKMmar = survfit( Surv(los,fail) ~ married,data = nurshome)

pdf("KMmarried.pdf",height = 5.5,width = 5.5)
plot(fitKMmar,mark.time = F,fun = "cloglog",main = "Evaluation of PH Assumption",
     xlab = "Length of Stay (days on a log-scale)",lty = 1:2,col = c("blue","red"),
     ylab = "Ln[-Ln(Survival Probabilities)]")
legend("topleft",lty = 1:2,col = c("blue","red"),bty="n",
       legend = c("Not married","Married"))
dev.off()

# Compare fitted with observed survival curves for gender
fitCoxgen = survfit( coxph( Surv(los,fail) ~ gender,data = nurshome),
                     newdata = data.frame(gender = c(0,1)) )

# First plot the fitted curves
pdf("ExpVsObs_gen.pdf",height = 5.5,width = 5.5)
plot(fitCoxgen,mark.time = F,xlab = "Length of Stay (days)",ylab = "Survival probability",
     main = "Observed KM vs Predicted Survival Curves \nBy Categories of Gender",
     lty = 1:2,col = c("blue","red"))
# Add the raw curves
lines(fitKMgen,lty = 3:4,col = c("green","orange"),mark.time = F)
legend("topright",bty = "n",lty = 1:4,col = c("blue","red","green","orange"),
       legend = c("Predicted: Female","Predicted: Male",
                  "Observed: Female","Observed: Male"),ncol = 2,cex = 0.9)
dev.off()

# Compare fitted with observed survival curves for married
fitCoxmar = survfit( coxph( Surv(los,fail) ~ married,data = nurshome),
                     newdata = data.frame(married = c(0,1)) )

# First plot the fitted curves
pdf("ExpVsObs_mar.pdf",height = 5.5,width = 5.5)
plot(fitCoxmar,mark.time = F,xlab = "Length of Stay (days)",ylab = "Survival probability",
     main = "Observed KM vs Predicted Survival Curves \nBy Categories of Marital Status",
     lty = 1:2,col = c("blue","red"))
# Add the raw curves
lines(fitKMmar,lty = 3:4,col = c("green","orange"),mark.time = F)
legend("topright",bty = "n",lty = 1:4,col = c("blue","red","green","orange"),
       legend = c("Predicted: Not married","Predicted: Married",
                  "Observed: Not married","Observed: Married"),ncol = 2,cex = 0.9)
dev.off()

# (ii): Generate a new variable
nurshome$hlthsex = c( 1*(nurshome$gender==0 & nurshome$health==2) + 
                      2*(nurshome$gender==1 & nurshome$health==2) +
                      3*(nurshome$gender==0 & nurshome$health==5) +
                      4*(nurshome$gender==1 & nurshome$health==5) )  
table(nurshome$hlthsex)

# Create a new data.frame containing the groups of interest only
nurshome2 = nurshome[nurshome$hlthsex!=0,]
table(nurshome2$hlthsex)

# Decode hlthsex as a factor variable
nurshome2$hlthsex = factor(nurshome2$hlthsex,levels = 1:4,
                           labels = c("Healthier Women","Healthier Men","Sicker Women" 
                                      ,"Sicker Men"))
# Fit KM curves
fitKMgr = survfit( Surv(los,fail) ~ hlthsex,data = nurshome2)

# Let's have a look at the median survival times
quantile(fitKMgr,probs = 0.5,conf.int = F)

pdf("KMhlthsex.pdf",height = 5.5,width = 5.5)
plot(fitKMgr,mark.time = F,fun = "cloglog",ylab = "-Ln[-Ln(Survival Probabilities)]",
     xlab = "Length of Stay (days on a log-scale)",lty = 1:4,
     col = c("blue","red","green","orange"),
     main = "Evaluation of PH Assumption \nBy Categories of health-sex")
legend("topleft",lty = 1:4,col = c("blue","red","green","orange"),bty = "n",
       legend = c("Healthier Women","Healthier Men","Sicker Women","Sicker Men"))
dev.off()

# Fit KM curves for males and females
fitKMgr.females = survfit( Surv(los,fail) ~ hlthsex,data = subset(nurshome2, gender==0))
fitKMgr.males = survfit( Surv(los,fail) ~ hlthsex,data = subset(nurshome2, gender==1))

# Compare fitted with observed survival curves for gender
# Men
fitCoxgen.males = survfit( coxph( Surv(los,fail) ~ hlthsex, data = nurshome2),
                      newdata = data.frame(hlthsex = as.factor(c("Healthier Men","Sicker Men")))) 

pdf("ch12kmph_males.pdf",height = 5.5,width = 5.5)
plot(fitKMgr.males, mark.time = F,
     xlab = "Length of Stay (days)",ylab = "Survival probability",
     lty = 1:2,col = c("blue","red"))
# Add the raw curves
lines(fitCoxgen.males, lty = 1:2,col = c("green","orange"),mark.time = F)
legend("topright",bty = "n",lty = 1:4,col = c("blue","red","green","orange"),
       legend = c("Predicted: Healthier Men","Sicker Men"),ncol = 1,cex = 0.9)
dev.off()

# Women
# Compare fitted with observed survival curves for gender
fitCoxgen.females = survfit( coxph( Surv(los,fail) ~ hlthsex, data = nurshome2),
                      newdata = data.frame(hlthsex = as.factor(c("Healthier Women","Sicker Women"))))

pdf("ch12kmph_females.pdf",height = 5.5,width = 5.5)
plot(fitKMgr.females, mark.time = F,
     xlab = "Length of Stay (days)",ylab = "Survival probability",
     lty = 1:2,col = c("blue","red"))
# Add the raw curves
lines(fitCoxgen.females, lty = 1:2,col = c("green","orange"),mark.time = F)
legend("topright",bty = "n",lty = 1:4,col = c("blue","red","green","orange"),
       legend = c("Predicted: Healthier Women","Sicker Women"),ncol = 1,cex = 0.9)
dev.off()

# (iii): Cox model stratified by gender
# We make the PH assumption for married and health, but NOT for gender
fit.strat <- coxph( Surv(los,fail) ~ married + health + strata(gender),data = nurshome)
summary(fit)

# Here we use the default, i.e., we obtain survival estimates for males and females 
# evaluated at the mean values of married + health.
# You can use any other values through the option newdata
fitCoxgen2 = survfit(fit.strat)

pdf("ch12stratph.stata",height = 5.5,width = 5.5)
plot(fitCoxgen2,mark.time = F,ylab = "-Ln[-Ln(Survival Probabilities)]",fun = "cloglog",
     xlab = "Length of Stay (days on a log-scale)",lty = 1:2,col = c("blue","red"),
     #main = "Evaluation of PH Assumption \nStratification by Gender"
)
legend("topleft",lty = 1:2,col = c("blue","red"),bty = "n",
       legend = c("Female","Male"))
dev.off()
