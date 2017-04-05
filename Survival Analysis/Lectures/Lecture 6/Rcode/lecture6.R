################################################################
### Code for lecture 6: Model Selection in Survival Analysis ###
################################################################

halibut = read.table("h:/teaching/athens/BIO223 (Survival-Yiannoutsos)/data/halibut.dat", quote="\"")
names(halibut) = c("obs","survtime","censor","towdur","depth","length","handling","logcatch")
head(halibut)

setwd("h:/teaching/athens/BIO223 (Survival-Yiannoutsos)/lectures/R/lect6/")
# Univariable analyses
KM.towdur<-survfit(Surv(survtime, censor)~as.numeric(towdur<100), 
                   data=halibut)
pdf("towdur.pdf",height = 6,width = 6)
plot(KM.towdur, lty=1:2, col=c("blue", "red"), 
     xlab="Minutes since capture", ylab='Survival probability')
     legend("topright", lty = 1:2, col = c("blue","red"), bty = "n",
            legend = c("Above 100 mins", "Below 100 mins"))
dev.off()
KM.length<-survfit(Surv(survtime, censor)~as.numeric(length<43), 
                        data=halibut)
pdf("length.pdf",height = 6,width = 6)
     plot(KM.length,  lty=1:2, col=c("blue", "red"), 
          xlab="Minutes since capture", ylab='Survival probability')
     legend("topright", lty = 1:2, col = c("blue","red"), bty = "n",
            legend = c("Above 43 cm", "Below 43 cm"))
dev.off()

KM.depth<-survfit(Surv(survtime, censor)~as.numeric(depth<8), 
                   data=halibut)
pdf("depth.pdf",height = 6,width = 6)
plot(KM.depth, lty=1:2, col=c("blue", "red"), 
     xlab="Minutes since capture", ylab='Survival probability')
legend("topright", lty = 1:2, col = c("blue","red"), bty = "n",
       legend = c("> 8 meters", "< 8 meters"))
dev.off()

KM.handling<-survfit(Surv(survtime, censor)~as.numeric(handling<10), 
                   data=halibut)
pdf("handling.pdf",height = 6,width = 6)
plot(KM.handling,  lty=1:2, col=c("blue", "red"), 
     xlab="Minutes since capture", ylab='Survival probability')
legend("topright", lty = 1:2, col = c("blue","red"), bty = "n",
       legend = c("More than 10 mins", "Less than 10 mins"))
dev.off()

KM.logcatch<-survfit(Surv(survtime, censor)~as.numeric(logcatch<4.203), 
                     data=halibut)
pdf("logcatch.pdf",height = 6,width = 6)
plot(KM.logcatch,  lty=1:2, col=c("blue", "red"), 
     xlab="Minutes since capture", ylab='Survival probability')
legend("topright", lty = 1:2, col = c("blue","red"), bty = "n",
       legend = c("Above the median", "Below the median"))
dev.off()

# The stepAIC function carries out stepwise model selection according to the AIC criterion

# Forward selection
library(survival)
fitmin = coxph(Surv(survtime,censor) ~ 1, data = halibut)
fitmax = coxph(Surv(survtime,censor) ~ towdur + depth + length + handling + logcatch
               ,data = halibut)
library(MASS)
fit.forward = stepAIC(fitmin,scope = list(lower = formula(fitmin),upper = formula(fitmax)),
                      direction = "forward",k = 3,trace = 1)
summary(fit.forward)

# Backward Selection
fit.backward = stepAIC(fitmax,scope = list(lower = formula(fitmin),
                                           upper = formula(fitmax)),
                      direction = "backward", k = 3,trace = 1)
summary(fit.backward)

# Stepwise Backward Selection 
fit.stepBack = stepAIC(fitmax,scope = list(lower = formula(fitmin),
                                           upper = formula(fitmax)),
                       direction = "both",k = 3,trace = 1)
summary(fit.stepBack)

# Stepwise Forward Selection 
fit.stepForw = stepAIC(fitmin,scope = list(lower = formula(fitmin),
                                           upper = formula(fitmax)),
                       direction = "both",k = 3,trace = 1)
summary(fit.stepForw)

#################################################
### Model checking by examining the residuals ###
#################################################

# Fit the model to be checked
fit = coxph( Surv(survtime,censor) ~ towdur + handling + length + logcatch, data = halibut)
summary(fit)

# Cox snell residuals
halibut$csres = halibut$censor - residuals(fit,type = "martingale")

fitcs = survfit(Surv(csres,censor) ~ 1,data = halibut)

pdf("coxsnell.pdf",height = 6,width = 6)
par(mfrow=c(1,2))
plot(fitcs, fun = "cloglog", conf.int = F, mark.time = F,
     xlab = "Cox-Snell residuals",ylab = "Log(-log Scs)")
abline(h = 0,col = "red")
abline(v = 1,col = "red") 
plot(log(fitcs$time), log(-log(fitcs$surv)), type="l",
     xlab = "Cox-Snell residuals (log scale)",ylab = "Log(-log Scs)")
abline(h = 0,col = "red")
abline(v = 0,col = "red") 
dev.off()
# Available residuals for a coxph object, see
?residuals.coxph

# Martingale residuals
par(mfrow=c(2,))
plot(halibut$towdur, residuals(fit, type="martingale"),
     xlim=c(0,120), ylab="Martingale residual", xlab="Towing duration")
abline(h = 0,col = "red")

plot(halibut$length, residuals(fit, type="martingale"),
     xlim=c(20,60), ylab="Martingale residual", xlab="Length of fish")
abline(h = 0,col = "red")

plot(halibut$logcatch, residuals(fit, type="martingale"),
     xlim=c(2,10), ylab="Martingale residual", xlab="Log(catch)")
abline(h = 0,col = "red")

plot(halibut$handling, residuals(fit, type="martingale"),
     xlim=c(0,40), ylab="Martingale residual", xlab="Handling duration")
abline(h = 0,col = "red")

par(mfrow=c(1,1))
plot(fit$linear.predictors, residuals(fit, type="martingale"),
     xlim=c(-2,2), ylab="Martingale residual", xlab="Predicted log(HR)")
abline(h = 0,col = "red")

# Deviance residuals
par(mfrow=c(2,2))
plot(halibut$towdur, residuals(fit, type="deviance"),
     xlim=c(0,120), ylab="Deviance residual", xlab="Towing duration")
abline(h = 0,col = "red")

plot(halibut$length, residuals(fit, type="deviance"),
     xlim=c(20,60), ylab="Deviance residual", xlab="Length of fish")
abline(h = 0,col = "red")

plot(halibut$logcatch, residuals(fit, type="deviance"),
     xlim=c(2,10), ylab="Deviance residual", xlab="Log(catch)")
abline(h = 0,col = "red")

plot(halibut$handling, residuals(fit, type="deviance"),
     xlim=c(0,40), ylab="Deviance residual", xlab="Handling duration")
abline(h = 0,col = "red")

par(mfrow=c(1,1))
plot(fit$linear.predictors, residuals(fit, type="deviance"),
     xlim=c(-2,2), ylab="Deviance residual", xlab="Predicted log(HR)")
abline(h = 0,col = "red")

# Schoenfeld residuals
par(mfrow=c(2,2))
plot(halibut$towdur[halibut$censor==1], residuals(fit, type="schoenfeld")[,1],
     xlim=c(0,120), ylab="Schoenfeld residual", xlab="Towing duration")
abline(h = 0,col = "red")

plot(halibut$length[halibut$censor==1], residuals(fit, type="schoenfeld")[,3],
     xlim=c(20,60), ylab="Schoenfeld residual", xlab="Length of fish")
abline(h = 0,col = "red")

plot(halibut$logcatch[halibut$censor==1], residuals(fit, type="schoenfeld")[,4],
     xlim=c(2,10), ylab="Schoenfeld residual", xlab="Log(catch)")
abline(h = 0,col = "red")

plot(halibut$handling[halibut$censor==1], residuals(fit, type="schoenfeld")[,2],
     xlim=c(0,40), ylab="Schoenfeld residual", xlab="Handling duration")
abline(h = 0,col = "red")

par(mfrow=c(1,1))
plot(fit$linear.predictors, residuals(fit, type="deviance"),
     xlim=c(-2,2), ylab="Deviance residual", xlab="Predicted log(HR)")
abline(h = 0,col = "red")

# Weighted Schoenfeld Residuals
temp = cox.zph(fit,transform = "identity")

# The smoothed curve approximates the manner in which the log(HR) changes over time

plot(temp[1],main = "Weighted Schoenfeld for tow duration")
plot(temp[2],main = "Weighted Schoenfeld for handling time")
plot(temp[3],main = "Weighted Schoenfeld for length of fish")
plot(temp[4],main = "Weighted Schoenfeld for log(catch)")
# abline(h = 0,col = "red")

pdf("wgt_sch_res.pdf",height = 6,width = 6)
par(mfrow=c(2,2))
plot(temp[1],main = "Tow duration")
plot(temp[2],main = "Handling time")
plot(temp[3],main = "Length of fish")
plot(temp[4],main = "log(catch)")
dev.off()

# Survival time on the log scale
temp = cox.zph(fit,transform = "log")
plot(temp[1],main = "Weighted Schoenfeld for tow duration")
plot(temp[2],main = "Weighted Schoenfeld for handling time")
plot(temp[3],main = "Weighted Schoenfeld for length of fish")
plot(temp[4],main = "Weighted Schoenfeld for log(catch)")

# Exploring interactions
#Find martingale residuals in the empty model; plot against each covariate
resid<-residuals(fitmin, type="martingale")
plot(halibut$depth, resid, ylab="Residual", xlab="Depth")
lines(lowess(halibut$depth,resid,iter=0))

plot(halibut$length, resid, ylab="Residual", xlab="Length of fish")
lines(lowess(halibut$length,resid,iter=0))

plot(halibut$handling, resid, ylab="Residual", xlab="Handling time")
lines(lowess(halibut$handling,resid,iter=0))

plot(halibut$logcatch, resid, ylab="Residual", xlab="Log(catch)")
lines(lowess(halibut$logcatch,resid,iter=0))

# Model selection with interactions
towdepth<-halibut$towdur*halibut$depth
lengthdepth<-halibut$length*halibut$depth
handlingdepth<-halibut$handling*halibut$depth
towlength<-halibut$towdur*halibut$length
towhandling<-halibut$towdur*halibut$handling
lengthhandling<-halibut$length*halibut$handling

fitmaxX = coxph(Surv(survtime,censor) ~ towdur + depth + length 
                + handling + logcatch + towdepth + lengthdepth 
                + handlingdepth + towlength + towhandling + lengthhandling,
                data = halibut)

fit.backwardX = stepAIC(fitmaxX, scope = list(lower = formula(fitmax),
                                           upper = formula(fitmaxX)),
                       direction = "backward", k = 3,trace = 1)
summary(fit.backwardX)
