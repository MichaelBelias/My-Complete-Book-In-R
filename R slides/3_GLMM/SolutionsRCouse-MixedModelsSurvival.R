#########################################################################################
# Description: Solutions to the practical exercises for the Mixed Models and Survival   #
#              Analysis part of the 1st R Summer School @ AUEB                          #
# Author: Dimitris Rizopoulos                                                           #
# Last update: 2014-06-11                                                               #
#########################################################################################


###############
# Practical 1 #
###############

# load the workspace the make the PBC dataset available
load(file = ...)

library("lattice")

# average longitudinal evolution per treatment group
xyplot(log(serBilir) ~ year, groups = drug, data = pbc2,    
    type = c("p", "smooth"), lwd = 2)

# average longitudinal evolution per gender
xyplot(log(serBilir) ~ year, groups = sex, data = pbc2,    
    type = c("p", "smooth"), lwd = 2)

# individual longitudinal trajectories (all in one plot)
xyplot(log(serBilir) ~ year, groups = id, data = pbc2, type = "l", lwd = 2)

# individual longitudinal trajectories separately for each subject
pdf() # needed if you do not use Rstudio
xyplot(log(serBilir) ~ year | id, data = pbc2, type = "l", 
    lwd = 2, layout = c(6, 6))
dev.off()

# AUC t-test
slist <- split(pbc2[c("serBilir", "year", "drug")], pbc2$id)
AUC <- do.call(rbind, lapply(slist, function (d) {
    d$logSerBilir <- log(d$serBilir)
    nd <- nrow(d)
    d$auc <- if (nd > 1) 
        sum(with(d, diff(year) * (logSerBilir[-1] + logSerBilir[-nd])/2))
    else NA
    d[1, c("auc", "drug")]
}))

t.test(auc ~ drug, data = AUC)

# variance function over time
library("splines")
res <- resid(lm(log(serBilir) ~ ns(year, 3) * (drug + sex) + age, data = pbc2))
res2 <- res * res
xyplot(res2 ~ year, data = pbc2, type = c("p", "smooth"), lwd = 2)

library("nlme")
# random intercepts
lme1 <- lme(log(serBilir) ~ ns(year, 3) * (drug + sex) + age, data = pbc2, 
    random = ~ 1 | id)

# detailed output, with parameter estimates and stadnard errors
summary(lme1)

# random intercepts & random slopes
lme2 <- lme(log(serBilir) ~ ns(year, 3) * (drug + sex) + age, data = pbc2, 
    random = ~ year | id)

# Likelihood ratio test
anova(lme1, lme2)

# random intercepts & spline random effects
lme3 <- lme(log(serBilir) ~ ns(year, 3) * (drug + sex) + age, data = pbc2, 
    random = ~ ns(year, 3) | id)

# Likelihood ratio test
anova(lme2, lme3)

# drop interaction terms
lme3 <- lme(log(serBilir) ~ ns(year, 3) * (drug + sex) + age, data = pbc2, 
    random = list(id = pdDiag(form = ~ ns(year, 3))))
lme4 <- lme(log(serBilir) ~ ns(year, 3) + drug + sex + age, data = pbc2, 
    random = list(id = pdDiag(form = ~ ns(year, 3))))

anova(lme4, lme3)

# refit models under ML
lme3 <- lme(log(serBilir) ~ ns(year, 3) * (drug + sex) + age, data = pbc2, 
    random = list(id = pdDiag(form = ~ ns(year, 3))), method = "ML")
lme4 <- lme(log(serBilir) ~ ns(year, 3) + drug + sex + age, data = pbc2, 
    random = list(id = pdDiag(form = ~ ns(year, 3))), method = "ML")

anova(lme4, lme3)

# drop age
lme5 <- lme(log(serBilir) ~ ns(year, 3) * (drug + sex), data = pbc2, 
    random = list(id = pdDiag(form = ~ ns(year, 3))), method = "ML")

anova(lme5, lme3)


###############
# Practical 2 #
###############

# load the workspace the make the PBC dataset available
load(file = ...)

library("survival")

# Kaplan-Meier estimator
KM <- survfit(Surv(Time, death) ~ 1, data = aids.id)
KM
plot(KM)

# Kaplan-Meier estimator per treatment group
KM.drug <- survfit(Surv(Time, death) ~ drug, data = aids.id)
KM.drug
plot(KM.drug)

# Kaplan-Meier estimator per gender
KM.gender <- survfit(Surv(Time, death) ~ gender, data = aids.id)
KM.gender
plot(KM.gender)

# Log-rank tests for treatment and gender
survdiff(Surv(Time, death) ~ drug, data = aids.id)
survdiff(Surv(Time, death) ~ gender, data = aids.id)

# Cox model with splines & interactions
cph1 <- coxph(Surv(Time, death) ~ ns(CD4, 4) * (drug + AZT + gender), data = aids.id)
summary(cph1)

# Cox model with main effects only
cph2 <- coxph(Surv(Time, death) ~ ns(CD4, 4) + drug + AZT + gender, data = aids.id)

# Likelihood ratio test
anova(cph2, cph1)

# detailed output, with HRs and p-values
summary(cph2)

# Schoenfeld residuals
schRes <- cox.zph(cph2)
par(mfrow = c(3, 3))
plot(schRes)

