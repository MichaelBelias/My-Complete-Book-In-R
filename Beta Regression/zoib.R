## package and data
if(!require(betareg)) install.packages(betareg)
if(!require(ggsci)) install.packages(ggsci)
if(!require(lmtest)) install.packages(lmtest)
if(!require(lme4)) install.packages(lme4)
if(!require(strucchange)) install.packages(strucchange)

data("GasolineYield", package = "zoib")
GasolineYield$batch <- as.factor(GasolineYield$batch)


#### betareg
gy <- betareg(yield ~ temp + batch, data = GasolineYield)
summary(gy)$coeff


#### zoib: fixed effect on batch.
d <- GasolineYield
eg1.fixed <-  zoib(model = yield ~ temp + as.factor(batch)| 1, data = GasolineYield, joint = FALSE,
       random = 0, EUID = 1:nrow(d), zero.inflation = FALSE,
       one.inflation = FALSE, n.iter = 1050, n.thin = 5, n.burn = 50,n.chain = 10 )
sample1 <- eg1.fixed$coeff
# check convergence of the MCMC chains
#### zoib: random effect on batch
eg1.random <- zoib(yield ~ temp | 1 | 1, data = GasolineYield, joint = FALSE,
                   random = 1, EUID = GasolineYield$batch, zero.inflation = FALSE,
                   one.inflation = FALSE, n.iter = 10200, n.thin = 50, n.burn = 200)
sample2 <- eg1.random$coeff


par(mfrow=c(3,4))
summary(sample2)
traceplot(sample2)
autocorr.plot(sample2)
gelman.diag(sample2)


traceplot(sample1)
autocorr.plot(sample1)
gelman.diag(sample1)

### posterior inferences from zoib: fixed
summ1 <- summary(sample1); summ1 <- cbind(summ1$stat[, 1], summ1$quant[, c(1, 5)])
### posterior inferences from zoib: random
summ2 <- summary(sample2); summ2 <- cbind(summ2$stat[, 1], summ2$quant[, c(1, 5)])
summ2<- summ2[-4, ]
### inferences from betareg
summ3 <- summary(gy)
summ3 <- cbind(c(summ3$mean[, 1], summ3$precision[, 1]), confint(gy))
summ3[12,] <- log(summ3[12, ]) # log-precision
### plot
names <- rownames(summ3); names[1] <- "intercept"; names[12] <- "log(precision)"
rownames(summ1) <- names
rownames(summ3) <- names
rownames(summ2) <- names[c(1, 2, 12)]
paraplot(summ1, summ2, summ3,
         legtext = c("zoib: fixed", "zoib: random", "betareg"))



data("AlcoholUse", package = "zoib")
AlcoholUse$Grade <- as.factor(AlcoholUse$Grade)
eg3 <- zoib(Percentage ~ Grade * Gender + MedDays|1|Grade * Gender + MedDays|1,
            data = AlcoholUse, random = 1, EUID = AlcoholUse$County,
            zero.inflation = TRUE, one.inflation = FALSE, joint = FALSE,
            n.iter = 5000, n.thin = 20, n.burn = 1000)
sample1 <- eg3$coeff
summary(sample1)



