## package and data
if(!require(betareg)) install.packages(betareg)
if(!require(ggsci)) install.packages(ggsci)
if(!require(lmtest)) install.packages(lmtest)
if(!require(lme4)) install.packages(lme4)
if(!require(strucchange)) install.packages(strucchange)

data("ReadingSkills", package = "betareg")


## ReadingSkills OLS and beta regression
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
  data = ReadingSkills, hessian = TRUE)
cl1 <- hcl(c(260, 0), 90, 40)
cl2 <- hcl(c(260, 0), 10, 95)
plot(accuracy ~ iq, data = ReadingSkills, col = cl2[as.numeric(dyslexia)],
  main = "Reading skills data", xlab = "IQ score", ylab = "Reading accuracy",
  pch = c(19, 17)[as.numeric(dyslexia)], cex = 1.5)
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5,
  pch = (1:2)[as.numeric(dyslexia)], col = cl1[as.numeric(dyslexia)])
nd <- data.frame(dyslexia = "no", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[1], lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = cl1[1], lty = 2, lwd = 2)
nd <- data.frame(dyslexia = "yes", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[2], lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = cl1[2], lty = 2, lwd = 2)
legend("topleft", c("control", "dyslexic", "betareg", "lm"),
  lty = c(NA, NA, 1:2), pch = c(19, 17, NA, NA), lwd = 2,
  col = c(cl2, 1, 1), bty = "n")
legend("topleft", c("control", "dyslexic", "betareg", "lm"),
  lty = c(NA, NA, 1:2), pch = c(1, 2, NA, NA),
  col = c(cl1, NA, NA), bty = "n")


## ReadingSkills bias correction/reduction regressions
rs_f <- accuracy ~ dyslexia * iq | dyslexia * iq
rs_ml <- betareg(rs_f, data = ReadingSkills, type = "ML")
rs_bc <- betareg(rs_f, data = ReadingSkills, type = "BC")
rs_br <- betareg(rs_f, data = ReadingSkills, type = "BR")


## ReadingSkills OLS and beta regression
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
                   data = ReadingSkills, hessian = TRUE)
cl1 <- hcl(c(260, 0), 90, 40)
cl2 <- hcl(c(260, 0), 10, 95)
plot(accuracy ~ iq, data = ReadingSkills, col = cl2[as.numeric(dyslexia)],
     main = "Reading skills data", xlab = "IQ score", ylab = "Reading accuracy",
     pch = c(19, 17)[as.numeric(dyslexia)], cex = 1.5)
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5,
       pch = (1:2)[as.numeric(dyslexia)], col = cl1[as.numeric(dyslexia)])
nd <- data.frame(dyslexia = "no", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[1], lty = 1, lwd = 2)
lines(nd$iq, predict(rs_bc, nd), col = cl1[1], lty = 2, lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = "purple", lty = 1, lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = "white", lty = 3, lwd = 2)
nd <- data.frame(dyslexia = "yes", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[2], lty = 1, lwd = 2)
lines(nd$iq, predict(rs_bc, nd), col = cl1[2], lty = 2, lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = "purple", lty = 1, lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = "white", lty = 4, lwd = 2)
legend("topleft", c("control", "dyslexic", "betareg","Bias corrected","logit"),
       lty = c(NA, NA, 1:3), pch = c(19, 17, NA, NA,NA,NA), lwd = 2,
       col = c(cl2, 1, 1,"purple"), bty = "n")
legend("topleft", c("control", "dyslexic","betareg","Bias corrected","logit"),
       lty = c(NA, NA, 1:3), pch = c(1, 2, NA, NA,NA,NA),
       col = c(cl1, NA, NA,NA,"white"), bty = "n")




## extract information
rs_list <- list(rs_ml, rs_bc, rs_br)
sapply(rs_list, coef)
sapply(rs_list, function(x) sqrt(diag(vcov(x))))
sapply(rs_list, logLik)

## visualization of fitted precisions
pr_phi <- sapply(list(
  "Maximum likelihood" = rs_ml,
  "Bias correction" = rs_bc,
  "Bias reduction" = rs_br), predict, type = "precision")
pairs(log(pr_phi), panel = function(x, y, ...) {
  panel.smooth(x, y, ...)
  abline(0, 1, lty = 2)
})


## additional random covariates without assocation to response
set.seed(1071)
n <- nrow(ReadingSkills)
ReadingSkills$x1 <- as.numeric(ReadingSkills$accuracy>0.9 &ReadingSkills$dyslexia==1 )
ReadingSkills$x2 <- runif(n)
ReadingSkills$x3 <- factor(sample(0:1, n, replace = TRUE))

## beta regression tree
rs_tree <- betatree(accuracy ~ iq | iq, ~ dyslexia ,
  data = ReadingSkills, minsplit = 2 )
plot(rs_tree)
coef(rs_tree)
rs_tree
sctest(rs_tree)

rs_tree <- betatree(accuracy ~ iq | iq | dyslexia  ,
                    data = ReadingSkills, minsize = 10)
plot(rs_tree)



## beta regression mixture model
rs_mix <- betamix(accuracy ~ iq, data = ReadingSkills, k = 3,
  extra_components = extraComponent(type = "uniform",
    coef = 0.99, delta = 0.01), nstart = 100)
rs_mix
summary(rs_mix)
table(clusters(rs_mix), ReadingSkills$dyslexia)

## visualization
par(mfrow = c(1, 2))
ix <- as.numeric(ReadingSkills$dyslexia)
prob <- 2 * (posterior(rs_mix)[cbind(seq_along(ix), clusters(rs_mix))] - 0.5)
col3 <- hcl(c(0, 260, 130), 65, 45, fixup = FALSE)
col1 <- col3[clusters(rs_mix)]
col2 <- hcl(c(0, 260, 130)[clusters(rs_mix)], 65 * abs(prob)^1.5, 95 - 50 * abs(prob)^1.5, fixup = FALSE)
plot(accuracy ~ iq, data = ReadingSkills, col = col2, pch = 19, cex = 1.5,
  xlim = c(-2, 2), main = "Mixture model (dyslexia unobserved)")
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5, pch = 1, col = col1)
iq <- -30:30/10
cf <- rbind(coef(rs_mix, model = "mean", component = 1:2), c(qlogis(0.99), 0))
for(i in 1:3) lines(iq, plogis(cf[i, 1] + cf[i, 2] * iq), lwd = 2, col = col3[i])
ix <- as.numeric(ReadingSkills$dyslexia)
col1 <- hcl(c(260, 0), 90, 40)[ix]
col2 <- hcl(c(260, 0), 10, 95)[ix]
plot(accuracy ~ iq, data = ReadingSkills, col = col2, pch = 19,
  cex = 1.5, xlim = c(-2, 2), main = "Partitioned model (dyslexia observed)")
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5, pch = 1, col = col1)
cf <- coef(rs_tree, model = "mean")
col3 <- hcl(c(260, 0), 90, 40)
for(i in 1:2) lines(iq, plogis(cf[i, 1] + cf[i, 2] * iq), lwd = 2, col = col3[i])


## GasolineYield bias correction/reduction regression
data("GasolineYield", package = "betareg")
gy <- lapply(c("ML", "BC", "BR"), function(x)
  betareg(yield ~ batch + temp, data = GasolineYield, type = x))
sapply(gy, coef, model = "precision")
sapply(gy, coef)
sapply(gy, function(x) sqrt(diag(vcov(x))))
sapply(gy, logLik)

## GasolineYield bias correction/reduction regression with precision log link
gy2 <- lapply(c("ML", "BC", "BR"), function(x)
  betareg(yield ~ batch + temp | 1, data = GasolineYield, type = x))
sapply(gy2, logLik)
sapply(gy2, coef)
sapply(gy2, function(x) sqrt(diag(vcov(x))))

