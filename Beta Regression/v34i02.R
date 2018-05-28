###################################################
## Preliminaries
###################################################


if(!require(betareg)) install.packages(betareg)
if(!require(ggsci)) install.packages(ggsci)
if(!require(lmtest)) install.packages(lmtest)
if(!require(lme4)) install.packages(lme4)
if(!require(strucchange)) install.packages(strucchange)

logitTransform <- function(p) { log(p/(1-p)) }
asinTransform <- function(p) { 2*asin(sqrt(p)) }

p <- seq(0.001, 0.999, length.out = 1000)
pLogit <- logitTransform(p)
plot(pLogit, p, type='l', lwd=2, col='red', las=1, ylab='p', xlab='logit(p)')
pAsin <- asinTransform(p)
plot(pAsin,p,  type='l', lwd=2, col='blue', 
     las=1, ylab='p', xlab='arcsine(p)', 
     main = "Two times arcsine - square root transformation of proportion")


rangeScale <- function(x) { (x-min(x)) / (max(x)-min(x)) }

pAsin.scaled <- rangeScale(pAsin)
pLogit.scaled <- rangeScale(pLogit)

plot( pAsin.scaled, p, las=1, type='l', lwd=2, col='blue', ylab='p', xlab='p transformed', 
      main = "Comparison of arcsine and logit transformation", sub = "Scaled into 0-1")
points(pLogit.scaled, p,  type='l', lwd=2, col='red')
text(0.8, 0.8, 'asin', col='blue')
text(0.5, 0.8, 'logit', col='red')


###################################################
## Section 1: Introduction
###################################################
## Visualization of various beta distributions
par(mfrow = c(1, 2), mar = c(4.1, 4.1, 4.1, 0.1))
dbeta2 <- function(x, mu, phi = 1) dbeta(x, mu * phi, (1 - mu) * phi)
x <- seq(from = 0.01, to = 0.99, length = 200)
xx <- cbind(x, x, x, x, x)

yy <- cbind(
  dbeta2(x, 0.10, 5),
  dbeta2(x, 0.25, 5),
  dbeta2(x, 0.50, 5),
  dbeta2(x, 0.75, 5),
  dbeta2(x, 0.90, 5)
)
matplot(xx, yy, type = "l", xlab = "y", ylab = "Density", main = expression(phi == 5),
  lty = 1, col = "black", ylim = c(0, 15))
text(0.05, 12  , "0.10")
text(0.95, 12  , "0.90")
text(0.22,  2.8, "0.25")
text(0.78,  2.8, "0.75")
text(0.50,  2.3, "0.50")

yy <- cbind(
  dbeta2(x, 0.10, 100),
  dbeta2(x, 0.25, 100),
  dbeta2(x, 0.50, 100),
  dbeta2(x, 0.75, 100),
  dbeta2(x, 0.90, 100)
)
matplot(xx, yy, type = "l", xlab = "y", ylab = "", main = expression(phi == 100),
  lty = 1, col = "black", ylim = c(0, 15))
text(0.10, 14.5, "0.10")
text(0.90, 14.5, "0.90")
text(0.25,  9.8, "0.25")
text(0.75,  9.8, "0.75")
text(0.50,  8.6, "0.50")


###################################################
### Section 4.1: Prater's gasoline yield data
###################################################
## Data
data("GasolineYield", package = "betareg")

## Beta regression with logit link
gy_logit <- betareg(yield ~ batch + temp, data = GasolineYield)
summary(gy_logit)

## Beta regression with log-log link (for visualization)
gy_loglog <- betareg(yield ~ batch + temp, data = GasolineYield,
  link = "loglog")

gy_logit2 <- lm(qlogis(yield) ~ batch + temp, data = GasolineYield)


## Visualization of data and models
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
redblue <- hcl(c(0, 260), 90, 40)
plot(yield ~ temp, data = GasolineYield, type = "n",
  ylab = "Proportion of crude oil converted to gasoline", 
  xlab = "Temperature at which all gasoline has vaporized",
  main = "Prater's gasoline yield data")
points(yield ~ temp, data = GasolineYield, cex = 1.75, 
  pch = 19, col = rev(gray.colors(10))[as.numeric(batch)])
points(yield ~ temp, data = GasolineYield, cex = 1.75)
legend("topleft", as.character(1:10), title = "Batch",
  col = rev(gray.colors(10)), pch = 19, bty = "n")
legend("topleft", as.character(1:10), title = "Batch", pch = 1, bty = "n")
lines(150:500, predict(gy_logit, 
  newdata = data.frame(temp = 150:500, batch = "6")),
  col = redblue[2], lwd = 2, lty = 2)
lines(150:500, predict(gy_loglog, 
  newdata = data.frame(temp = 150:500, batch = "6")),
  col = redblue[1], lwd = 2)
lines(150:500, plogis(predict(gy_logit2, 
                       newdata = data.frame(temp = 150:500, batch = "6"))),col = "green", lwd = 2)
legend("bottomright", c("log-log", "logit", "logit-transformation"),
  col = c("red","blue","green"), lty = 1:2, lwd = 2, bty = "n")

## Omitting influential observation #4
gy_logit4 <- update(gy_logit, subset = -4)
coef(gy_logit, model = "precision")
coef(gy_logit4, model = "precision")

## Diagnostic displays
par(mfrow = c(3, 2))
set.seed(123)
plot(gy_logit, which = 1:4, type = "pearson")
plot(gy_logit, which = 5, type = "deviance", sub.caption = "")
plot(gy_logit, which = 1, type = "deviance", sub.caption = "")


###################################################
### Section 4.1: Household food expenditures
###################################################
## Data
data("FoodExpenditure", package = "betareg")

## OLS regression
fe_lm <- lm(I(food/income) ~ income + persons, data = FoodExpenditure)
par(mfrow=c(2,2))
plot(fe_lm)
## OLS regression

fe_lmer <- lmer(I(food/income) ~ income + persons + (1|persons), data = FoodExpenditure)
summary(fe_lmer)


## Beta regression
fe_beta <- betareg(I(food/income) ~ income + persons,
  data = FoodExpenditure)

## Beta regression with variable dispersion (for visualization)
fe_beta2 <- betareg(I(food/income) ~ income + persons | persons,
  data = FoodExpenditure)
plot(fe_beta2)



## Visualization of data and models
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
redblueblack <- hcl(c(0, 260, 0), c(90, 90, 0), c(40, 40, 0))
plot(I(food/income) ~ income, data = FoodExpenditure,
  xlab = "Household income", ylab = "Proportion of food expenditures",
  main = "Food expenditures data", type = "n", ylim = c(0.04, 0.57))
points(I(food/income) ~ income, data = FoodExpenditure, cex = 1.75, pch = 19,
  col = rev(gray.colors(7))[persons])
points(I(food/income) ~ income, data = FoodExpenditure, cex = 1.75)
legend("bottomleft", rev(as.character(sort(unique(FoodExpenditure$persons)))),
  title = "Persons", col = gray.colors(7), pch = 19, bty = "n")
legend("bottomleft", rev(as.character(sort(unique(FoodExpenditure$persons)))),
  title = "Persons", pch = 1, bty = "n")
lines(10:100, predict(fe_lm, 
  newdata = data.frame(income = 10:100, persons = mean(FoodExpenditure$persons))),
  col = redblueblack[3], lwd = 2, lty = 2)
lines(10:100, predict(fe_beta, 
  newdata = data.frame(income = 10:100, persons = mean(FoodExpenditure$persons))),
  col = redblueblack[2], lwd = 2, lty = 5)
lines(10:100, predict(fe_beta2, 
  newdata = data.frame(income = 10:100, persons = mean(FoodExpenditure$persons))),
  col = redblueblack[1], lwd = 2)
legend("topright", c("logit, var. disp.", "logit, fix. disp.", "lm"),
  col = redblueblack, lty = c(1, 5, 2), lwd = 2, bty = "n")

## Breusch-Pagan test for heteroskedasticity in OLS regression

bptest(fe_lm)


###################################################
### Section 4.2: Variable dispersion model
###################################################
## For Prater's gasoline yield data
gy_logit2 <- betareg(yield ~ batch + temp | temp, data = GasolineYield)
summary(gy_logit2)

## LR test against fixed dispersion model
lrtest(gy_logit, gy_logit2)

## Fixed vs. variable dispersion in food expenditure data
lrtest(fe_beta, fe_beta2)
AIC(fe_beta, fe_beta2, k = log(nrow(FoodExpenditure)))


###################################################
### Section 4.3: Selection of different link functions
###################################################
## Comparisons for Prater's gasoline yield data
summary(gy_logit)$pseudo.r.squared
summary(gy_loglog)$pseudo.r.squared
AIC(gy_logit, gy_logit2, gy_loglog)

## RESET-type test examples
lrtest(gy_logit, . ~ . + I(predict(gy_logit, type = "link")^2))
lrtest(gy_loglog, . ~ . + I(predict(gy_loglog, type = "link")^2))

## Graphical comparison of model fits
plot(abs(residuals(gy_loglog, type = "response")),
  abs(residuals(gy_logit, type = "response")))
abline(0, 1, lty = 2)

## Convergence comparisons
gy_loglog2 <- update(gy_loglog, link.phi = "log")
summary(gy_loglog2)$iterations
summary(gy_loglog)$iterations

## Goodness-of-fit comparisons for food expenditure data
sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"),
  function(x) logLik(update(fe_beta2, link = x)))


###################################################
### Section 5.1: Dyslexia and IQ predicting reading accuracy
###################################################
## Data
data("ReadingSkills", package = "betareg")

## OLS regression
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)
coeftest(rs_ols)

rs_arsine <- lm(asin(sqrt(accuracy)) ~ dyslexia * iq, data = ReadingSkills)
coeftest(rs_arsine)
sin_square <- function(x) sin(x)^2

## Beta regression with variable dispersion
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia * iq ,
                   data = ReadingSkills, hessian = TRUE)
coeftest(rs_beta)


## Beta regression without variable dispersion
rs_beta_fixed <- betareg(accuracy ~ dyslexia * iq ,
                   data = ReadingSkills, hessian = TRUE)
coeftest(rs_beta_fixed)


rs_ols_lmer <- lmer(qlogis(accuracy) ~ dyslexia * iq  +(iq|dyslexia) , data = ReadingSkills)
summary(rs_ols_lmer)
coeftest(rs_beta)




## Visualization of data and models
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
cl1 <- hcl(c(260, 0), 90, 40)
cl2 <- hcl(c(260, 0), 10, 95)
plot(accuracy ~ iq, data = ReadingSkills, col = cl2[as.numeric(dyslexia)],
  main = "Reading skills data", xlab = "IQ score", ylab = "Reading accuracy",
  pch = c(19, 17)[as.numeric(dyslexia)], cex = 1.5 , sub =  "Figure 1")
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5,
  pch = (1:2)[as.numeric(dyslexia)], col = cl1[as.numeric(dyslexia)])
nd1 <- data.frame(dyslexia = "no", iq = -30:30/10)
lines(nd1$iq, predict(rs_beta, nd1), col = cl1[1], lwd = 2)
lines(nd1$iq, plogis(predict(rs_ols, nd1)), col = cl1[1], lty = 2, lwd = 2)
lines(nd1$iq, sin_square(predict(rs_arsine, nd1)), col = cl1[1],lty = 3, lwd = 2)
nd2 <- data.frame(dyslexia = "yes", iq = -30:30/10)
lines(nd2$iq, predict(rs_beta, nd2), col = cl1[2], lwd = 2)
lines(nd2$iq, plogis(predict(rs_ols, nd2)), col = cl1[2], lty = 2, lwd = 2)
lines(nd2$iq, sin_square(predict(rs_arsine, nd2)), col = cl1[2], lty = 3, lwd = 2)
legend("topleft", c("control", "dyslexic", "betareg", "logit","arcsine"),
  lty = c(NA, NA, 1:3), pch = c(19, 17, NA, NA,NA), lwd = 2,
  col = c(cl1[1],cl2[2], "black", "black","black"), bty = "n")
legend("topleft", c("control", "dyslexic","betareg", "logit","arcsine"),
       lty = c(NA, NA, 1:3), pch = c(1, 2, NA, NA,NA),
       col = c(cl1, NA, NA,NA), bty = "n")


ggplot(ReadingSkills, aes(x = iq, y = accuracy)) +
  geom_point(size = 4, aes(fill = dyslexia), shape = 21) +
  scale_fill_lancet() +
  geom_line(data = ReadingSkills[ReadingSkills$dyslexia=="yes",],aes(y = predict(rs_beta,ReadingSkills[ReadingSkills$dyslexia=="yes",]),
                colour = "red", linetype = "beta-reg"),size = 1) +
  geom_line(data = ReadingSkills[ReadingSkills$dyslexia=="yes",],
            aes(y = plogis(predict(rs_ols, ReadingSkills[ReadingSkills$dyslexia=="yes",])), 
                colour = "red", linetype = "OLS"),size = 1)+
  geom_line(data = ReadingSkills[ReadingSkills$dyslexia=="no",],
            aes(y = predict(rs_beta,ReadingSkills[ReadingSkills$dyslexia=="no",]),
            colour = "blue", linetype = "beta-reg"),size = 1) +
  geom_line(data = ReadingSkills[ReadingSkills$dyslexia=="no",],
            aes(y = plogis(predict(rs_ols, ReadingSkills[ReadingSkills$dyslexia=="no",])), 
            colour = "blue", linetype = "OLS"),size = 1)  +
  scale_colour_manual("", values = c("blue", "red")) +
  scale_linetype_manual("", values = c("solid", "dashed")) +
  theme_bw()



###################################################
### Section 5.2:  Structural change testing in beta regressions
###################################################
## Artificial data generation
set.seed(123)
y1 <- c(rbeta(150, 0.3 * 4, 0.7 * 4), rbeta(50, 0.5 * 4, 0.5 * 4))
y2 <- c(rbeta(100, 0.3 * 4, 0.7 * 4), rbeta(100, 0.3 * 8, 0.7 * 8))

## Generalized empirical fluctuation processes

y1_gefp <- gefp(y1 ~ 1, fit = betareg)
y2_gefp <- gefp(y2 ~ 1, fit = betareg)

## Graphical test
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
plot(y1_gefp, aggregate = FALSE)
plot(y2_gefp, aggregate = FALSE)

