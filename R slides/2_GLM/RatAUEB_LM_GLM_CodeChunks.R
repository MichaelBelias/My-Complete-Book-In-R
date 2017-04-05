## This script reproduces the analyses on the slides of the module
## "Linear and Generalized Linear Models" of the R Summer School at
## AUEB. The code chunks in this script should be run serially.

## Author: Ioannis Kosmidis
## Date: 23/06/2014


# You should have ran the chunk PACKAGES below before attending the
# course. If not uncomment and paste into R

## ----PACKAGES--------------------------------------------
# update.packages(ask = FALSE, repos = 'http://cran.rstudio.com/')
# RatAUEBpackages <- c('alr3', 'brglm', 'car', 'caret', 'elasticnet', 'fortunes', 'ggplot2', 'gridExtra', 'GGally', 'leaps', 'lmtest', 'mfp', 'relimp')
# install.packages(RatAUEBpackages, repos = 'http://cran.rstudio.com/')


## ----libraries--------------------------------------------
library(MASS)
library(alr3)
library(car)
library(relimp)
library(fortunes)
library(mfp)
library(ggplot2)
library(GGally)
library(gridExtra)
library(leaps)
library(elasticnet)
library(caret)


## ----brains--------------------------------------------------------------
## Get the brain weight data out of the alr3 package
Brains <- brains
head(Brains)


## ----brains-asis-plot------------
WeightsPlot <- qplot(BodyWt, BrainWt, data = Brains, alpha = I(0.5))
print(WeightsPlot)


## ----brains-log-----------------
LogWeightsPlot <- qplot(log(BodyWt), log(BrainWt), data = Brains, alpha = I(0.5))
print(LogWeightsPlot)


## ----brains-log-fit------------------------------------------------------
# Fit the simple linear model to the brain weight data
LogWeightsLM <- lm(log(BrainWt) ~ log(BodyWt), data = Brains)
# Extract the least squares estimates for beta1 and beta2
(LogWeightsCoefs <- coef(LogWeightsLM))


## ----brains-log-fit-plot-----
LogWeightsPlotLine <- LogWeightsPlot + geom_abline(intercept = LogWeightsCoefs[1], slope = LogWeightsCoefs[2])
print(LogWeightsPlotLine)


## ----brains-fit-plot---------
BodyWtRange <- range(Brains$BodyWt)
BodyWtGrid <- seq(BodyWtRange[1], BodyWtRange[2], length.out = 100)
BrainWtGrid <- exp(LogWeightsCoefs[1])*BodyWtGrid^LogWeightsCoefs[2]
WeightsPlot +  geom_path(data= data.frame(BodyWt = BodyWtGrid, BrainWt = BrainWtGrid))


## ----brains-log-fit-plot-offsets----
LogWeightsFitted <- fitted(LogWeightsLM)
LogWeightsPlotLine + geom_segment(aes(x = log(BodyWt), xend = log(BodyWt), y = log(BrainWt), yend = LogWeightsFitted), col = "blue", alpha = 0.5)


## ----brains-log-fit-summary----------------------------------------------
coef(summary(LogWeightsLM))


## ----brains-log-fit-confint----------------------------------------------
confint(LogWeightsLM, level = 0.95)


## ----exercise-brains-sqrt-fit----------------
## # Plot the square root of the average brain weight versus the sqrt of the average of the body weight
## SqrtWeightsPlot <- qplot(sqrt(BodyWt), sqrt(BrainWt), data = Brains, alpha = I(0.5))
## print(SqrtWeightsPlot)
## # Fit the model on the square root of the average weights
## SqrtWeightsLM <- lm(sqrt(BrainWt) ~ sqrt(BodyWt), data = Brains)
## # Get the coefficients
## (SqrtWeightsCoefs <- coef(SqrtWeightsLM))
## # Test for zero slope
## summary(SqrtWeightsLM)
## # Calculate a 99% confidence interval for the slope
## confint(SqrtWeightsLM, level = 0.99)
## ## Add the fitted line on the scatterplot of the square roots
## SqrtWeightsPlotLine <- SqrtWeightsPlot + geom_abline(intercept = SqrtWeightsCoefs[1], slope = SqrtWeightsCoefs[2])
## print(SqrtWeightsPlotLine)


## ----tree-data-----------------------------------------------------------
load("RatAUEB_LM_GLM_data.RData")
str(Trees)


## ----tree-data-plot----
(TreesPairs <- ggpairs(Trees, alpha = I(0.5)))


## ----tree-model-matrix---------------------------------------------------
FormulaTrees <- volume ~ diameter4 + diameter16 + height
XTrees <- model.matrix(FormulaTrees, data = Trees)
head(XTrees)


## ----tree-model-response-------------------------------------------------
(YTrees <- Trees$volume)


## ----trees-fit-model-----------------------------------------------------
drop(solve(t(XTrees)%*%XTrees)%*%(t(XTrees)%*%YTrees))
TreesLM1 <- lm(volume ~ diameter4 + diameter16 + height, data = Trees)
coef(TreesLM1)


## ----tree-fit-summary----------------------------------------------------
(TreesLM1Summary <- summary(TreesLM1))


## ----tree-fit-confint----------------------------------------------------
confint(TreesLM1, level = 0.95)


## ----trees-fit-expected-value--------------------------------------------
newdata <- data.frame(diameter4 = 14, diameter16 = 10, height = 90)
predict(TreesLM1, newdata = newdata, interval = "confidence")


## ----trees-fit-expected-value-1------------------------------------------
(newdata1 <- rbind(newdata, c(13.5, 8.2, 82)))
predict(TreesLM1, newdata = newdata1, interval = "confidence")


## ----trees-fit-predicted-value-------------------------------------------
newdata <- data.frame(diameter4 = 12.35, diameter16 = 11.77, height = 89.2)
predict(TreesLM1, newdata = newdata, interval = "prediction")


## ----predict-vs-confidence1-----
BodyWtGrid <- seq(BodyWtRange[1], BodyWtRange[2], length.out = 100)
CLogBrainWt <- predict(LogWeightsLM, data.frame(BodyWt = BodyWtGrid), interval = "confidence", level = 0.99)
CLogBrainWt <- data.frame(BodyWt = BodyWtGrid, CLogBrainWt)
PLogBrainWt <- predict(LogWeightsLM, data.frame(BodyWt = BodyWtGrid), interval = "prediction", level = 0.99)
PLogBrainWt <- data.frame(BodyWt = BodyWtGrid, PLogBrainWt)
Pregion <- geom_ribbon(data = PLogBrainWt, aes(y = fit, ymin = lwr, ymax = upr, fill = 'prediction'), alpha = 0.3)
Cregion <- geom_ribbon(data = CLogBrainWt, aes(y = fit, ymin = lwr, ymax = upr, fill = 'confidence'), alpha = 0.3)
LogWeightsPlotLine + Pregion + Cregion


## ----trees-compare-------------------------------------------------------
# Remove the diameters from the full model and refit.
TreesLM2 <- update(TreesLM1, ~ . - diameter4 - diameter16)
# Same as TreesLM2 <- lm(volume ~ height, data = Trees). Now,
# use the anova function to compare full vs nested
anova(TreesLM2, TreesLM1)


## ----brains-fit-ttestftest-----------------------------------------------
#   E.g. Test whether the slope is zero in LogWeightsLM
coef(summary(LogWeightsLM))["log(BodyWt)",]
anova(update(LogWeightsLM, ~ . - log(BodyWt)), LogWeightsLM)


## ----tree-fit-summary-2--------------------------------------------------
TreesLM1Summary
# Also try anova(TreesLM1, update(TreesLM1, . ~ 1))


## ----exercise-trees--------------------------
## TreesLM1
## ## Let's remove diameter4 given that it is not significant (t-test)
## TreesLM3 <- lm(volume ~ diameter16 + height, data = Trees)
## summary(TreesLM3)
## # diameter16 and height seem to explain as much of the variability in volume as the model with both diameters does (the decrease in R^2 is only marginal).
## # Furthermore each of diameter16 and height is highly significant in explaining volume given that the other is in the model.
## # An F-test for comparing TreesLM1 and TreesLM3 is the same as the t-test we used earlier (only one parameter removed)
## anova(TreesLM1, TreesLM3)


## ----duncan-data---------------------------------------------------------
str(Duncan)
head(Duncan)


## ----duncan-data-pairs----
DuncanPairs <- ggpairs(Duncan, color = "type", alpha = I(0.5))


## ----duncan-dummy--------------------------------------------------------
## The model matrix is
head(model.matrix(prestige ~ income + education + type, data = Duncan), n = 12)


## ----duncan-fit----------------------------------------------------------
DuncanLM <- lm(prestige ~ income + education + type, data = Duncan)
summary(DuncanLM)


## ----exercise-dummys-------------------------
## # Fit the model
## DuncanLM2 <- lm(prestige ~ type + education, data = Duncan)
## # Extract the coefficients
## DuncanLM2Coefs <- coef(DuncanLM2)
## ## Extract the implied lines
## lines <- data.frame(
##     intercepts  = DuncanLM2Coefs[1] + c(0, DuncanLM2Coefs[2], DuncanLM2Coefs[3]),
## slopes = c(DuncanLM2Coefs[4], DuncanLM2Coefs[4], DuncanLM2Coefs[4]),
## types <- c('bc', 'prof', 'wc'))
## qplot(education, prestige, data = Duncan, color = type, alpha = I(0.5)) +
##     geom_abline(data = lines, aes(intercept = intercepts, slope = slopes, colour = types))


## ----Duncan-drop1--------------------------------------------------------
drop1(DuncanLM, test = "F")


## ----exercise-line1----------------------------------------
## drop1(DuncanLM, test = 'F')


## ----exercise-line2----------------------------------------
## summary(DuncanLM)


## ----exercise-line3----------------------------------------
## anova(update(DuncanLM, ~ . - type), DuncanLM)


## ----setoptions1----------------------------
cwidth <- getOption("width")
options(width = cwidth + 20)


## ----duncan-interaction-model-----------------------------
head(model.matrix(prestige ~ (income + education)*type, data = Duncan), n = 12)


## ----duncan-interaction-----------------------------------
DuncanLMinteraction <- lm(prestige ~ (income + education)*type, data = Duncan)
summary(DuncanLMinteraction)


## ----setoptions2----------------------------
options(width = cwidth)


## ----duncan-interaction-drop1--------------------------------------------
drop1(DuncanLMinteraction, test = 'F')


## ----exercise-interactions-------------------
## # Fit the model
## DuncanLM3 <- lm(prestige ~ type*education, data = Duncan)
## # Extract the coefficients
## DuncanLM3Coefs <- coef(DuncanLM3)
## ## Extract the implied lines
## lines <- data.frame(
##     intercepts  = DuncanLM3Coefs[1] + c(0, DuncanLM3Coefs[2], DuncanLM3Coefs[3]),
## slopes = DuncanLM3Coefs[4] + c(0, DuncanLM3Coefs[5], DuncanLM3Coefs[6]),
## types <- c('bc', 'prof', 'wc'))
## qplot(education, prestige, data = Duncan, color = type, alpha = I(0.5)) +
##     geom_abline(data = lines, aes(intercept = intercepts, slope = slopes, colour = types))


## ----trees-relimp--------------------------------------------------------
 # Perhaps also check the reference in ?relimp
relimp(TreesLM1, set1 = c(2, 3), set2 = c(4)) # 1 is the Intercept!!


## ----Trees-sequential----------------------------------------------------
anova(TreesLM1)


## ----Trees-sequential-2--------------------------------------------------
  anova(lm(volume ~ diameter16 + diameter4 + height, data = Trees))


## ----exercise-anova--------------------------
## anova(lm(volume ~ diameter4 + diameter16 + height, data = Trees))
## anova(lm(volume ~ diameter16 + diameter4 + height, data = Trees))
## anova(lm(volume ~ diameter4 + height + diameter16, data = Trees))
## anova(lm(volume ~ height + diameter4 + diameter16, data = Trees))
## anova(lm(volume ~ height + diameter16 + diameter4, data = Trees))


## ----fortune-TypeIvsTypeIII----------------------------------------------
fortune("have been fed this nonsense")


## ----trees-response-vs-covariate----
DuncanPairs


## ----wool-data-----------------------------
str(Wool)


## ----wool-residuals1-----------------------
WoolLM <- lm(cycles ~ ., data=Wool)
residualPlots(WoolLM)


## ----wool-residuals2-----------------------
residualPlots(WoolLM, ~ 1)


## ----trees-residual-omitted----
# Recall that TreesLM2 has no diameter covariates
zeroLine <- geom_abline(intercept = 0, slope = 0, lty = 3)
smoothLoess <- stat_smooth(method = "loess", se = FALSE)
q1 <- qplot(Trees$height, resid(TreesLM2)) + zeroLine + smoothLoess
q2 <- qplot(Trees$diameter4, resid(TreesLM2)) + zeroLine  + smoothLoess
q3 <- qplot(Trees$diameter16, resid(TreesLM2)) + zeroLine  + smoothLoess
grid.arrange(q1, q2, q3, ncol = 3)


## ----wool-residuals4----
boxCox(WoolLM)
summary(powerTransform(WoolLM))


## ----wool-residuals5----
LogWoolLM <- update(WoolLM, log(cycles) ~ .)
residualPlots(LogWoolLM)


## ----wool-residuals6----
# Using the plot.lm method of R to construct a Q-Q plot of statndardized residuals
plot(LogWoolLM, which = 2)

## ----wool-residuals6--------------------------------------
## # This is the same as
## qqnorm(rstandard(LogWoolLM), ylim = c(-2.5, 2))
## qqline(rstandard(LogWoolLM), col = "red")
## # We can also use qqPlot method of the car package
## qqPlot(LogWoolLM)
## # qqPLot compares studentized residuals against a t


## ----exercise-Ornstein-1-------------------------------------------------
OrnsteinLM <- lm(interlocks + 1 ~ assets + sector + nation, data = Ornstein)


## ----exercise-Ornstein-2---------------------
## # Checking OrnsteinLM -- looks ugly
## residualPlots(OrnsteinLM)
## plot(OrnsteinLM, which = 2)
## # Pairs plot
## ggpairs(Ornstein)
## # Looks like the relationship between interlocks and assets is closer to logarithmic than linear
## OrnsteinLMLog <- lm(interlocks + 1 ~ log(assets) + sector + nation, data = Ornstein)
## # Make residual plots
## residualPlots(OrnsteinLMLog) ## homoscedasticity and linearity still in doubt...
## plot(OrnsteinLMLog, which = 2) ## Very heavy tails...
## # Try a power transformation of the response
## pT <- powerTransform(OrnsteinLMLog)
## summary(pT) # lambda is between 0.12 and 0.32
## # Get optimal lambda and calculate transformed response
## LambdaOpt <- coef(pT)
## OrnsteinTrans <- transform(Ornstein,
##                        ynew=bcPower(interlocks + 1, LambdaOpt))
## # Fit model on the transformed response
## OrnsteinLMLogTrans <- lm(ynew ~ log(assets) + sector + nation, data = OrnsteinTrans)
## # Model checking -- much better
## residualPlots(OrnsteinLMLogTrans)
## plot(OrnsteinLMLogTrans, which = 2)


## ----hartnagel-serial-1--------------------------------------------------
N <- nrow(Hartnagel)
# Fit the linear model in ?durbinWatsonTest
HartnagelLM <- lm(fconvict ~ tfr + partic + degrees + mconvict, data=Hartnagel)
# Extract the residuals
HResiduals <- residuals(HartnagelLM)
# Observations are ordered here by Hartnage$year so plot the residuals versus order
p1 <- qplot(y = HResiduals) + zeroLine
# Plot the ith versus the (i-1)th residual
p2 <- qplot(HResiduals[-1], HResiduals[-N]) + smoothLoess


## ----hartnagel-serial-2----
grid.arrange(p1, p2, ncol = 2)
durbinWatsonTest(HartnagelLM)


## ----simulate2----
 # simulate some x values and some y values
set.seed(11)
N <- 100; beta1 <- 0.5; beta2 <- -0.2; sigma <- 0.5
simuDat <- data.frame(x = rnorm(N, 2, 0.5))
simuDat <- within(simuDat, {
    y <- beta1 + beta2*x + rnorm(N, 0, sigma)})
trueLine <- geom_abline(intercept = beta1, slope = beta2, lty = 2, lwd = 1.1)
xlims <- xlim(0, 4)
ylims <- ylim(-2.5, 2)
# plot the simulated data
s1 <- qplot(x, y, data = simuDat, alpha = I(0.4)) + xlims + ylims + ggtitle("Data")
# Add two outliers in simuDat
outliers <- data.frame(x = c(1.5, 1.25), y = c(-2, -2.25))
s2 <- s1 + geom_point(data = outliers, aes(x, y), color = "red") + ggtitle("Data and outliers")
# Add a leverage point in simuDat
simuDatL <- rbind(simuDat)
leverages <- data.frame(x = 4,  y = 0.5)
s3 <- s1 + geom_point(data = leverages, aes(x, y), color = "red") + ggtitle("Data and leverage")
# Add a leverage + outlier
extremes <- data.frame(x = 0,  y = -2.25)
s4 <- s1 + geom_point(data = extremes, aes(x, y), color = "red") + ggtitle("Data and (leverage + outlier)")
s1lm <- s1 + stat_smooth(method = "lm", fullrange = TRUE, lwd = 1.1) + trueLine
s2lm <- s2 + stat_smooth(data = rbind(simuDat, outliers), method = "lm", fullrange = TRUE, lwd = 1.1) + trueLine
s3lm <- s3 + stat_smooth(data = rbind(simuDat, leverages), method = "lm", fullrange = TRUE, lwd = 1.1) + trueLine
s4lm <- s4 + stat_smooth(data = rbind(simuDat, extremes), method = "lm", fullrange = TRUE, lwd = 1.1) + trueLine
grid.arrange(s1lm, s2lm, s3lm, s4lm, nrow = 2, ncol = 2)


## ----influence1----
## Let's check influence manually here
hthres <- 2*length(coef(DuncanLM))/nobs(DuncanLM)
r1 <- qplot(y = hatvalues(DuncanLM)) + geom_abline(intercept = hthres, slope = 0, lty = 2) + labs(x = "Observation", y = "Leverage", title = "Leverages")
rthres <- 2
r2 <- qplot(y = abs(rstandard(DuncanLM))) + geom_abline(intercept = rthres, slope = 0, lty = 2) + labs(x = "Observation", y = "|Stand. residual|", title = "Absolute standardized residuals")
grid.arrange(r1, r2, ncol = 2)


## ----influence2----------------------------
# Cook's distances for DuncanLM
cthres <- 8/(nobs(DuncanLM) - 2*length(coef(DuncanLM)))
qplot(y = cooks.distance(DuncanLM)) + geom_abline(intercept = cthres, slope = 0, lty = 2) + labs(x = "Observation", y = "Cook's distance", title = "Cook's distances")


## ----diagnostics----
# plot.lm produces the following 4 plots for model checking by default
par(mfrow = c(2, 2))
plot(DuncanLM)


## ----perfect-collinearity------------------------------------------------
# Add an extra covariate in DuncanLM that is the linear combination of income and education. The addition of that covariate causes perfect collinearity. lm detects this, drops the covariate and proceeds with the rest.
update(DuncanLM, . ~ . + I(2*income + 3*education))


## ----tree-collinearity1----
TreesPairs


## ----tree-collinearity2--------------------------------------------------
# The variance inflation factors for diameter4 and diameter6 are larger than 7 which means that the estimated standard errors for diameter4 and diameter15 are inflated by more than sqrt(7) compared to the case of linear independence
vif(TreesLM1)


## ----tree-collinearity3--------------------------------------------------
# Let's remove diameter4 and recalculate VIF's
TreesLM3 <- lm(volume ~ diameter16 + height, data = Trees)
# Everything looks reasonable in terms of VIF's
vif(TreesLM3)


## ----tree-collinearity4--------------------------------------------------
summary(TreesLM3)


## ----tree-collinearity0--------------------------------------------------
TreesLM1Summary


## ----Chile----------------
# A national survey conducted in April/May 1988 by FLASCO/Chile on the intention to vote for Pinochet. Status quo is a scale of support for the status quo.
ChilePlot <- ggplot(data = Chile, aes(x = statusquo, y = ifelse(vote == "Y", 1, 0)))
AddJitter <- geom_jitter(position = position_jitter(height = 0.05), alpha = I(0.1)) # we add jitter to avoid overplotting
(ChilePlotLM <- ChilePlot + AddJitter + stat_smooth(method = "lm", fullrange = TRUE) + xlim(-2, 2.5))


## ----logitlink---------------
logits <- seq(-10, 10, length.out = 100)
logitlink <- data.frame(etas = logits, probs = plogis(logits))
tmp <- data.frame(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1)
ggplot(logitlink, aes(x = logits, y = probs))  + geom_line(color = "blue") + labs(x = expression(paste(beta^T,x)), y = "probability")


## ----Chile1---------------
# Add the fitted logistic curve to ChilePlot
ChilePlotLM + stat_smooth(method = "glm", family = "binomial", col ="red", fullrange = TRUE)


## ----setoptions3----------------------------
options(width = 45)


## ----family-objects1--------------------------------------
# Binomial with canonical link (logit)
binomial()
# Now Binomial with cloglog link
binomial("cloglog")
# Normal with canonical link (identity)
gaussian()


## ----family-objects2--------------------------------------
# Poisson with canonical link (log)
poisson()
# Gamma with canonical link (inverse)
Gamma()
# The components of a family object
names(poisson("1/mu^2"))


## ----setoptions4----------------------------
options(width = cwidth, cache = FALSE)


## ----Chile-ML------------------------------------------------------------
Chile1 <- na.omit(Chile)
ChileGLM <- glm(ifelse(vote == "Y", 1, 0) ~ statusquo + age + sex + education + income, data = Chile1, family =binomial("logit"))
coef(ChileGLM)


## ----Chile-summary----------------------------------------
# Individual hypotheses tests using the above large sample results. Interpretation is the same as for linear models.
summary(ChileGLM)


## ----Chile-confint-----------------------
# Profile likelihood confidence intervals are used by default for GLM objects
confint(ChileGLM)
# For the familiar "estimate +- quantile * se" confidence intervals use
confint.default(ChileGLM)


## ----Chile-anode------------------------------------------
# R provides a wide range of tests
anova(ChileGLM, test = "LRT")


## ----Chile-Residuals-----------------------------------------------------
head(residuals(ChileGLM, type = "response"))
head(residuals(ChileGLM, type = "pearson"))
head(residuals(ChileGLM, type = "deviance"))


## ----Chile-Residuals1----------------------
# Residual plots for a GLM
residualPlots(ChileGLM)


## ----Chile-Residuals2----------------------
# A Normal Q-Q plot of the standardized deviance residuals
plot(ChileGLM, which = 2)


## ----Chile-influence-----------------------
plot(ChileGLM, which = 5)


## ----bodyfat-leaps-------------------------
data(bodyfat)
Bodyfat <- bodyfat[-match(c("case", "brozek", "density"), names(bodyfat))]
yBodyfat <- c(Bodyfat[, 1])
xBodyfat <- as.matrix(Bodyfat[, -1])
# nvmax is the maximum number of covariates to be considered
BodyfatBS <- regsubsets(x = xBodyfat, y = yBodyfat, nvmax = 13)
plot(BodyfatBS)


## ----bodyfat-stepwise--------------------------------------
## BodyfatLM <- lm(siri ~ ., data = Bodyfat)
## step(BodyfatLM)


## ----bodyfat-train-------------------------
# Use the caret package that can do both LASSO and elastic net.
# Chose 10-fold cross-validation with 10 repetitions
fitControl <- trainControl(method = "repeatedcv", repeats = 10, number = 10,
verboseIter = FALSE)
# We let the fraction (J(beta)/J(hat{beta})) to vary from almost 0 to 1. This is equivalent to setting a grid for lambda1
tuneParsLASSO <- expand.grid(fraction = seq(0.0001, 1, length = 5))
# For elastic net we set lambda (this is lambda2 in previous slides) to vary from 0 to 10
tuneParsEnet <- expand.grid(fraction = seq(0.0001, 1, length = 5),
                            lambda = seq(0, 10, length = 5))


## ----bodyfat-lasso1---------------------------------------
# Apply LASSO to the bodyfat setting
BodyfatLASSO <- train(x = xBodyfat, y = yBodyfat, method = "lasso",
                      trControl = fitControl, tuneGrid = tuneParsLASSO)
BodyfatLASSO


## ----bodyfat-lasso2------------------------
# Plot the root mean squared error as a function of the fraction
BestTuning <- BodyfatLASSO$bestTune$fraction
ggplot(BodyfatLASSO) + geom_vline(xintercept = BestTuning, col = "blue")


## ----bodyfat-lasso3------------------------
# Solution path
plot(BodyfatLASSO$finalModel)
abline(v = BestTuning, lty = 2)


## ----bodyfat-lasso4------------------------
# enet output (from elasticnet package)
BodyfatLASSO$finalModel


## ----bodyfat-lasso5------------------------
# coefficients at optimal fraction (in case you want to see those)
predict(BodyfatLASSO$finalModel, mode = "fraction", s = BestTuning, type = "coefficient")$coef


## ----bodyfat-enet----
# Apply Elastic net to the bodyfat setting: this may take some time
BodyfatEnet <- train(x = xBodyfat, y = yBodyfat, method = "enet",
                     trControl = fitControl, tuneGrid = tuneParsEnet)
# Plot the root mean squared error as a function of the tuning parameters. Here we replace the awkward label for lambda2 with lambda2
ggplot(BodyfatEnet) + labs(x = expression(lambda[2]))
# Again the lasso fit is selected (lambda2 = zero)
