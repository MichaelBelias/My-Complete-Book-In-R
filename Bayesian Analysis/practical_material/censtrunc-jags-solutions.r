library(rjags)

################################################################################
### 3.1 CENSORING


### Question 1(a). Rounded chickens implemented using interval censoring.

## Model code.
cchick.mod <- "model {
 for (i in 1:9) {
  y[i] ~ dnorm(mu, 1)
  is.censored[i] ~ dinterval(y[i], ycut)
 }
 mu ~ dunif(0, 20)
}"

## Data:
# "is.censored" gives the index into "ycut" giving lower bound for the corresponding y value
cchick.dat <- list(
    y           = rep(NA, 9),
    is.censored = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
    ycut = c(5.5, 6.5, 7.5, 8.5)
)

## Initial values: each unknown y value initialised to the midpoint of the interval.
cchick.in <- list(mu=7, y=c(6,6,6,7,7,7,8,8,8))

cchick.jag <- jags.model(textConnection(cchick.mod), data=cchick.dat, inits=cchick.in)
update(cchick.jag, 1000) # burn-in
sam <- coda.samples(cchick.jag, var=c("y","mu"), n.iter=100000)
summary(sam)

## The posterior summmary statistics for y[1] to y[3] are the same up to Monte Carlo error, as they have exactly the same prior and likelihood.  Likewise for y[4] to y[6] and y[7] to y[9].

## The posterior mean and median of mu are 7, as expected, since the mean of the data (integrated over the intervals) is exactly 7.   Similarly, the point estimates of y[4],y[5],y[6] are all 7.

## The point estimates of y[1],y[2],y[3] are around 6.1, and y[7],y[8],y[9] are around 7.9.   This is because the y[i] are assumed to be drawn from a normal distribution with a mean of mu, estimated as 7, and SD of exactly 1.  Thus the posterior distribution of each y[i] combines a N(mu,1) prior with the observation that y[i] is censored on an interval centered at the integer 6, 7, or 8.

## If the prior standard deviation is made smaller, then there will be more "shrinkage" of the highest and lowest y[i] towards the prior mean:  the posteriors of y[1],y[2],y[3] will be pulled further above 6, and the posteriors of y[7],y[8],y[9] will be pulled further below 8.

## On the other hand, as the standard deviation increases, the data will dominate the prior, there will be no shrinkage, and the point estimate of each y[i] tend towards exactly the integer 6, 7 or 8.


### Question 1(b).  Estimating the posterior SD

cchick.sd.mod <- "model {
 for (i in 1:9) {
  y[i] ~ dnorm(mu, tau)
  is.censored[i] ~ dinterval(y[i], ycut)
 }
 mu ~ dunif(0, 20)
 tau <- 1 / (sig*sig)
 sig ~ dunif(0, 1000)
}"
cchick.sd.in <- list(mu=7, sig=1, y=c(6,6,6,7,7,7,8,8,8))
cchick.jag <- jags.model(textConnection(cchick.sd.mod), data=cchick.dat, inits=cchick.sd.in)
update(cchick.jag, 1000) # burn-in
sam <- coda.samples(cchick.jag, var=c("y","mu","sig"), n.iter=100000)
summary(sam)

## The SD is estimated to be around 1, though this will be mildly sensitive to the prior given only 9 data points.


### Question 2.

### Model, inits: cchick.mod and cchick.in as defined in part 1 can be used again.
### Only the data changes, since weights over 8 are right-censored on (7.5, )
cchick.rcens.dat <- list(
    y           = rep(NA, 9),
    is.censored = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
    ycut = c(5.5, 6.5, 7.5)
)

cchick.jag <- jags.model(textConnection(cchick.mod), data=cchick.rcens.dat, inits=cchick.in)

sam <- coda.samples(cchick.jag, var=c("y","mu"), n.iter=100000)

summary(sam)

### Now that we assume values over 8.5 are possible for weights recorded as 8, the posterior mean and median for y[7],y[8] and y[9] are above 8, and the credible intervals are wider compared to part 1(a).  Though since the y[i] are assumed to come from a normal with SD 1, values over 9 will not be very common.



################################################################################
### 3.2 TRUNCATION

### Truncated Normal model for the six chickens data using T()

tchick.mod <- "model {
 for (i in 1:6) {
  y[i] ~ dnorm(mu, 1)T(,8)
 }
 mu ~ dunif(0, 10)
}"

tchick.dat <- list(y = c(6, 6, 6, 7, 7, 7))
tchick.in <- list(mu=6)
tchick.jag <- jags.model(textConnection(tchick.mod), data=tchick.dat, inits=tchick.in)
update(tchick.jag, 1000) # burn-in
sam <- coda.samples(tchick.jag, var=c("mu"), n.iter=10000)
summary(sam)

### Note that mu can theoretically be outside the truncation interval,
### so we should allow the prior to extend above 8. Though we might
### expect sampling problems if mu is too far outside these bounds.



################################################################################
### 3.3 ZEROES TRICK

## 1. Truncated normal model for chickens data again

tchick.mod <- "
data {
  for (i in 1:6){
    z[i] <- 0
  }
  U <- 8
  C <- 1000
}
model {
 for (i in 1:6) {
  z[i] ~ dpois(phi[i])
  phi[i] <- -logL[i] + C
  logL[i] <- - 0.5*(y[i] - mu)*(y[i] - mu) -
             log(phi((U - mu)))
 }
 mu ~ dunif(0, 10)
}"
tchick.jag <- jags.model(textConnection(tchick.mod), data=tchick.dat, inits=tchick.in)
update(tchick.jag, 1000) # burn-in
sam <- coda.samples(tchick.jag, var=c("mu"), n.iter=10000)
summary(sam)

## Should match the results from using T(), up to Monte Carlo error.


## 2. Zeroes trick to implement censored data

zcens.mod <- "
data {
  for (i in 1:3){
    z[i] <- 0
  }
}
model {
 for (i in 1:6) {
  y[i] ~ dnorm(mu, 1)
 }
 for (i in 1:3){
  z[i] ~ dpois(phi[i])
  phi[i] <- -logL[i] + 1000
  logL[i] <- log(1 - phi(8 - mu))
 }
 ypred ~ dnorm(mu, 1)T(8,)
 mu ~ dunif(0, 20)
}"

zcens.dat <- list(
    y = c(6.1, 6.3, 6.4, 7.1, 7.2, 7.3)
)
zcens.in <- list(mu=7)

zcens.jag <- jags.model(textConnection(zcens.mod), data=zcens.dat, inits=zcens.in)
update(zcens.jag, 1000) # burn-in
sam <- coda.samples(zcens.jag, var=c("mu","ypred"), n.iter=100000)
summary(sam)
traceplot(sam)

### equivalent model from slides
cens.mod <- "
model {
    for (i in 1:9) {
        y[i] ~ dnorm(mu, 1)
        is.censored[i] ~ dinterval(y[i], 8)
    }
    mu ~ dunif(0, 20)
}"
cens.data <- list(
    y = c(6.1, 6.3, 6.4, 7.1, 7.2, 7.3, NA, NA, NA),
    is.censored = c(0, 0, 0, 0, 0, 0, 1, 1, 1)
)
cens.in <- list(mu=7, y=c(NA,NA,NA, NA, NA, NA, 8.1, 8.1, 8.1))

cens.jag <- jags.model(textConnection(cens.mod), data=cens.data, inits=cens.in)
update(cens.jag, 1000) # burn-in
sam <- coda.samples(cens.jag, var=c("mu","y[7]"), n.iter=1000000)
summary(sam)

## Results should agree



################################################################################
### 3.4 SURVIVAL ANALYSIS

install.packages("JM")  ### For the AIDS data
install.packages("flexsurv")   ## For parametric survival models with left truncation in R
library(JM)
library(survival)  ## For standard parametric survival models in R
library(flexsurv)

### Basic survival model in R for the AIDS data
sr <- survreg(Surv(Time, death) ~ 1, data=aids.id, dist="weibull")
summary(sr)

## Fit the same model as Bayesian, in JAGS.
## JAGS/BUGS parameterisation is different : S(x) = exp(- lambda * x^nu)
## R survreg uses S(x) = exp(- (x/beta) ^ nu)
## where (1/beta)^nu = lambda

aids.data <- list(
    ## set censored survival times to be missing
    y = ifelse(aids.id$death, aids.id$Time, NA),
    ## maximum survival time: point beyond which nobody in the data is censored
    ## defined so that is.censored[i] = 0 for death times, and is.censored[i] = 1 for censoring times
    c = ifelse(aids.id$death, max(aids.id$Time)+1, aids.id$Time),
    is.censored = 1 - aids.id$death,
    n = nrow(aids.id)
)

aids.model <- "
model {
 for (i in 1:n) {
  y[i] ~ dweib(nu, lambda)
  is.censored[i] ~ dinterval(y[i], c[i])
 }

### SILLY PRIORS

## lambda: precision of 0.0001 = SD 100: so 95% prior credible interval for mean = 1/lambda is
## (exp(-200), exp(200)), which gives prior probability of about of 0.5 to implausible values for mean survival time for humans.

# lambda <- exp(loglambda)
# loglambda ~ dnorm(0, 0.0001)

## Similarly for nu, prior SD of 100 for log(nu) translates to implausibly big odds ratios for doubled time.

# nu <- exp(lognu)
# lognu ~ dnorm(0, 0.0001)

### MORE SENSIBLE PRIORS

 meansurv ~ dunif(0, 60*12) # mean survival in months
 lambda <- 1 / meansurv
 loghr ~ dunif(log(1/100), log(100))
 nu <- loghr / log(2) + 1

## for comparison with R survreg()
 beta <- pow(lambda, -1/nu)
 rintercept <- log(beta) # (Intercept) reported by survreg
 rlogscale <- log(1/nu)  # Log(scale) reported by survreg summary() output
}
"

library(rjags)
## Sensible to initialise the survival time y here, or else JAGS could
## auto-generate an initial value which is less than the censoring
## time, giving the "inconsistent with observed parents" error
aids.in <- list(lognu=0, loglambda=0,
                y=ifelse(is.na(aids.data$y),
                         aids.data$c + 1, NA))
surv.jag <- jags.model(textConnection(aids.model), aids.data, aids.in)
sam <- coda.samples(surv.jag, n.iter=3000, variable.names=c("nu", "lambda", "rintercept", "rlogscale"))
###
summary(sam)
### Check it matches R results, but many more JAGS iterations probably
### needed for sure convergence and precision.


### With a time dependent covariate: CD4 count

aids.tdc.data <- list(
    start = aids$start,
    stop = ifelse(aids$event, aids$stop, NA),
    c = ifelse(aids$event, max(aids$stop)+1, aids$stop),
    is.censored = 1 - aids$event,
    CD4 = aids$CD4,
    n = nrow(aids)
)

aids.tdc.model <- "
model {
 for (i in 1:n) {
  stop[i] ~ dweib(nu, lambda[i])T(start[i], )
  is.censored[i] ~ dinterval(stop[i], c[i])
  log(lambda[i]) <- a[1] + a[2]*CD4[i]
 }

 ## Priors
 meansurv ~ dunif(0, 60*12) # mean survival in months
 a[1] <- log(1/meansurv)
 loghr ~ dunif(log(1/100), log(100)) # log HR for doubled time
 nu <- loghr / log(2) + 1
 a[2] ~ dunif(log(1/100), log(100)) # log HR for 1 unit of CD4

 scale <- exp(a[1]) # reported by flexsurvreg with dist=\"weibullPH\"
}
"

aids.tdc.in <- list(meansurv=5, loghr=0.2, a=c(NA, 0),
                stop=ifelse(is.na(aids.tdc.data$stop), aids.tdc.data$c + 1, NA))
surv.tdc.jag <- jags.model(textConnection(aids.tdc.model), aids.tdc.data, aids.tdc.in) # slow again
sam <- coda.samples(surv.tdc.jag, n.iter=300, variable.names=c("nu","a","scale")) # slow, need more for convergence
plot(sam)
summary(sam)

## Equivalent (non-Bayesian) model in flexsurv
fs <- flexsurvreg(Surv(start, stop, event) ~ CD4, data=aids, dist="weibullPH") # can ignore the warnings
fs ## OK matches roughly, but the JAGS run needs many more iterations

