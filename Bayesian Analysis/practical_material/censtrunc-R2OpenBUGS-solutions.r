library(R2OpenBUGS)

################################################################################
### 3.1 CENSORING


### Question 1(a). Rounded chickens implemented using interval censoring.
## Model code in models/cchick_mod_solution.txt

## Data:
cchick.dat <- list(
    low = c(5.5, 5.5, 5.5, 6.5, 6.5, 6.5, 7.5, 7.5, 7.5), 
    upp = c(6.5, 6.5, 6.5, 7.5, 7.5, 7.5, 8.5, 8.5, 8.5)
)

## Initial values: each unknown y value initialised to the midpoint of the interval.
cchick.in <- list(list(mu=7, y=c(6,6,6,7,7,7,8,8,8)))

cchick.sim <- bugs(data=cchick.dat,
                   inits=cchick.in,
                   model="cchick_mod_solution.txt",
                   working.directory="models",
                   n.chains=1, n.burnin=1000, n.iter=10000,
                   parameters.to.save=c("y","mu"), DIC=FALSE)
print(cchick.sim, digits=3)

## Current: 1 chains, each with 10000 iterations (first 1000 discarded)
## Cumulative: n.sims = 9000 iterations saved
##       mean    sd  2.5%   25%   50%   75% 97.5%
## y[1] 6.077 0.276 5.545 5.859 6.117 6.312 6.481
## y[2] 6.083 0.276 5.547 5.865 6.121 6.320 6.483
## y[3] 6.078 0.280 5.544 5.854 6.112 6.323 6.484
## y[4] 6.996 0.284 6.527 6.751 6.991 7.241 7.471
## y[5] 7.000 0.284 6.530 6.757 7.000 7.244 7.471
## y[6] 7.000 0.283 6.528 6.758 7.001 7.238 7.475
## y[7] 7.922 0.276 7.518 7.685 7.887 8.140 8.456
## y[8] 7.924 0.274 7.519 7.686 7.892 8.139 8.454
## y[9] 7.918 0.277 7.517 7.677 7.878 8.136 8.455
## mu   6.999 0.349 6.319 6.762 6.999 7.234 7.691

## The posterior summmary statistics for y[1] to y[3] are the same up to Monte Carlo error, as they have exactly the same prior and likelihood.  Likewise for y[4] to y[6] and y[7] to y[9].

## The posterior mean and median of mu are 7, as expected, since the mean of the data (integrated over the intervals) is exactly 7.   Similarly, the point estimates of y[4],y[5],y[6] are all 7.

## The point estimates of y[1],y[2],y[3] are around 6.1, and y[7],y[8],y[9] are around 7.9.   This is because the y[i] are assumed to be drawn from a normal distribution with a mean of mu, estimated as 7, and SD of exactly 1.  Thus the posterior distribution of each y[i] combines a N(mu,1) prior with the observation that y[i] is censored on an interval centered at the integer 6, 7, or 8.

## If the prior standard deviation is made smaller, then there will be more "shrinkage" of the highest and lowest y[i] towards the prior mean:  the posteriors of y[1],y[2],y[3] will be pulled further above 6, and the posteriors of y[7],y[8],y[9] will be pulled further below 8.

## On the other hand, as the standard deviation increases, the data will dominate the prior, there will be no shrinkage, and the point estimate of each y[i] tends towards exactly the integer 6, 7 or 8.


### Question 1(b).  Estimating the posterior SD
## Model code in models/cchick_sd_mod_solution.txt

cchick.sd.in <- list(mu=7, sig=1, y=c(6,6,6,7,7,7,8,8,8))
cchick.sd.sim <- bugs(data=cchick.dat, inits=cchick.in, model="cchick_sd_mod_solution.txt",
                   working.directory="models",
                   n.burnin=1000, n.iter=100000, n.chains=1,
                   param=c("y","mu","sig"), DIC=FALSE)
print(cchick.sd.sim, digits=3)

##            mean    sd   2.5%    25%    50%    75%  97.5%
## y[1]      6.101 0.274  5.551  5.892  6.144  6.338  6.485
## y[2]      6.101 0.275  5.551  5.889  6.144  6.339  6.485
## y[3]      6.102 0.276  5.550  5.891  6.147  6.341  6.485
## y[4]      7.000 0.282  6.528  6.760  7.002  7.238  7.472
## y[5]      7.001 0.283  6.527  6.761  7.001  7.240  7.473
## y[6]      7.000 0.283  6.529  6.759  7.000  7.240  7.472
## y[7]      7.898 0.275  7.515  7.660  7.856  8.108  8.450
## y[8]      7.899 0.276  7.516  7.661  7.854  8.112  8.451
## y[9]      7.898 0.275  7.515  7.662  7.855  8.108  8.450
## mu        7.001 0.361  6.277  6.786  7.001  7.218  7.720
## sig       0.985 0.346  0.530  0.752  0.917  1.138  1.844

## The SD is estimated to be around 1, though this will be mildly sensitive to the prior given only 9 data points.


### Question 2.

### Weights over 8 are right-censored on (7.5, )
### Model code in models/cchick_rcens_mod_solution.txt

cchick.rcens.dat <- list(
    low = c(5.5, 5.5, 5.5, 6.5, 6.5, 6.5, 7.5, 7.5, 7.5), 
    upp = c(6.5, 6.5, 6.5, 7.5, 7.5, 7.5)
)

cchick.rcens.sim <- bugs(data=cchick.rcens.dat, inits=cchick.in,
                         model="cchick_rcens_mod_solution.txt",
                         working.directory="models",
                         n.chains=1, n.burnin=1000, n.iter=10000,
                         parameters.to.save=c("y","mu"),
                         DIC=FALSE)
print(cchick.rcens.sim, digits=3)

##       mean    sd  2.5%   25%   50%   75% 97.5%
## y[1] 6.086 0.278 5.542 5.868 6.124 6.327 6.481
## y[2] 6.087 0.277 5.546 5.869 6.124 6.327 6.485
## y[3] 6.083 0.277 5.550 5.864 6.120 6.324 6.483
## y[4] 7.010 0.286 6.527 6.766 7.014 7.257 7.474
## y[5] 7.011 0.284 6.528 6.771 7.017 7.256 7.475
## y[6] 7.003 0.283 6.529 6.764 7.006 7.242 7.475
## y[7] 8.179 0.555 7.525 7.748 8.045 8.467 9.557
## y[8] 8.176 0.549 7.523 7.747 8.040 8.472 9.544
## y[9] 8.182 0.551 7.523 7.749 8.050 8.489 9.525
## mu   7.090 0.362 6.382 6.847 7.086 7.331 7.808

### Now that we assume values over 8.5 are possible for weights recorded as 8, the posterior mean and median for y[7],y[8] and y[9] are above 8, and the credible intervals are wider compared to part 1(a).  Though since the y[i] are assumed to come from a normal with SD 1, values over 9 will not be very common.



################################################################################
### 3.2 TRUNCATION

### Truncated Normal model for the six chickens data using T() 
## Model code in models/tchick_mod_solution.txt

tchick.dat <- list(y = c(6, 6, 6, 7, 7, 7))
tchick.in <- list(mu=6)
tchick.in <- list(list(mu=7, sig=1, y=c(6,6,6,7,7,7,8,8,8)))
tchick.sim <- bugs(data=tchick.dat, inits=tchick.in, model="tchick_mod_solution.txt",
                   working.directory="models",
                   n.burnin=1000, n.iter=10000, n.chains=1,
                   param=c("mu"), DIC=FALSE)
print(tchick.sim, digits=2)

## Current: 1 chains, each with 10000 iterations (first 1000 discarded)
## Cumulative: n.sims = 9000 iterations saved
##    mean    sd  2.5%   25%   50%   75% 97.5% 
## mu 6.72  0.49  5.80  6.39  6.71  7.04  7.73

### Note that mu can theoretically be outside the truncation interval,
### so we should allow the prior to extend above 8. Though we might
### expect sampling problems if mu is too far outside these bounds.



################################################################################
### 3.3 ZEROES TRICK

## 1. Truncated normal model for chickens data again
## Model code in models/tchick_zeros_mod_solution.txt

tchick.sim <- bugs(data=tchick.dat, inits=tchick.in, model="tchick_zeros_mod_solution.txt",
                   working.directory="models",
                   n.burnin=1000, n.iter=10000, n.chains=1,
                   param=c("mu"), DIC=FALSE)
print(tchick.sim, digits=2)

## Current: 1 chains, each with 10000 iterations (first 1000 discarded)
## Cumulative: n.sims = 9000 iterations saved
##  mean    sd  2.5%   25%   50%   75% 97.5% 
##  6.72  0.49  5.80  6.39  6.71  7.04  7.73 

## Should match the results from using T(), up to Monte Carlo error. 


## 2. Zeroes trick to implement censored data
## Model code in models/cchick_zerocens_mod_solution.txt

zcens.dat <- list(y = c(6.1, 6.3, 6.4, 7.1, 7.2, 7.3, NA, NA, NA))
zcens.in <- list(mu=7)

zcens.sim <- bugs(data=zcens.dat, inits=zcens.in, model="chick_zerocens_mod_solution.txt",
                  working.directory="models",
                  n.burnin=1000, n.iter=100000, n.chains=1,
                  param=c("mu","ypred"), DIC=FALSE)
print(zcens.sim, digits=2)

##       mean   sd 2.5%  25%  50%  75% 97.5%
## mu    7.36 0.35 6.68 7.12 7.36 7.59  8.04
## ypred 8.62 0.51 8.02 8.22 8.49 8.89  9.88

### equivalent model from slides in models/cchick_cens_slides_mod.txt
cens.dat <- list(y = c(6.1, 6.3, 6.4, 7.1, 7.2, 7.3, NA, NA, NA))
cens.in <- list(list(mu=7, y=c(NA,NA,NA, NA, NA, NA, 8.1, 8.1, 8.1)))
cens.sim <- bugs(data=cens.dat, inits=cens.in, model="cchick_cens_slides_mod.txt",
                  working.directory="models",
                  n.burnin=1000, n.iter=100000, n.chains=1,
                  param=c("mu","y[7]"), DIC=FALSE)
print(cens.sim, digits=2)

##      mean   sd 2.5%  25%  50%  75% 97.5%
## mu   7.36 0.35 6.68 7.12 7.36 7.60  8.05
## y[7] 8.61 0.51 8.02 8.22 8.49 8.89  9.89

## same results up to Monte Carlo error 

## Note: to see Monte Carlo error, need the coda package 
library(coda)
cens.coda <- as.mcmc.list(cens.sim) # convert R2OpenBUGS object to CODA format. 
summary(cens.coda) # "Time-series SE" is the Monte Carlo error 





################################################################################
### 3.4 SURVIVAL ANALYSIS

install.packages("JM")  ### For the AIDS data
library(JM)
library(survival)  ## For standard parametric survival models in R

### Basic survival model in R for the AIDS data
sr <- survreg(Surv(Time, death) ~ 1, data=aids.id, dist="weibull")
summary(sr)
##             Value Std. Error    z        p
## (Intercept)  3.23     0.0657 49.2 0.00e+00
## Log(scale)  -0.31     0.0674 -4.6 4.25e-06

## Fit the same model as Bayesian, in OpenBUGS
## JAGS/BUGS parameterisation is different : S(x) = exp(- lambda * x^nu)
## R survreg uses S(x) = exp(- (x/beta) ^ nu)
## where (1/beta)^nu = lambda

aids.data <- list(
    y = ifelse(aids.id$death, aids.id$Time, NA),
    c = ifelse(aids.id$death, 0, aids.id$Time),
    n = nrow(aids.id)
)
aids.in <- list(list(meansurv=10, loghr=0))

aids.sim <- bugs(data=aids.data,
                 inits=aids.in,
                 model="survival_mod.txt",
                 working.directory="models",
                 n.chains=1, n.burnin=1000, n.iter=10000,
                 parameters.to.save=c("lambda", "meansurv", "hr", "rintercept", "rlogscale"),
                 DIC=FALSE)
aids.sim

## Current: 1 chains, each with 10000 iterations (first 1000 discarded)
## Cumulative: n.sims = 9000 iterations saved
##            mean   sd 2.5%  25%  50%   75% 97.5%
## lambda      0.0  0.0  0.0  0.0  0.0   0.0   0.0
## meansurv   92.9 23.6 57.1 75.5 89.1 107.1 146.3
## hr          1.3  0.1  1.2  1.3  1.3   1.4   1.5
## rintercept  3.2  0.1  3.1  3.2  3.2   3.3   3.4
## rlogscale  -0.3  0.1 -0.5 -0.4 -0.3  -0.3  -0.2
