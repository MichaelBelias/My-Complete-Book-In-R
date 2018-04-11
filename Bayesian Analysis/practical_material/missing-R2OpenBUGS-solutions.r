library(R2OpenBUGS)

##################################################
### 1. OUTCOMES MISSING NOT AT RANDOM
##################################################

## MMSE data
dat <- list(t = c(0, 5, 10, 15, 20),
            y = c(28, 26, 27, 25, NA) )
plot(dat$t, dat$y, ylim=c(0, 30)) # quick visualisation
ini <- list(list(alpha=27, beta=-0.1, sigma=2))   # OpenBUGS is more fussy than JAGS about having good initial values here


### Part 1.  Priors
### See code in models/mmse_mod_solution.txt

### alpha: MMSE is bounded on this range, so the intercept (equal to the baseline value)  is bounded by definition on this range.  We have no other prior information to restrict it, so a uniform is reasonable.   Alternatively we could use a distribution truncated on this range, or a linearly-translated beta distribution, if we wanted a peak at somewhere between 0 and 30.

### beta: beta is the expected change in MMSE after one year.  "97.5% certain that MMSE cannot increase over time, and 97.5% certain that after 1 year it cannot drop more than 20 points" implies a prior mean halfway between 0 and -20 for the expected change after one year.
### Then interpreting the range of (-20, 0) as the 2.5% and 97.5% quantiles, or +- 2 standard deviations, the prior standard deviation is then the range width divided by 4, equal to 5 units of MMSE, equivalently a prior precision of 1/(5*5) = 0.04.

### sigma:  as noted in the question sheet, the normal error model is not strictly realistic given that MMSE is bounded on (0,30).  Even so, we do not want to allow a-priori standard deviations for MMSE which permit MMSE values outside (0,30) with a large probability.  Thus it is unlikely to be much over 10.  A uniform is reasonable, alternatively we might use a positive distribution (sugh as the exponential or gamma) for the SD, variance or precision, concentrated on values below about 10.


### Part 2.  R commands to run the model and monitor variables of interest

mmse.sim <- bugs(data=dat, inits=ini, model="mmse_mod_solution.txt",
                 working.directory="./models",
                 n.chains=1, n.burnin=1000, n.iter=10000,
                 param=c("sigma","alpha","beta","y[5]","p20"))
print(mmse.sim, digits=2)

## Current: 1 chains, each with 10000 iterations (first 1000 discarded)
## Cumulative: n.sims = 9000 iterations saved
##           mean   sd  2.5%   25%   50%   75% 97.5%
## sigma     2.14 1.74  0.58  1.02  1.52  2.53  7.57
## alpha    27.27 1.95 22.08 26.74 27.62 28.35 29.69
## beta     -0.13 0.22 -0.47 -0.24 -0.16 -0.06  0.42
## y[5]     24.68 4.19 16.69 22.90 24.50 26.21 34.09
## p20       0.08 0.26  0.00  0.00  0.00  0.00  1.00
## deviance 14.59 4.56  8.67 11.09 13.56 17.09 25.63

## The posterior distributions are concentrated compared to their priors, but they are all likely to be sensitive to the exact choice of prior, given only 4 data points.
## For sigma, the posterior is particularly skewed, so that the estimate of the 97.5 quantile is unstable, requiring many more iterations to estimate precisely.


### Part 3.  Code including a non-random missingness mechanism

### Needed to add the missingness model, priors for the parameters in
### the missingness model, and to add the missing data indicator to
### the data list.  Also sensible to add initial values for the
### missingness model parameters.


### Add a new variable "miss" to the data
dat <- list(t = c(0, 5, 10, 15, 20),
            y = c(28, 26, 27, 25, NA),
            miss = c(0, 0 , 0, 0, 1) )
### Extend the initial values
ini <- list(list(alpha=27, beta=-0.1, sigma=2, a=0))

mmse.sim <- bugs(data=dat, inits=ini, model="mmse_miss_mod_solution.txt",
                 working.directory="./models",
                 n.chains=1, n.burnin=1001, n.iter=100000,
                 param=c("sigma","alpha","beta","y[5]","p20"))
print(mmse.sim, digits=2)

## Current: 1 chains, each with 1e+05 iterations (first 1001 discarded)
## Cumulative: n.sims = 98999 iterations saved
##           mean   sd  2.5%   25%   50%   75% 97.5%
## sigma     2.42 2.00  0.58  1.07  1.68  2.98  8.42
## alpha    27.38 1.79 22.68 26.79 27.68 28.45 29.75
## beta     -0.20 0.22 -0.71 -0.29 -0.19 -0.10  0.22
## y[5]     22.38 5.38  6.91 21.22 23.61 25.13 29.47
## p20       0.19 0.39  0.00  0.00  0.00  0.00  1.00
## deviance 20.76 4.75 13.89 17.04 19.92 23.82 31.29

## Note the estimate of beta is lower, that is the MMSE decline is
## estimated to be quicker.

## The estimated probability of the missing measurement being below 20
## has increased from 0.08 to 0.19.


## Part 4: t-distributed errors.

## Simply comment out the normally-distributed errors in the previous
## model, uncomment the line with t-distributed errors, and re-run.

## With normal errors: p(y < 20) = 0.19, and beta estimated as -0.19 (-0.72, 0.23)

## With t distributed errors -- much lower values of MMSE still are
## plausible: we estimate p(y < 20) = 0.22 after 100000 samples, though
## beta is substantively the same.





##################################################
### 1.2.  MISSING COVARIATES
##################################################

### 1. Add an imputation model for BEDNET to the code.

malaria <- read.table("malaria_data.txt", col.names=c("Y","AGE","BEDNET","GREEN","PHC"), skip=1, nrows=805)

mal.in <- list(list(alpha=0, beta=c(0,0,0,0), qmu=0, gamma=c(0,0,0)))

mal.sim <- bugs(data=malaria, inits=mal.in, model="malaria_mod_solution.txt",
                working.directory="models",
                n.chains=1, n.burnin=1000, n.iter=9000,
                param=c("alpha","or","beta","qmu", "gamma","BEDNET[1:10]","BEDNET[531:540]"), DIC=FALSE)

## This is relatively slow compared to other examples -- may take a couple of minutes to run.
## If in doubt, reduce n.burnin and n.iter to smaller value and do a test run to predict how long a full 10000 will take

print(mal.sim, digits=4) ## note in R2OpenBUGS typing the name of the object just shows the quantiles.  To get the posterior mean and sd, explicitly call the print method like this

##                 mean     sd      2.5%       25%       50%       75%     97.5%
## or[1] (AGE)    1.3044 0.0908    1.1390    1.2410    1.3000    1.3640    1.4910
## or[2] (BEDNET) 1.1361 0.2702    0.7099    0.9429    1.0990    1.2940    1.7460
## or[3] (GREEN)  0.9828 0.0230    0.9385    0.9671    0.9823    0.9980    1.0290
## or[4] (PHC)    0.6059 0.1062    0.4260    0.5310    0.5976    0.6698    0.8373

## The OR for AGE is estimated as 1.30 (1.14, 1.49) after accounting for missing covariates by imputation, so it appears that the complete-case analysis was biased.  The estimates are also more precise after including the partially-observed data.  For the other coefficients, the same behaviour is seen to a lesser degree.

## BEDNET[1]      0.8029
## BEDNET[2]      0.8336
## BEDNET[3]      0.8081
## BEDNET[4]      0.8212
## BEDNET[5]      0.8080
## BEDNET[6]      0.7975
## BEDNET[7]      0.8144
## BEDNET[8]      0.8280
## BEDNET[9]      0.8141
## BEDNET[10]     0.8094
## BEDNET[531]    0.5530
## BEDNET[532]    0.5331
## BEDNET[533]    0.5671
## BEDNET[534]    0.5356
## BEDNET[535]    0.5352
## BEDNET[536]    0.5240
## BEDNET[537]    0.5686
## BEDNET[538]    0.5622
## BEDNET[539]    0.5375
## BEDNET[540]    0.5450

## The posterior means for the imputed BEDNET values in rows 1-10 and 531-540 are shown above.

## Firstly we note that they are on average 0.8 in rows 1-10 and 531-540 because the value of PHC is all 1 and all 0 in these blocks respectively.  PHC is a particularly strong predictor of BEDNET, as shown by the estimates of the imputation model coefficients:

## gamma[1]      -0.0300 0.1016 (AGE)
## gamma[2]      -0.0759 0.0167 (GREEN)
## gamma[3]       1.2478 0.2044 (PHC)

## Secondly: we might expect the mean imputed values for BEDNET[1] and BEDNET[7] to be slightly lower than their neighbours, because Y=0 for these cases, and BEDNET has a positive association with Y.  Conversely the posterior means of BEDNET[533] and BEDNET[537] might be expected to be higher than their neighbours, since Y=1 in these cases.   These patterns are seen to some extent, but are not so striking since the association of BEDNET with Y is weak.

## Note also the Monte Carlo SE of the BEDNET means is about 0.006
library(coda)
mal.coda <- as.mcmc.list(mal.sim)
summary(mal.coda)
