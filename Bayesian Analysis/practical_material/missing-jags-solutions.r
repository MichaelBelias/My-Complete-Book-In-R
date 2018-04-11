library(rjags)


##################################################
### 1. OUTCOMES MISSING NOT AT RANDOM
##################################################

## MMSE data
dat <- list(t = c(0, 5, 10, 15, 20),
            y = c(28, 26, 27, 25, NA) )
plot(dat$t, dat$y, ylim=c(0, 30)) # quick visualisation
ini <- list(alpha=20, beta=-10, sigma=1)


### Part 1.  Priors

### alpha: MMSE is bounded on this range, so the intercept (equal to the baseline value)  is bounded by definition on this range.  We have no other prior information to restrict it, so a uniform is reasonable.   Alternatively we could use a distribution truncated on this range, or a linearly-translated beta distribution, if we wanted a peak at somewhere between 0 and 30.

### beta: beta is the expected change in MMSE after one year.  "97.5% certain that MMSE cannot increase over time, and 97.5% certain that after 1 year it cannot drop more than 20 points" implies a prior mean halfway between 0 and -20 for the expected change after one year.
### Then interpreting the range of (-20, 0) as the 2.5% and 97.5% quantiles, or +- 2 standard deviations, the prior standard deviation is then the range width divided by 4, equal to 5 units of MMSE, equivalently a prior precision of 1/(5*5) = 0.04.

### sigma:  as noted in the question sheet, the normal error model is not strictly realistic given that MMSE is bounded on (0,30).  Even so, we do not want to allow a-priori standard deviations for MMSE which permit MMSE values outside (0,30) with a large probability.  Thus it is unlikely to be much over 10.  A uniform is reasonable, alternatively we might use a positive distribution (sugh as the exponential or gamma) for the SD, variance or precision, concentrated on values below about 10.

mmse.mod <- "
model {
 for (i in 1:5){
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha + beta*t[i]
 }
 p20 <- step(20 - y[5])

 ### PRIOR DISTRIBUTIONS
 alpha ~ dunif(0, 30)
 beta ~ dnorm(-10, 0.04)
 sigma ~ dunif(0, 10)

 tau <- 1/(sigma*sigma)
}
"

### Part 2.  rjags commands to run the model and monitor variables of interest

mmse.jag <- jags.model(textConnection(mmse.mod), dat, ini)
sam <- coda.samples(mmse.jag, var=c("sigma","alpha","beta","y[5]","p20"), n.iter=10000)
sam <- window(sam, 1001, 10000)
summary(sam)
plot(sam, ask=TRUE)


## The posterior distributions are concentrated compared to their priors, but they are all likely to be sensitive to the exact choice of prior, given only 4 data points.
## For sigma, the posterior is particularly skewed, so that the estimate of the 97.5 quantile is unstable, requiring many more iterations to estimate precisely.


### Part 3.  Code including a non-random missingness mechanism

### Needed to add the missingness model, priors for the parameters in
### the missingness model, and to add the missing data indicator to
### the data list.  Also sensible to add initial values for the
### missingness model parameters.

mmse.miss.mod <- "
model {
 for (i in 1:5){
  y[i] ~ dnorm(mu[i], tau)
#  y[i] ~ dt(mu[i], tau, 4) # alternative t errors for Part 4
  mu[i] <- alpha + beta*t[i]

  ## Model for missingness
  miss[i] ~ dbern(p[i])
  logit(p[i]) <- a - b*y[i]/5 # OR doubled (b=log(2)) if MMSE 5 units lower
 }
 p20 <- step(20 - y[5])

 alpha ~ dunif(0, 30)
 beta ~ dnorm(-10, 0.04)
 sigma ~ dunif(0, 10)
 tau <- 1/(sigma*sigma)

 ## Priors for missingness model parameters
 a ~ dlogis(0, 1)
 b <- log(2) #
}
"

### Add a new variable "miss" to the data
dat <- list(t = c(0, 5, 10, 15, 20),
            y = c(28, 26, 27, 25, NA),
            miss = c(0, 0 , 0, 0, 1) )
### Extend the initial values
ini <- list(alpha=20, beta=-10, sigma=1, a=0)

mmse.jag <- jags.model(textConnection(mmse.miss.mod), dat, ini)
sam.miss <- coda.samples(mmse.jag, var=c("sigma","alpha","beta","y[5]","p20"), n.iter=100000)
sam.miss <- window(sam.miss, 1001, 100000)
summary(sam.miss)

##          Mean     SD  Naive SE Time-series SE
## alpha 27.3940 1.8332 0.0058264       0.020827
## beta  -0.2000 0.2173 0.0006905       0.002236
## p20    0.1858 0.3890 0.0012363       0.005049
## sigma  2.4073 1.9902 0.0063252       0.034163
## y[5]  22.4297 5.2977 0.0168370       0.074133

## 2. Quantiles for each variable:
##          2.5%     25%     50%     75%   97.5%
## alpha 22.5335 26.8278 27.7042 28.4711 29.7456
## beta  -0.6994 -0.2857 -0.1874 -0.1004  0.2199
## p20    0.0000  0.0000  0.0000  0.0000  1.0000
## sigma  0.5790  1.0656  1.6708  2.9528  8.3580
## y[5]   6.9667 21.2529 23.6052 25.1362 29.5126

## Note the estimate of beta is lower, that is the MMSE decline is
## estimated to be quicker.

## The estimated probability of the missing measurement being below 20
## has increased from 0.076 to 0.19.


## Part 4: t-distributed errors.

## Simply comment out the normally-distributed errors in the previous
## model, uncomment the line with t-distributed errors, and re-run.

## With normal errors: p(y < 20) = 0.19, and beta estimated as -0.19 (-0.72, 0.23)

## With t distributed errors -- much lower values of MMSE still are
## plausible: we estimate p(y < 20) = 0.24 after 100000 samples, though
## beta is substantively the same.





##################################################
### 1.2.  MISSING COVARIATES
##################################################

### Add an imputation model for BEDNET to the code.

malaria <- read.table("malaria_data.txt", col.names=c("Y","AGE","BEDNET","GREEN","PHC"), skip=1, nrows=805)

mal.mod <- "
model{
   for(i in 1:805) {
      Y[i] ~ dbern(p[i])
      logit(p[i]) <- alpha + beta[1]*(AGE[i] - mean(AGE[])) + beta[2]*BEDNET[i] +
                                    beta[3]*(GREEN[i] - mean(GREEN[])) + beta[4]*PHC[i]
     BEDNET[i] ~ dbern(q[i])
     logit(q[i]) <- qmu + gamma[1]*AGE[i] + gamma[2]*GREEN[i] + gamma[3]*PHC[i]
   }
   # vague priors on regression coefficients of analysis model
   alpha ~ dlogis(0, 1)
   for (i in 1:4){
    beta[i] ~ dt(0, 0.16, 1)
    or[i] <- exp(beta[i])
   }
   qmu ~ dlogis(0,1)
   for (i in 1:3){
     gamma[i] ~ dnorm(0, 0.1)  # fairly arbitrary choice of prior here -- sensitivity analysis advised!
   }
}
"

mal.in <- list(alpha=0, beta=c(0,0,0,0), qmu=0, gamma=c(0,0,0))
mal.jag <- jags.model(textConnection(mal.mod), malaria, mal.in)
sam <- coda.samples(mal.jag, c("alpha","or","beta","qmu", "gamma","BEDNET[1:10]","BEDNET[531:540]"), n.iter=20000)
sam <- window(sam, 1001, 20000)
summary(sam)


## gives quantiles:
##                2.5%      25%      50%        75%    97.5%
## or[1]        1.1313  1.24062  1.30121  1.3626021  1.49158
## or[2]        0.7105  0.96046  1.11920  1.3052839  1.74721
## or[3]        0.9405  0.96831  0.98384  0.9992591  1.03098
## or[4]        0.4286  0.53722  0.60455  0.6795076  0.84905

## The OR for AGE is estimated as 1.30 (1.14, 1.49) after accounting for missing covariates by imputation, so it appears that the complete-case analysis was biased.  The estimates are also more precise after including the partially-observed data.  For the other coefficients, the same behaviour is seen to a lesser degree.

## BEDNET[1]    0.80763
## BEDNET[2]    0.81600
## BEDNET[3]    0.81884
## BEDNET[4]    0.81400
## BEDNET[5]    0.80979
## BEDNET[6]    0.80668
## BEDNET[7]    0.80389
## BEDNET[8]    0.81979
## BEDNET[9]    0.81926
## BEDNET[10]   0.81726
##  ...
## BEDNET[531]  0.56005
## BEDNET[532]  0.53889
## BEDNET[533]  0.56658
## BEDNET[534]  0.53763
## BEDNET[535]  0.53842
## BEDNET[536]  0.53432
## BEDNET[537]  0.56579
## BEDNET[538]  0.54742
## BEDNET[539]  0.54426
## BEDNET[540]  0.54274

## The posterior means for the imputed BEDNET values in rows 1-10 and 531-540 are shown above.

## Firstly we note that they are on average 0.8 in rows 1-10 and 531-540 because the value of PHC is all 1 and all 0 in these blocks respectively.  PHC is a particularly strong predictor of BEDNET, as shown by the estimates of the imputation model coefficients:

##                 Mean      SD
## gamma[1]    -0.03075 0.09950 (AGE)
## gamma[2]    -0.08256 0.04125 (GREEN)
## gamma[3]     1.24217 0.21412 (PHC)

## Secondly: we might expect the mean imputed values for BEDNET[1] and BEDNET[7] to be slightly lower than their neighbours, because Y=0 for these cases, and BEDNET has a positive association with Y.  Conversely the posterior means of BEDNET[533] and BEDNET[537] might be expected to be higher than their neighbours, since Y=1 in these cases.   These patterns are seen to some extent, but are not so striking since the association of BEDNET with Y is weak.
