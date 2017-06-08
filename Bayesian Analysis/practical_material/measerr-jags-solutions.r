library(rjags)


##################################################
### QUESTION 1:  BINARY MISCLASSIFICATION

cervix_mod <- "
model {
  for (i in 1:n) {
    d[i]         ~ dbern(p[i])
    logit(p[i]) <- beta0 + beta*z[i]
    z[i]         ~ dbern(psi)
    x[i]         ~ dbern(phi[z[i]+1, d[i]+1])
  }
  for (j in 1:2) {
    for (k in 1:2) {
      phi[j, k]  ~ dunif(0, 1)
    }
  }
  psi            ~ dunif(0, 1)
  beta0          ~ dlogis(0, 1)
  beta           ~ dt(0, 0.16, 1) # Cauchy with scale 2.5 (Gelman 2008)
  or.beta        <- exp(beta)
}
"

dat <- data.frame(
d = c(1,1,1,1,0,0,0,0,1,1,0,0),
z = c(0,0,1,1,0,0,1,1,NA,NA,NA,NA),
x = c(0,1,0,1,0,1,0,1,0,1,0,1)
)
count = c(13,3,5,18,33,11,16,16,318,375,701,535)

## full data
cervix_data <- c(as.list(dat[rep(1:12,count),]), n = sum(count))
## gold standard subset only
cervix_data <- c(as.list(dat[rep(1:8,count[1:8]),]), n = sum(count[1:8]))

cervix_inits <-
  list(psi = 0.5, beta0 = 0, beta = 0,
       phi = structure(.Data = c(0.5,0.5,0.5,0.5), .Dim = c(2,2)))

cervix.jag <- jags.model(textConnection(cervix_mod), cervix_data, cervix_inits)
sam <- coda.samples(cervix.jag, c("beta","or.beta","phi"), 20000)
summary(sam)


### With complete data only: the posterior precision of the OR is not much worse: 95% CI (0.89, 4.11).  Though the measurement error is big, and the precision improvements would have been greater if the measurement error were smaller.

### While we do get slightly more precise estimates from the full model, the potential disadvantages are the assumptions that it relies on.  Specifically, the assumption that the complete cases have the same parameter values as the incomplete cases.  Particularly the misclassification probabilities and odds ratios.   These assumptions are uncheckable from data.


##################################################
### QUESTION 2:  CLASSICAL MEASUREMENT ERROR (dugongs)

dugongs_mod <- "
model {
  for(j in 1:n) {
    y[j]       ~ dnorm(mu[j], tau)
    mu[j]     <- alpha - beta*pow(gamma, z[j])
    x[j]       ~ dnorm(z[j], 1)
    z[j]       ~ dunif(0, 100)
  }
  alpha        ~ dunif(0, 100)
  beta         ~ dunif(0, 100)
  gamma        ~ dunif(0, 1)
  tau         <- 1/sigma2
  log(sigma2) <- 2*log.sigma
  log.sigma    ~ dunif(-10, 10)
  for (j in 1:n) {resx[j] <- x[j] - z[j]}
}
"

inits <- list(alpha = 3, beta = 2, gamma = 0.9, log.sigma = -5)

dat <- list(x = c(1.0,  1.5,  1.5,  1.5, 2.5,   4.0,  5.0,  5.0,  7.0,
                  8.0,  8.5,  9.0,  9.5, 9.5,  10.0, 12.0, 12.0, 13.0,
                  13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5),
            y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
                  2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
                  2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57),
            n = 27)

dugongs <- jags.model(textConnection(dugongs_mod), data=dat, inits=inits)
sam <- coda.samples(dugongs, c("mu","z","alpha","beta","gamma","sigma2"), 50000)
sam <- window(sam, 1001, 50000) # discard burn-in
summ <- summary(sam)

## Extract posterior median of true age
zmedian <- summ$quantiles[paste("z[",1:27,"]",sep=""),"50%"]
plot(zmedian, dat$y, xlab="Observed age (years)",
     ylab="Length (m)", pch=19, col="blue", xlim=c(0,35), ylim=c(1.5, 2.75))
# points(zmedian, dat$y, col="red") # after increasing the error...

## Solution:
## You should find that after increasing the measurement error, the dugong length is a smooth nonlinear function of the fitted true ages z, with less noise than in the original model.   The reason is that the estimates of z are now influenced much more by their likelihood contribution from the length data y through the nonlinear model, compared to the likelihood contribution from the observed ages x.   Increasing the measurement error makes the information provided by x about z weaker.  Also notice that the uncertainty around the fitted true ages increases when the measurement error is increased (compare the posterior standard deviations of z).



##################################################
### QUESTION 3:  BERKSON MEASUREMENT ERROR (pollution)

### As presented in the lectures, without calibration data
poll_mod <- "
model{
 for (j in 1:3) {
   y[j]         ~ dbin(p[j], n[j])
   logit(p[j]) <- theta[1] + theta[2]*z[j]
   z[j]         ~ dnorm(mu[j], 0.01232)
   mu[j]       <- alpha + beta*x[j]
 }
 theta[1]       ~ dlogis(0, 1)
 theta[2]       ~ dnorm(0, 0.2)
 or10 <- exp(theta[2]*10)
}
"

poll_data <- list(y = c(21, 20, 15), n = c(48, 34, 21),
                  x = c(10, 30, 50), alpha = 4.48, beta = 0.76)
poll_inits <- list(theta = c(0.0, 0.0), z = c(0.0, 0.0, 0.0))
polljag <- jags.model(textConnection(poll_mod), poll_data, poll_inits)
sam <- coda.samples(polljag, var=c("theta","or10"), n.iter=100000)
summary(sam)


## Solution: model including calibration data relating observed exposures xcal to true exposures zcal

poll_cal_mod <- "
model{
 for (j in 1:3) {
   y[j]         ~ dbin(p[j], n[j])
   logit(p[j]) <- theta[1] + theta[2]*z[j]
   z[j]         ~ dnorm(mu[j], prec)
   mu[j]       <- alpha + beta*x[j]
 }
 for (j in 1:ncal) {
   zcal[j] ~ dnorm(mucal[j], prec)
   mucal[j] <- alpha + beta*xcal[j]
 }
 theta[1]       ~ dlogis(0, 1)
 theta[2]       ~ dnorm(0, 0.2)
 or10 <- exp(theta[2]*10)
 alpha ~ dnorm(0, 0.001)
 beta ~ dnorm(0, 0.001)
 sig ~ dunif(0, 100)
 prec <- 1/(sig*sig)
}
"

poll_data <- list(y = c(21, 20, 15), n = c(48, 34, 21),
                  x = c(10, 30, 50),
                  xcal = c(3, 6, 9, 10, 11, 13, 15, 20, 20.5, 21,
                           22, 23, 24, 24.5, 28, 30, 33, 41, 47, 53, 60, 70, 90),
                  zcal = c(1.1, 10.7, 3.8, 26.4, 15.8, 7, 20.3, 26.3, 25.2, 17.7, 34.8, 25.5,
                           17.1, 3.2, 35.9, 26.9, 29.4, 44.1, 47.6, 50.1, 58.4, 64.7, 73.6),
                  ncal = 23
                  )

poll_inits <- list(theta = c(0.0, 0.0), alpha=0.2, beta=1, sig=1)
polljag <- jags.model(textConnection(poll_cal_mod), poll_data, poll_inits)
sam <- coda.samples(polljag, var=c("theta","or10","alpha","beta","sig"), n.iter=100000)
summary(sam)

### Including the calibration data does not make much difference to the posterior distribution of the odds ratio, compared to using known error coefficients alpha,beta,sigma.  (which are similar to the posterior estimates of alpha, beta and sigma from the data).

## If anything, the estimates appear more precise with the calibration data included - it is unclear why - perhaps because the calibration model permits the chance of smaller (as well as larger) measurement errors, and smaller errors would lead to more precise estimates.

## Though the point of this exercise is just to get used to coding the evidence synthesis in BUGS!
