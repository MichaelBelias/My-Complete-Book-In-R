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
cervix_inits <-
  list(psi = 0.5, beta0 = 0, beta = 0,
       phi = structure(.Data = c(0.5,0.5,0.5,0.5), .Dim = c(2,2)))

## full data
cervix_data <- c(as.list(dat[rep(1:12,count),]), n = sum(count))
## gold standard subset only
cervix_data <- c(as.list(dat[rep(1:8,count[1:8]),]), n = sum(count[1:8]))

cervix.jag <- jags.model(textConnection(cervix_mod), cervix_data, cervix_inits)

sam <- coda.samples(cervix.jag, c("beta","or.beta","phi"), 20000)
summary(sam)

### Just run on both datasets, compare the results and discuss.



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
mus <- sam2[,paste("mu[",1:27,"]",sep="")]
zmedian <- summ$quantiles[paste("z[",1:27,"]",sep=""),"50%"]
plot(zmedian, dat$y, xlab="Observed age (years)",
     ylab="Length (m)", pch=19, col="blue", xlim=c(0,35), ylim=c(1.5, 2.75))

## Exercise:
## Increase the prior measurement error SD (i.e. lower precision),
## run the model again, extracting the posterior median true ages,
## and overlay on the previous plot, something like this:
# points(zmedian, dat$y, col="red")


dugongs_mod2 <- "
model {
for(j in 1:n) {
y[j]       ~ dnorm(mu[j], tau)
mu[j]     <- alpha - beta*pow(gamma, z[j])
x[j]       ~ dnorm(z[j], 5)
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

inits2 <- list(alpha = 3, beta = 2, gamma = 0.9, log.sigma = -5)

dat2 <- list(x = c(1.0,  1.5,  1.5,  1.5, 2.5,   4.0,  5.0,  5.0,  7.0,
                  8.0,  8.5,  9.0,  9.5, 9.5,  10.0, 12.0, 12.0, 13.0,
                  13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5),
            y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
                  2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
                  2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57),
            n = 27)

dugongs2 <- jags.model(textConnection(dugongs_mod2), data=dat, inits=inits2)
sam2 <- coda.samples(dugongs2, c("mu","z","alpha","beta","gamma","sigma2"), 50000)
sam2 <- window(sam2, 1001, 50000) # discard burn-in
summ2 <- summary(sam2)

## Extract posterior median of true age
mus <- sam2[,paste("mu[",1:27,"]",sep="")]
zmedian2 <- summ2$quantiles[paste("z[",1:27,"]",sep=""),"50%"]
plot(zmedian, dat$y, xlab="Observed age (years)",
     ylab="Length (m)", pch=19, col="blue", xlim=c(0,35), ylim=c(1.5, 2.75))
points(zmedian2, dat2$y, col="red")



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

## Data which includes the calibration data relating observed exposures xcal to true exposures zcal
poll_data <- list(y = c(21, 20, 15), n = c(48, 34, 21),
                  x = c(10, 30, 50),
                  xcal = c(3, 6, 9, 10, 11, 13, 15, 20, 20.5, 21,
                           22, 23, 24, 24.5, 28, 30, 33, 41, 47, 53, 60, 70, 90),
                  zcal = c(1.1, 10.7, 3.8, 26.4, 15.8, 7, 20.3, 26.3, 25.2, 17.7, 34.8, 25.5,
                           17.1, 3.2, 35.9, 26.9, 29.4, 44.1, 47.6, 50.1, 58.4, 64.7, 73.6),
                  ncal = 23
                  )

## Exercise: adapt the model to include this data

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


