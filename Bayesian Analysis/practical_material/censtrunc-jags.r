library(rjags)

################################################################################
### 3.1. CENSORING


### Question 1(a). Rounded chickens implemented using interval censoring.

cchick.mod = "model {
    for( i in 1:9){
    y[i] ~ dnorm(mu,1)
    is.censored[i] ~ dinterval(y[i],ycut)
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


### Question 1(b).  Estimating the posterior SD

### Write model code etc. as exercise


round.cchick.mod = "model {
    for( i in 1:9){
    y[i] ~ dnorm(mu,1)
    yround[i] ~ dround(y[i],0)
    }
 mu ~ dunif(0, 20)
}"

round.cchick.dat <- list(
  yround= c(6,6,6,7,7,7,8,8,8)
)

## Initial values: each unknown y value initialised to the midpoint of the interval.
round.cchick.in <- list(mu=7, y=c(6,6,6,7,7,7,8,8,8))

round.cchick.jag <- jags.model(textConnection(round.cchick.mod), data=round.cchick.dat, 
                         inits=round.cchick.in)
update(round.cchick.jag, 1000) # burn-in
round.sam <- coda.samples(round.cchick.jag, var=c("y","mu"), n.iter=100000)
summary(round.sam)


### Question 2(a). Same model where the scales can not show values over 8

### Write model code etc. as exercise



################################################################################
### 3.2. TRUNCATION

### Truncated Normal model for the six chickens data using T() 

### Write model code etc. as exercise




################################################################################
### 3.3. ZEROES TRICK

## 1. Truncated normal model for chickens data again

## 2. Zeroes trick to implement censored data

### Write model code etc. as exercise
