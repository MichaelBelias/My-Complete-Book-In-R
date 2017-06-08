library(rjags)


##################################################
### 1.1. OUTCOMES MISSING NOT AT RANDOM 
##################################################

## MMSE data 
dat <- list(t = c(0, 5, 10, 15, 20),
            y = c(28, 26, 27, 25, NA) )
plot(dat$t, dat$y, ylim=c(0, 30)) # quick visualisation 
ini <- list(alpha=20, beta=-10, sigma=1)


### Part 1.  Priors 

mmse.mod <- "
model {
 for (i in 1:5){
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- alpha + beta*t[i]
 }

 p20 <- step(20 - y[5])

 ### INSERT PRIOR DISTRIBUTIONS HERE
 alpha ~ dunif(-20,20) 
 beta ~  dnorm(-10, 10)
 sigma ~ dunif( 0,10)  

 tau <- 1/(sigma*sigma)
}
"

### Part 2.  rjags commands to run the model and monitor variables of interest

mmse.jag <- jags.model(textConnection(mmse.mod), dat, ini)
sam <- coda.samples(mmse.jag, var=c("sigma","alpha","beta","y[5]","p20"), n.iter=10000)
sam <- window(sam, 1001, 10000) # discard burn-in (convergence assumed before 1000) 
summary(sam)
dev.new()
plot(sam, ask=TRUE)


### Part 3.  Adapt the code above to include a non-random missingness mechanism


### Part 4.  Change the normal to a t error distribution 




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

### INSERT IMPUTATION MODEL HERE 
    BEDNET[i] ~ dbern(q[i])
    logit(q[i]) <- gamma[1] + gamma[2]*AGE[i] +gamma[3]*GREEN[i] + gamma[4]*PHC[i]
   }

   # vague priors on regression coefficients of analysis model
   alpha ~ dlogis(0, 1)
   for (i in 1:4){
    beta[i] ~ dt(0, 0.16, 1)
    or[i] <- exp(beta[i])
   }

### PRIORS FOR IMPUTATION MODEL COEFFICIENTS HERE 
   for (i in 1:4){
    gamma[i] ~ dnorm(0, 1)
   }
}
"

mal.in <- list(alpha=0, beta=c(0,0,0,0), gamma=c(0,0,0,0))


### Run model, monitoring and summarising variables indicated in the questions 

mal.jag <- jags.model(textConnection(mal.mod), malaria, mal.in)
sam <- coda.samples(mal.jag, c("alpha","or","beta","gamma","BEDNET[1:10]","BEDNET[531:540]"), n.iter=10000)
traceplot(sam[,c("beta[1]","beta[2]","beta[3]","beta[4]")])
summary(sam)
