library(R2OpenBUGS)


##################################################
### 1.1. OUTCOMES MISSING NOT AT RANDOM
##################################################

## MMSE data
dat <- list(t = c(0, 5, 10, 15, 20),
            y = c(28, 26, 27, 25, NA),
            miss = c(0, 0 , 0, 0, 1) )
plot(dat$t, dat$y, ylim=c(0, 30)) # quick visualisation
ini <- list(list(alpha=27, beta=-0.1, sigma=2))   # OpenBUGS is more fussy than JAGS about having good initial values here 


### Part 1.  Priors

## [ code template in models/mmse_mod.txt ]


### Part 2.  R commands to run the model and monitor variables of interest

mmse.sim <- bugs(data=dat, inits=ini, model="models/mmse_mod.txt",
                 n.chains=2, n.burnin=1000, n.iter=10000,
                 param==c("sigma","alpha","beta","y[5]","p20"))
mmse.sim


### Part 3.  Adapt the code above to include a non-random missingness mechanism


### Part 4.  Change the normal to a t error distribution




##################################################
### 1.2.  MISSING COVARIATES
##################################################

## [ code template in models/malaria_mod.txt ]

### Add an imputation model for BEDNET to the code.

malaria <- read.table("malaria_data.txt", col.names=c("Y","AGE","BEDNET","GREEN","PHC"), skip=1, nrows=805)

mal.in <- list(list(alpha=0, beta=c(0,0,0,0), qmu=0, gamma=c(0,0,0)))


### Run model, monitoring and summarising variables indicated in the questions

mal.sim <- bugs(data=malaria, inits=mal.in, model="malaria_mod.txt",
                working.directory="models",
                n.chains=1, n.burnin=1000, n.iter=9000,
                param=c("alpha","or","beta","qmu", "gamma","BEDNET[1:10]","BEDNET[531:540]"), DIC=FALSE)

print(mal.sim, digits=4)
## note in R2OpenBUGS typing the name of the object just shows the quantiles.  To get the posterior mean and sd, explicitly call the print method like this.
