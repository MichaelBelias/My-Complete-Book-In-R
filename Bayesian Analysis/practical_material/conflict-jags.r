library(rjags)

## Set of functions for processing posterior samples in R
## including functions for calculating deviance summaries
source("./fns.R")

#############################################################
## SEPSIS META-ANALYSIS MODEL                              ##
## QUESTION 1 - cross-validatory mixed-predictive p-values ##
#############################################################

## Leave-One-Out cross-validation for Sepsis model
sepCV_model <- "
model
{
  ## Cross-validation: repeat data set leaving one out each time
  for(i in 1:Ns)
  {
    ## For each study, Ns = total number of studies
    for(j in 1:Ns)
    {
      ## for each arm (Na = 2 for standard MA)
      for(k in 1:Na)
      {
        ## Binomial likelihood
        ycv[i,j,k] ~ dbin(p[i,j,k], ncv[i,j,k])
        
        ## on logit scale, proportion is 
        ## probability of success in terms of study baselines mu and
        ## study-specific treatment contrasts delta (log odds ratios,
        ## relative to study baseline), if not left-out
        logit(p[i,j,k]) <- ((1 - equals(i,j)) * (mu[i,j] + delta[i,j,k]))
      }
      
      ## for each non-baseline arm
      for(k in 2:Na)
      {
        ## contrasts (log odds ratios, whether direct or indirect)
        ## are independent effects
        delta[i,j,k] ~ dnorm(d[i,k], prec.d[i])
      }
      
      ## study-specific baselines, vague priors
      mu[i,j] ~ dnorm(0, 0.01)
      
      ## contrasts at baseline are 0
      delta[i,j,1] <- 0
    }
    
    ## for each arm
    for(k in 1:Na)
    {
      ## logit proportions for left-out study equal to MLE
      loo[i,k] <- (ycv[i,i,k] / ncv[i,i,k])
      delta.loo[i,k] <- logit(loo[i,k]) - logit(loo[i,1])
      
      ## replicated data for left-out study
      yrep[i,k] ~ dbin(prep[i,k], ncv[i,i,k])
      logit(prep[i,k]) <- mu[i,i] + delta[i,i,k]
      
      ## cross-validatory mixed-predictive p-values
      ## one-sided mid-p-value, testing yrep >= ycv for both treatment arms
      pval.cvmix[i,k] <- ??
    }
    
    ## for each non-baseline arm
    for(k in 2:Na)
    {
      ## cross-validatory mixed-predictive p-value
      ## one-sided, testing logit(delta) > logit(delta.loo)
      ## (not integer so not a mid-p-value)
      pval.cvmix.delta[i,k] <- ??
    }
    
    ## Priors for basic parameter (log odds ratio of treatment vs control)
    ## d[i,1] is overall baseline treatment (as opposed to study-specific baselines)
    d[i,1] <- 0
    for(t in 2:Nt)
    {
      d[i,t] ~ dnorm(0, 0.01)
    }
    prec.d[i] <- 1 / (sd.d[i] * sd.d[i])
    sd.d[i] ~ dunif(0,10)
  }
}
"

## data list
(sep_dat <- list(
  Ns = 10,
  Nt = 2,
  Na = 2,
  ## copy data Ns times, once for each cross-validation
  ycv  = aperm(structure(.Data = rep(c(
    23,   20,
    8,    2,
    5,    0,
    14,    8,
    209,  186,
    5,    4,
    13,   10,
    13,   19,
    8,    3,
    39,   40
  ), 10), .Dim = c(2,10,10)), c(3,2,1)),
  ncv  = aperm(structure(.Data = rep(c(
    65,   61,
    43,   43,
    59,   56,
    32,   34,
    1212, 1204,
    50,  100,
    34,   68,
    41,   40,
    40,   40,
    381,  372
  ), 10), .Dim = c(2,10,10)), c(3,2,1))
))


## Inits
initsSep<- list(
  list(d = t(structure(.Data = rep(c(NA,0),10),
                       .Dim = c(2,10))),
       mu = t(structure(.Data = rep(0,100),
                        .Dim = c(10,10)))
  ),
  list(d = t(structure(.Data = rep(c(NA,0.1),10),
                       .Dim = c(2,10))),
       mu = t(structure(.Data = rep(c(1,-1,-2,0,0,-2,1,0,2,2),10),
                        .Dim = c(10,10)))
  )
)


## Parameters to monitor
paramsSep <- c("p","mu","delta","d","sd.d","pval.cvmix","pval.cvmix.delta")


## Initialise model
sep.jm <- jags.model(textConnection(sepCV_model),
                     data = sep_dat,
                     inits = initsSep,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(sep.jm, n.iter = 1000)
## samples to keep
outSep <- coda.samples(sep.jm,
                       variable.names = paramsSep,
                       n.iter = 9000,
                       n.thin = 3)
print(proc.time() - ptm)

## check traces

## summary stats

## cross-validatory mixed-predictive p-values, comparing yrep to ycv:

## compared to cross-validatory mixed-predictive p-values, comparing delta to 
## T(ycv[i,]) = logit(ycv[i,2]/ncv[i,2]) - logit(ycv[i,1]/ncv[i,1]):

## plot comparisons for Sandberg and Fanaroff studies


#########################################################################################


############################################
## SEPSIS META-ANALYSIS MODEL             ##
## QUESTION 2 - systematic bias modelling ##
############################################



## RANDOM EFFECTS + BIAS MODEL
sepBias_model <- "
model
{
    ## For each study, Ns = total number of studies
    for(j in 1:Ns)
        {
            ## for each arm (Na = 2 for standard MA)
            for(k in 1:Na)
                {
                    ## Binomial likelihood
                    y[j,k] ~ dbin(p[j,k], n[j,k])

                    ## For model criticism - mean of sampling distribution and deviances
                    yhatp[j,k] <- p[j,k] * n[j,k]
                    devp[j,k] <- 2 * (y[j,k] * (log(y[j,k])-log(yhatp[j,k])) + (n[j,k]-y[j,k]) * (log(n[j,k]-y[j,k]) - log(n[j,k]-yhatp[j,k])))

                    ## probability of success in terms of study baselines mu and
                    ## study-specific treatment contrasts delta (log odds ratios,
                    ## relative to study baseline) + study-specific bias adjustment
                    logit(p[j,k]) <- mu[j] + deltaB[j,k]
               }

            ## for each non-baseline arm
            for(k in 2:Na)
                {
                    ## contrasts (log odds ratios, whether direct or indirect)
                    ## are independent effects
                    deltaB[j,k] <- delta[j,k] + bias[j]
                    delta[j,k] ~ dnorm(d[k], prec.d)
                }
                    
            ## study-specific baselines, vague priors
            mu[j] ~ dnorm(0, 0.01)

            ## contrasts at baseline are 0
            deltaB[j,1] <- 0
        }
  
    ## Priors for basic parameter (log odds ratio of treatment vs control)
    ## d[1] is overall baseline treatment (as opposed to study-specific baselines)
    d[1] <- 0
    for(i in 2:Nt)
        {
            d[i] ~ dnorm(0, 0.01)
        }
    prec.d <- 1 / (sd.d * sd.d)
    sd.d ~ dunif(0,10)


    ## informative priors for studies with at least one high risk of bias
    ## otherwise no bias
    for(j in 1:Nsafe)
        {
            bias[safe[j]] <- 0
        }
    for(j in 1:Nrisky)
        {
            bias[risky[j]] ~ dnorm(bias.mu, bias.prec)
        }
}
"



## data list
(sepBias_dat <- list(
     Ns = 10,
     Nt = 2,
     Na = 2,
     y  = t(structure(.Data = c(
                           23,   20,
                            8,    2,
                            5,    0,
                           14,    8,
                          209,  186,
                            5,    4,
                           13,   10,
                           13,   19,
                            8,    3,
                           39,   40
                      ), .Dim = c(2,10))),
     n  = t(structure(.Data = c(
                            65,   61,
                            43,   43,
                            59,   56,
                            32,   34,
                          1212, 1204,
                            50,  100,
                            34,   68,
                            41,   40,
                            40,   40,
                           381,  372
                      ), .Dim = c(2,10))),
     Nsafe = 3,
     Nrisky = 7,
     safe = c(3,5,8),
     risky = c(1:2,4,6:7,9:10),
     bias.mu = -0.2,
     bias.prec = 1 / (0.1^2)
     ))


## Inits
initsSepBias <- list(
            list(
                d=c(NA,0),
                mu=rep(0, 10),
                bias = c(0.5,0.5,NA,0.5,NA,0.5,0.5,NA,0.5,0.5)
                ),
  
            list(
                d=c(NA,0.1),
                mu=c(1,-1,-2,0,0,-2,1,0,2,2),
                bias = c(-0.5,-0.5,NA,-0.5,NA,-0.5,-0.5,NA,-0.5,-0.5)
                )
            )




## Parameters to monitor
## Key parameters N (number at each severity level); c (conditional probabilities)
## Plus nodes monitored to use in calculating deviance summaries
paramsSepBias <- c("p","mu","delta","d","bias","devp","yhatp")


## Initialise model
sepBias.jm <- jags.model(textConnection(sepBias_model),
                       data = sepBias_dat,
                       inits = initsSepBias,
                       n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(sepBias.jm, n.iter = 1000)
## samples to keep
outSepBias <- coda.samples(sepBias.jm,
               variable.names = paramsSepBias,
               n.iter = 9000,
               n.thin = 1)
print(proc.time() - ptm)



## summary stats

## extract treatment effects
## the mean of the random effects, odds ratio scale


## Deviance summaries
round(devsSepBias <- getDevs(post = outSepBias,
                             daty = list(y = sepBias_dat$y),
                             datn = list(n = sepBias_dat$n),
                             pNames = c("^p\\["),
                             likTypes = c("binomial"),
                             devstem = "dev",
                             yhatstem = "yhat",
                             plugin = "mean"), digits = 2)





##################################
## SEPSIS META-ANALYSIS MODEL   ##
## QUESTION 3 - additive biases ##
##################################


## Add means on log-scale:

## Add variances on log-scale:


## Change prior values to reflect the added means/variances



## Initialise model
sepBiasAdd.jm <- jags.model(textConnection(sepBias_model),
                       data = sepBiasAdd_dat,
                       inits = initsSepBias,
                       n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(sepBiasAdd.jm, n.iter = 1000)
## samples to keep
outSepBiasAdd <- coda.samples(sepBiasAdd.jm,
               variable.names = paramsSepBias,
               n.iter = 9000,
               n.thin = 1)
print(proc.time() - ptm)



## summary stats

## Check traces

## extract treatment effects
## the mean of the random effects, odds ratio scale
