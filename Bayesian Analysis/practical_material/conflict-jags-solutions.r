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
      pval.cvmix[i,k] <- step(yrep[i,k] - ycv[i,i,k] - 0.5) + (0.5 * equals(yrep[i,k], ycv[i,i,k]))
    }
    
    ## for each non-baseline arm
    for(k in 2:Na)
    {
      ## cross-validatory mixed-predictive p-value
      ## one-sided, testing logit(delta) > logit(delta.loo)
      ## (not integer so not a mid-p-value)
      pval.cvmix.delta[i,k] <- step(delta[i,i,k] - delta.loo[i,k])
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
dev.new()
chkTraces(outSep, "^delta\\[.,.,2\\]")


## summary stats
(statsSep <- summary(outSep))

## cross-validatory mixed-predictive p-values, comparing yrep to ycv:
pval.cvmix <- getNode(outSep, "pval.cvmix\\[")
apply(pval.cvmix, 2, mean)

## compared to cross-validatory mixed-predictive p-values, comparing delta to 
## T(ycv[i,]) = logit(ycv[i,2]/ncv[i,2]) - logit(ycv[i,1]/ncv[i,1])
pval.cvmix.delta <- getNode(outSep, "pval.cvmix.delta")
apply(pval.cvmix.delta, 2, mean)

## plot comparisons for Sandberg and Fanaroff studies
delta <- getNode(outSep, "^delta\\[")
delta.loo <- logit(sep_dat$ycv[,,2] / sep_dat$ncv[,,2]) - logit(sep_dat$ycv[,,1] / sep_dat$ncv[,,1])

study = c("Sandberg vs rest", "Fanaroff vs rest")
studyNum = c(8,5)

par(mfrow = c(1,2), mar = c(5.7,4.5,2,0.2))
for(k in 1:2)
{
  ti <- density(delta[,grep(paste("\\[",studyNum[k],",",studyNum[k],",2\\]",sep = ""),colnames(delta))], bw = 0.5)
  tl <- delta.loo[studyNum[k]]
  
  plot(ti$x, ti$y, type = "l", xlim = c(-3,3), ylim = c(0,1), xlab = "", ylab = "", las = 1, lwd = 3, cex.axis = 2, cex.lab = 2, main = study[k], cex.main = 2)
  
  abline(v = tl, col = "red")
  pm1 <- round(mean(pval.cvmix[,grep(paste(studyNum[k], ",1", sep=""),colnames(pval.cvmix))]),3)
  pm2 <- round(mean(pval.cvmix[,grep(paste(studyNum[k], ",2", sep=""),colnames(pval.cvmix))]),3)
  pc <- round(mean(pval.cvmix.delta[,grep(paste(studyNum[k], ",2", sep=""),colnames(pval.cvmix.delta))]),3)
  legend("topleft",
         c(as.expression(bquote(p[c] == .(pc))), as.expression(bquote(p[m1] == .(pm1))), as.expression(bquote(p[m2] == .(pm2))), "Left-out", "Rest"),
         col = c("transparent","transparent","transparent","red","black"),
         lty = c(rep(NULL,3),1,1), bty = "n", cex = 1.5)
}


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
(statsSepBias <- summary(outSepBias))

## extract treatment effects
## the mean of the random effects, odds ratio scale
dRE <- exp(getNode(outSepBias, "^d\\["))
dRE <- as.matrix(dRE, dim = c(length(dRE),1))
sumStats(dRE)
##          d[1]      d[2]
## Mean        1 0.7039143
## SD          0 0.2026137
## Median      1 0.6978790
## 2.5%ile     1 0.3378744
## 97.5%ile    1 1.1135639


## Deviance summaries
round(devsSepBias <- getDevs(post = outSepBias,
                             daty = list(y = sepBias_dat$y),
                             datn = list(n = sepBias_dat$n),
                             pNames = c("^p\\["),
                             likTypes = c("binomial"),
                             devstem = "dev",
                             yhatstem = "yhat",
                             plugin = "mean"), digits = 2)
##              y    n pobs yexpbar thetabar  Dbar Dhat    pD   DIC
## devp[1,1]   23   65 0.35   23.97     0.37  0.91 0.06  0.85  1.76
## devp[1,2]   20   61 0.33   19.01     0.31  0.93 0.07  0.86  1.79
## devp[2,1]    8   43 0.19    6.84     0.16  1.20 0.22  0.97  2.17
## devp[2,2]    2   43 0.05    3.12     0.07  1.00 0.49  0.51  1.51
## devp[3,1]    5   59 0.08    3.50     0.06  1.77 0.61  1.17  2.94
## devp[3,2]    0   56 0.00    1.48     0.03  3.03 3.01  0.02  3.04
## devp[4,1]   14   32 0.44   13.22     0.41  0.93 0.08  0.86  1.79
## devp[4,2]    8   34 0.24    8.68     0.26  0.84 0.07  0.77  1.61
## devp[5,1]  209 1212 0.17  209.66     0.17  0.97 0.00  0.97  1.94
## devp[5,2]  186 1204 0.15  185.30     0.15  0.99 0.00  0.99  1.98
## devp[6,1]    5   50 0.10    4.46     0.09  0.89 0.07  0.82  1.71
## devp[6,2]    4  100 0.04    4.53     0.05  0.76 0.07  0.69  1.44
## devp[7,1]   13   34 0.38   11.56     0.34  1.25 0.27  0.98  2.23
## devp[7,2]   10   68 0.15   11.42     0.17  1.08 0.22  0.87  1.95
## devp[8,1]   13   41 0.32   15.07     0.37  1.30 0.46  0.84  2.14
## devp[8,2]   19   40 0.48   16.96     0.42  1.35 0.42  0.93  2.28
## devp[9,1]    8   40 0.20    7.25     0.18  1.01 0.09  0.91  1.92
## devp[9,2]    3   40 0.08    3.74     0.09  0.80 0.17  0.62  1.42
## devp[10,1]  39  381 0.10   40.96     0.11  1.05 0.11  0.94  1.99
## devp[10,2]  40  372 0.11   37.97     0.10  1.13 0.12  1.01  2.15
## TOTAL       NA   NA   NA      NA       NA 23.19 6.62 16.57 39.76





##################################
## SEPSIS META-ANALYSIS MODEL   ##
## QUESTION 3 - additive biases ##
##################################


## Add means on log-scale:
(nu <- log(0.88) + log(0.82) + log(0.77) + log(1.00))
#[1] -0.5876491
exp(nu)
#[1] 0.555632

## Add variances on log-scale:
(tausq <- 0.07^2 + 0.075^2 + 0.11^2 + 0.05^2)
(tau <- tausq^(0.5))
#[1] 0.1585087


## Change prior values to reflect the added means/variances
sepBiasAdd_dat <- sepBias_dat
sepBiasAdd_dat$bias.mu <- nu
sepBiasAdd_dat$bias.prec <- 1 / tausq




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
(statsSepBiasAdd <- summary(outSepBiasAdd))

## Check traces
chkTraces(outSepBiasAdd, "bias")

## extract treatment effects
## the mean of the random effects, odds ratio scale
dREAdd <- exp(getNode(outSepBiasAdd, "^d\\["))
dREAdd <- as.matrix(dREAdd, dim = c(length(dREAdd),1))
sumStats(dREAdd)
##          d[1]      d[2]
## Mean        1 0.9031052
## SD          0 0.2519055
## Median      1 0.8963890
## 2.5%ile     1 0.4408431
## 97.5%ile    1 1.4357678
