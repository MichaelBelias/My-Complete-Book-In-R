library(R2OpenBUGS)
library(coda)

## Set of functions for processing posterior samples in R
## including functions for calculating deviance summaries
source("./fns.R")


#############################################################
## SEPSIS META-ANALYSIS MODEL                              ##
## QUESTION 1 - cross-validatory mixed-predictive p-values ##
#############################################################

## Leave-One-Out cross-validation for Sepsis model
sepCV_model <- function()
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
        ## copy data Ns times, once for each cross-validation
        ycv[i,j,k] <- y[j,k]
        ncv[i,j,k] <- n[j,k]
        
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


# write the model code out to a file
write.model(sepCV_model, con = "sepCV_model.txt")


## data list
(sep_dat <- list(
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
  ), .Dim = c(2,10)))
))


## initial values
initsSep <- list(
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


## Run in R2OpenBUGS
print(ptm <- proc.time())
outSep <- as.mcmc.list(
  bugs(model = "sepCV_model.txt",
       data = sep_dat,
       inits = initsSep,
       n.chains = 2,
       parameters.to.save = c("p","mu","delta","d","pval.cvmix","pval.cvmix.delta"),
       n.burnin = 1000,
       n.iter = 10000,
       n.thin = 3,
       DIC = FALSE))
print(proc.time() - ptm)

## check traces
chkTraces(outSep, "delta")

## summary stats
(statsSep <- summary(outSep))

## cross-validatory mixed-predictive p-values, comparing yrep to ycv:
pval.cvmix <- getNode(outSep, "pval.cvmix\\[")
apply(pval.cvmix, 2, mean)
#pval.cvmix[10,1] pval.cvmix[10,2]  pval.cvmix[1,1]  pval.cvmix[1,2] 
#       0.5924167        0.5673333        0.5250000        0.5006667 
# pval.cvmix[2,1]  pval.cvmix[2,2]  pval.cvmix[3,1]  pval.cvmix[3,2] 
#       0.5566667        0.6005000        0.5986667        0.8335833 
# pval.cvmix[4,1]  pval.cvmix[4,2]  pval.cvmix[5,1]  pval.cvmix[5,2] 
#       0.4990000        0.5170833        0.5589167        0.5365000 
# pval.cvmix[6,1]  pval.cvmix[6,2]  pval.cvmix[7,1]  pval.cvmix[7,2] 
#       0.5867500        0.6059167        0.5100000        0.5475000 
# pval.cvmix[8,1]  pval.cvmix[8,2]  pval.cvmix[9,1]  pval.cvmix[9,2] 
#       0.5219167        0.4685000        0.5600000        0.5830833 

## Compared to cross-validatory mixed-predictive p-values, comparing delta to 
## T(ycv[i,]) = logit(ycv[i,2]/ncv[i,2]) - logit(ycv[i,1]/ncv[i,1])
pval.cvmix.delta <- getNode(outSep, "pval.cvmix.delta")
apply(pval.cvmix.delta, 2, mean)
#pval.cvmix.delta[10,2]  pval.cvmix.delta[1,2]  pval.cvmix.delta[2,2] 
#            0.17633333             0.26616667             0.90866667 
# pval.cvmix.delta[3,2]  pval.cvmix.delta[4,2]  pval.cvmix.delta[5,2] 
#            1.00000000             0.70100000             0.24766667 
# pval.cvmix.delta[6,2]  pval.cvmix.delta[7,2]  pval.cvmix.delta[8,2] 
#            0.70983333             0.86183333             0.03933333 
# pval.cvmix.delta[9,2] 
#            0.78600000 

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


####################################################################################################



############################################
## SEPSIS META-ANALYSIS MODEL             ##
## QUESTION 2 - systematic bias modelling ##
############################################



## RANDOM EFFECTS + BIAS MODEL
sepBias_model <- function()
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


# write the model code out to a file
write.model(sepBias_model, con = "sepBias_model.txt")




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



## Run in R2OpenBUGS
print(ptm <- proc.time())
outSepBias <- as.mcmc.list(
    bugs(model = "sepBias_model.txt",
               data = sepBias_dat,
               inits = initsSepBias,
               n.chains = 2,
               parameters.to.save = c("p","mu","delta","d","bias","devp","yhatp"),
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = FALSE))
print(proc.time() - ptm)


## summary stats
(statsSepBias <- summary(outSepBias))

## extract treatment effects
## the mean of the random effects, odds ratio scale
dRE <- exp(getNode(outSepBias, "^d\\["))
dRE <- as.matrix(dRE, dim = c(length(dRE),1))
sumStats(dRE)
##               [,1]
## Mean     0.6967664
## SD       0.2113277
## Median   0.6906998
## 2.5%ile  0.3194914
## 97.5%ile 1.1066307




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
## devp[1,1]   23   65 0.35   23.96     0.37  0.90 0.06  0.84  1.75
## devp[1,2]   20   61 0.33   19.02     0.31  0.95 0.07  0.88  1.82
## devp[2,1]    8   43 0.19    6.94     0.16  1.16 0.18  0.97  2.13
## devp[2,2]    2   43 0.05    3.09     0.07  0.99 0.47  0.52  1.51
## devp[3,1]    5   59 0.08    3.56     0.06  1.75 0.56  1.19  2.93
## devp[3,2]    0   56 0.00    1.47     0.03  2.99 2.97  0.02  3.00
## devp[4,1]   14   32 0.44   13.32     0.42  0.92 0.06  0.86  1.78
## devp[4,2]    8   34 0.24    8.66     0.25  0.85 0.07  0.78  1.63
## devp[5,1]  209 1212 0.17  209.68     0.17  0.99 0.00  0.99  1.97
## devp[5,2]  186 1204 0.15  185.18     0.15  0.98 0.00  0.98  1.96
## devp[6,1]    5   50 0.10    4.50     0.09  0.88 0.06  0.82  1.71
## devp[6,2]    4  100 0.04    4.52     0.05  0.76 0.06  0.69  1.45
## devp[7,1]   13   34 0.38   11.61     0.34  1.22 0.25  0.97  2.20
## devp[7,2]   10   68 0.15   11.38     0.17  1.08 0.21  0.87  1.95
## devp[8,1]   13   41 0.32   15.02     0.37  1.30 0.44  0.86  2.16
## devp[8,2]   19   40 0.48   16.96     0.42  1.33 0.42  0.91  2.24
## devp[9,1]    8   40 0.20    7.30     0.18  0.98 0.08  0.89  1.87
## devp[9,2]    3   40 0.08    3.73     0.09  0.77 0.17  0.60  1.37
## devp[10,1]  39  381 0.10   41.11     0.11  1.07 0.12  0.95  2.02
## devp[10,2]  40  372 0.11   37.96     0.10  1.15 0.12  1.03  2.17
## TOTAL       NA   NA   NA      NA       NA 23.01 6.39 16.62 39.63






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


## Run in R2OpenBUGS
print(ptm <- proc.time())
outSepBiasAdd <- as.mcmc.list(
    bugs(model = "sepBias_model.txt",
               data = sepBiasAdd_dat,
               inits = initsSepBias,
               n.chains = 2,
               parameters.to.save = c("p","mu","delta","d","bias","devp","yhatp"),
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = FALSE))
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
##               [,1]
## Mean     0.8994474
## SD       0.2616204
## Median   0.8902971
## 2.5%ile  0.4424114
## 97.5%ile 1.4493316


