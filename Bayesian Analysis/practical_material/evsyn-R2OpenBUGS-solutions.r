library(R2OpenBUGS)
library(coda)

## Set of functions for processing posterior samples in R
## including functions for calculating deviance summaries
source("./fns.R")

#############################################
## QUESTION 1 - HIV model                  ##
##            - explore evidence synthesis ##
#############################################

## Model A: using single datum (y1), flat priors
hivModelA <- function()
{
  ## Flat priors, i.e. each a, b = 1 (Unif(0,1))
  ## Set prior parameters in data list
  pi ~ dbeta(a.pi,b.pi)
  delta ~ dbeta(a.delta,b.delta)
  
  ## Likelihoods
  ## prevalence data only
  for(i in 1:1)
  {
    y[i] ~ dbin(p[i], n[i])
  }
  
  ## Proportions in terms of basic and functional parameters
  p[1] <- pi
  p[2] <- pi * (1 - delta)
  p[3] <- delta
}



## Write model code to a text file
write.model(hivModelA, con = "hivModelA.txt")


## DATA - datum implies maximum likelihood estimate of 0.05 for pi
hivDatA <- list(
  y = c(  5, NA),
  n = c(100, NA),
  a.pi = 1, b.pi = 1,
  a.delta = 1, b.delta = 1
)

## INITIAL VALUES FOR 2 CHAINS
hivInits <- list(
  ## chain 1
  list(
    pi = 0.1,
    delta = 0.9
  ),
  ## chain 2
  list(
    pi = 0.2,
    delta = 0.2
  )
)
  

## Parameters to monitor
## Basic & functional parameters p
hivParams <- c("p")



## Run in R2OpenBUGS
print(ptm <- proc.time())
outHIVA1 <- as.mcmc.list(bugs(model = "hivModelA.txt",
               data = hivDatA,
               inits = hivInits,
               n.chains = 2,
               parameters.to.save = hivParams,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE))
print(proc.time() - ptm)


## Check traces for convergence (function provided in fns.R) - the string in the second 
## argument is one that can be passed to the grep command to extract particular 
## parameters from the mcmc.list
## Note how delta (p[3]) is unidentified, i.e. posterior follows the prior,
## whereas pi (p[1]) is identified by the data and pi(1-delta) is partially identified
## (it has pi as its upper bound, since delta lies in [0,1] a priori)
chkTraces(outHIVA1, node = "p")

## Summary statistics - using rjags provided command to look at summaries of all
## monitored parameters
summary(outHIVA1)
## Iterations = 1001:9000
## Thinning interval = 1 
## Number of chains = 2 
## Sample size per chain = 8000 

## 1. Empirical mean and standard deviation for each variable,
##    plus standard error of the mean:

##             Mean      SD  Naive SE Time-series SE
## deviance 4.43209 1.43866 0.0113736      0.0113739
## p[1]     0.05868 0.02305 0.0001823      0.0001781
## p[2]     0.02913 0.02136 0.0001689      0.0001709
## p[3]     0.50269 0.28865 0.0022820      0.0022820

## 2. Quantiles for each variable:

##              2.5%     25%     50%     75%   97.5%
## deviance 3.430000 3.52900 3.87800 4.74600 8.49002
## p[1]     0.022229 0.04201 0.05580 0.07225 0.11180
## p[2]     0.001148 0.01225 0.02528 0.04160 0.08039
## p[3]     0.024119 0.25490 0.50465 0.75075 0.97660


## Summary statistics for a particular parameter only, using functions provided
## in fns.R - second argument of getNode() function is a regular expression that can 
## be passed to the grep command, e.g. to summarise pi = p[1]:
sumStats(as.matrix(getNode(outHIVA1, node = "p\\[1\\]"), ncol = 1))
##                [,1]
## Mean     0.05868274
## SD       0.02305398
## Median   0.05580000
## 2.5%ile  0.02222900
## 97.5%ile 0.11180000



## Add in informative prior for delta (mean 0.75) and re-run
hivDatA$a.delta <- 75
hivDatA$b.delta <- 25


## Run in R2OpenBUGS
print(ptm <- proc.time())
outHIVA2 <- as.mcmc.list(bugs(model = "hivModelA.txt",
               data = hivDatA,
               inits = hivInits,
               n.chains = 2,
               parameters.to.save = hivParams,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE))
print(proc.time() - ptm)


## Summary statistics - note that all parameters (basic + functional) are now identified
sumStats(as.matrix(getNode(outHIVA2, node = "p"), ncol = 1))
##                p[1]        p[2]       p[3]
## Mean     0.05885497 0.014718560 0.75010570
## SD       0.02341074 0.006475577 0.04265271
## Median   0.05581000 0.013705000 0.75150000
## 2.5%ile  0.02185975 0.005064975 0.66170000
## 97.5%ile 0.11230000 0.030060250 0.82920250


## Add in informative prior for pi (mean 0.15) and re-run
hivDatA$a.pi <- 15
hivDatA$b.pi <- 85



## Run in R2OpenBUGS
print(ptm <- proc.time())
outHIVA3 <- as.mcmc.list(bugs(model = "hivModelA.txt",
               data = hivDatA,
               inits = hivInits,
               n.chains = 2,
               parameters.to.save = hivParams,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE))
print(proc.time() - ptm)


## Summary statistics - the informative prior for pi implies a different region of 
## support to the data, so note the posterior for pi = p[1] is a compromise between
## prior and likelihood
sumStats(as.matrix(getNode(outHIVA3, node = "p"), ncol = 1))
##                p[1]        p[2]       p[3]
## Mean     0.09978753 0.024945907 0.75006029
## SD       0.02089474 0.006844129 0.04322352
## Median   0.09864000 0.024255000 0.75140000
## 2.5%ile  0.06254975 0.013520000 0.66050000
## 97.5%ile 0.14380250 0.040210000 0.83040250


## Plot posteriors from the three models against each other:
psA <- list(getNode(outHIVA1,"p"), getNode(outHIVA2,"p"), getNode(outHIVA3,"p"))
psA.dens <- lapply(psA, function(x) { apply(x, 2, density, from = 0, to = 1) })

## Informative delta prior (Model A2) vs flat priors (Model A1)
par(mfrow = c(1,3))
txtMain <- c(expression(pi), expression(pi(1-delta)), expression(delta))
xlims <- list(c(0,0.25), c(0,0.25), c(0,1))
txtXlab = c("Prevalence","Undiagnosed prevalence","Proportion diagnosed")
for(i in 1:3)
{
  plot(psA.dens[[1]][[i]]$x, psA.dens[[1]][[i]]$y, type = "l", main = txtMain[i], xlim = xlims[[i]], xlab = txtXlab[i], ylim = c(0,85), ylab = "", lwd = 3, las = 1, cex = 2, cex.main = 2, cex.lab = 2, cex.axis = 2)
  for(m in 2:2)
    lines(psA.dens[[m]][[i]]$x, psA.dens[[m]][[i]]$y, col = m, lwd = 3)
  legend("topright", c("Flat priors", expression(delta %~% Beta(75,25))), bty = "n", lty = rep(1,2), col = 1:2, cex = 2, lwd = rep(3,2))
}

## Info pi & delta priors (Model A3) vs Informative delta prior (Model A2) vs flat priors (Model 1)
par(mfrow = c(1,3))
for(i in 1:3)
{
  plot(psA.dens[[1]][[i]]$x, psA.dens[[1]][[i]]$y, type = "l", main = txtMain[i], xlim = xlims[[i]], xlab = txtXlab[i], ylim = c(0,85), ylab = "", lwd = 3, las = 1, cex = 2, cex.main = 2, cex.lab = 2, cex.axis = 2)
  for(m in 2:3)
    lines(psA.dens[[m]][[i]]$x, psA.dens[[m]][[i]]$y, col = m, lwd = 3)
  legend("topright", c("Flat priors", expression(delta %~% Beta(75,25)), expression(+ pi %~% Beta(15,85))), bty = "n", lty = rep(1,3), col = 1:3, cex = 2, lwd = rep(3,3))
}




## Model B: using two data points (y1, y2), flat priors
hivModelB <- function()
{
  ## Flat priors, i.e. each a, b = 1 (Unif(0,1))
  ## Set prior parameters in data list
  pi ~ dbeta(a.pi,b.pi)
  delta ~ dbeta(a.delta,b.delta)
  
  ## Likelihoods
  ## prevalence and undiagnosed prevalence data
  for(i in 1:2)
  {
  y[i] ~ dbin(p[i], n[i])
  }
  
  ## Proportions in terms of basic and functional parameters
  p[1] <- pi
  p[2] <- pi * (1 - delta)
  p[3] <- delta
}

## Write model code to a text file
write.model(hivModelB, con = "hivModelB.txt")

## DATA
hivDatB <- list(
  y = c(  5,   3000),
  n = c(100, 100000),
  a.pi = 1, b.pi = 1,
  a.delta = 1, b.delta = 1
)


## Run in R2OpenBUGS
print(ptm <- proc.time())
outHIVB1 <- as.mcmc.list(bugs(model = "hivModelB.txt",
               data = hivDatB,
               inits = hivInits,
               n.chains = 2,
               parameters.to.save = hivParams,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE))
print(proc.time() - ptm)

## Summary statistics - all parameters are now identifiable even with flat priors,
## although the information on delta = p[3] is quite uncertain, since it is indirect
sumStats(as.matrix(getNode(outHIVB1, node = "p"), ncol = 1))
##                p[1]         p[2]       p[3]
## Mean     0.05003542 0.0300013200 0.33770835
## SD       0.01776491 0.0005397389 0.18899735
## Median   0.04479000 0.0300000000 0.33000000
## 2.5%ile  0.03067000 0.0289600000 0.02182975
## 97.5%ile 0.09799050 0.0310800000 0.69460250


## Check traces - note that pi = p[1] and delta = p[3] are now not mixing well. This is because the introduction
## of the data y2 informing pi(1-delta) induces correlation between pi and delta.
chkTraces(outHIVB1, node = "p")


## Scatter plot of pi against delta, to see correlation
par(mfrow = c(1,1))
plot(getNode(outHIVB1, "p\\[1\\]"), getNode(outHIVB1, "p\\[3\\]"), pch = ".", xlab = "pi", ylab = "delta")



## Plot posteriors from A1 vs B1 (both flat priors)
psB <- list(getNode(outHIVA1,"p"), getNode(outHIVB1,"p"))
psB.dens <- lapply(psB, function(x) { apply(x, 2, density, from = 0, to = 1) })

par(mfrow = c(1,3))
for(i in 1:3)
{
  plot(psB.dens[[1]][[i]]$x, psB.dens[[1]][[i]]$y, type = "l", main = txtMain[i], xlim = xlims[[i]], xlab = txtXlab[i], ylim = c(0,85), ylab = "", lwd = 3, las = 1, cex = 2, cex.main = 2, cex.lab = 2, cex.axis = 2)
  lines(psB.dens[[2]][[i]]$x, psB.dens[[2]][[i]]$y, col = 2, lwd = 3)
  legend("topright", expression(list(y[1] == 5, n[1] == 100), list(y[2] == 3000, n[2] == 100000)), bty = "n", lty = rep(1,2), col = 1:2, cex = 1.5, lwd = rep(3,2))
}




## Model C: using 3 data points (y1, y2, y3), flat priors
hivModelC <- function()
{
  ## Flat priors, i.e. each a, b = 1 (Unif(0,1))
  ## Set prior parameters in data list
  pi ~ dbeta(a.pi,b.pi)
  delta ~ dbeta(a.delta,b.delta)
  
  ## Likelihoods
  ## prevalence, undiagnosed prevalence and prop'n diagnosed data
  for(i in 1:3)
  {
  y[i] ~ dbin(p[i], n[i])
  }
  
  ## Proportions in terms of basic and functional parameters
  p[1] <- pi
  p[2] <- pi * (1 - delta)
  p[3] <- delta
}

## Write model code to a text file
write.model(hivModelC, con = "hivModelC.txt")


## DATA
hivDatC <- list(
  y = c(  5,   3000,  90),
  n = c(100, 100000, 100),
  a.pi = 1, b.pi = 1,
  a.delta = 1, b.delta = 1
)



## Run in R2OpenBUGS
print(ptm <- proc.time())
outHIVC1 <- as.mcmc.list(bugs(model = "hivModelC.txt",
               data = hivDatC,
               inits = hivInits,
               n.chains = 2,
               parameters.to.save = hivParams,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE))
print(proc.time() - ptm)

## Summary statistics - all parameters are now identifiable even with flat priors,
## although the information on delta = p[3] is quite uncertain, since it is indirect
sumStats(as.matrix(getNode(outHIVC1, node = "p"), ncol = 1))
##                p[1]         p[2]       p[3]
## Mean     0.15547792 0.0298971869 0.80338561
## SD       0.02345987 0.0005287353 0.02948302
## Median   0.15400000 0.0299000000 0.80540000
## 2.5%ile  0.11469750 0.0288800000 0.73960000
## 97.5%ile 0.20460000 0.0309800000 0.85510000


## Check traces - the introduction of the third datum on delta allows the mixing to improve again, despite the correlation between pi and delta
chkTraces(outHIVC1, node = "p")


## Scatter plot of pi against delta, to see correlation
par(mfrow = c(1,1))
plot(getNode(outHIVC1, "p\\[1\\]"), getNode(outHIVC1, "p\\[3\\]"), pch = ".", xlab = "pi", ylab = "delta")


## Plot posteriors from A1 vs B1 vs C1 (all flat priors)
psC <- list(getNode(outHIVA1,"p"), getNode(outHIVB1,"p"), getNode(outHIVC1,"p"))
psC.dens <- lapply(psC, function(x) { apply(x, 2, density, from = 0, to = 1) })

par(mfrow = c(1,3))
for(i in 1:3)
{
  plot(psC.dens[[1]][[i]]$x, psC.dens[[1]][[i]]$y, type = "l", main = txtMain[i], xlim = xlims[[i]], xlab = txtXlab[i], ylim = c(0,85), ylab = "", lwd = 3, las = 1, cex = 2, cex.main = 2, cex.lab = 2, cex.axis = 2)
  for(m in 2:3)
    lines(psC.dens[[m]][[i]]$x, psC.dens[[m]][[i]]$y, col = m, lwd = 3)
  legend("topright", expression(list(y[1] == 5, n[1] == 100), list(y[2] == 3000, n[2] == 100000), list(y[3] == 90, n[3] == 100)), bty = "n", lty = rep(1,3), col = 1:3, cex = 1.5, lwd = rep(3,3))
}




################################################################################################

###################################################
## QUESTION 2 - 'Flu severity model              ##
##            - exploration & deviance summaries ##
###################################################

## Model code
fluModel <- function()
{
    ## For each age group
    for(a in 1:A)
    {
        ## Number of infections at each severity level l
        ## Population size is a constant
        for(l in INF:INF)
        {
            Ncont[a,l] <- p[a,l] * Npop[a]
            N[a,l] <- round(Ncont[a,l])
        }
        ## Each number of infections is (mean) pr{l|l-1}*number at level below
        ## SYM, GP, HOS, ICU, DEA
        for(l in SYM:DEA)
        {
            Ncont[a,l] <- Ncont[a,l-1] * p[a,l-1]
            N[a,l] <- round(Ncont[a,l])
        }

        ## Four conditional probabilities with flat priors
        ## (I|P, H|G, ICU|H, D|ICU)
        ## Two conditional probabilities with informative priors
        ## S|I, G|S
        for(l in INF:DEA)
        {
            p[a,l] ~ dbeta(pA[l], pB[l])
        }

        ## Observed prevalences indirectly informing incidence
        ## pre- and post-pandemic
        for(t in 1:2)
        {
            ypi[a,t] ~ dbin(pi[a,t], npi[a,t])

            ## for model criticism - on original scale
            yhatpi[a,t] <- pi[a,t] * npi[a,t]
            devpi[a,t] <- 2 * (ypi[a,t] * (log(ypi[a,t])-log(yhatpi[a,t])) + (npi[a,t]-ypi[a,t]) * (log(npi[a,t]-ypi[a,t]) - log(npi[a,t]-yhatpi[a,t])))

            ## on logit scale
            lgtpi[a,t] <- logit(pi[a,t])
            devlgtpi[a,t] <- -2 * ((ypi[a,t] * (lgtpi[a,t] - log(ypi[a,t]/npi[a,t]) + log((npi[a,t]-ypi[a,t])/npi[a,t]))) - (npi[a,t] * (log(1 + exp(lgtpi[a,t])) + log((npi[a,t]-ypi[a,t])/npi[a,t]))))
        }

        ## Flat prior for pre-pandemic prevalence,
        ## post-pandemic prevalence defined in terms of pre-pandemic and incidence I|P
        pi[a,1] ~ dbeta(1,1)
        pi[a,2] <- pi[a,1] + p[a,INF]

        ## Observed number at 3 levels, NegBin likelihood, allowing for over-dispersion
        ## G, H, D
        for(l in 1:3)
        {
            y[a,l] ~ dnegbin(psi[a,l],r[l])
            psi[a,l] <- r[l] / (N[a,counts[l]] + r[l])

            ## for model criticism - on original scale
            yhatpsi[a,l] <- r[l] * (1 - psi[a,l]) / psi[a,l]
            devpsi[a,l] <- 2 * (r[l]*(log(r[l]) - log(psi[a,l] * (r[l]+y[a,l]))) + y[a,l]*(log(y[a,l]) - log((1 - psi[a,l])*(r[l]+y[a,l]))))

            ## on logit/log scale
            lgtpsi[a,l] <- logit(psi[a,l])
            devlgtpsi[a,l] <- 2 * ((exp(logr[l])*(logr[l] - lgtpsi[a,l])) - ((exp(logr[l]) + y[a,l]) * (log(exp(logr[l]) + y[a,l]) - log(1 + exp(lgtpsi[a,l])))) + (y[a,l] * log(y[a,l])))
        }

        ## Observations informing 2 conditional probabilities
        ## ICU|HOS, DEA|ICU
        for(l in 1:2)
        {
            yp[a,l] ~ dbin(p[a,probs[l]], np[a,l])

            ## for model criticism - on original scale
            yhatp[a,l] <- p[a,probs[l]] * np[a,l]
            devp[a,l] <- 2 * (yp[a,l] * (log(yp[a,l])-log(yhatp[a,l])) +  (np[a,l]-yp[a,l]) * (log(np[a,l]-yp[a,l]) - log(np[a,l]-yhatp[a,l])))

            ## on logit/log scale
            lgtp[a,l] <- logit(p[a,probs[l]])
            devlgtp[a,l] <- -2 * ((yp[a,l] * (lgtp[a,l] - log(yp[a,l]/np[a,l]) + log((np[a,l]-yp[a,l])/np[a,l]))) - (np[a,l] * (log(1 + exp(lgtp[a,l])) + log((np[a,l]-yp[a,l])/np[a,l]))))
        }
    }

    ## over-dispersion parameter for negative binomial counts of severe outcomes
    for(l in 1:3)
    {
        r[l] <- round(exp(logr[l]))
        logr[l] <- 1 / invlogr[l]
        invlogr[l] ~ dunif(odL,odU)
    }
}

## Write model code to a text file
write.model(fluModel, con = "fluModel.txt")


## Set up data list
INF <- 1; SYM <- 2; GP <- 3; HOS <- 4; ICU <- 5; DEA <- 6;

fluDat <- list(
    A   = 3,  ## Number of age groups
    INF = INF,  ## Severity levels
    SYM = SYM,
    GP  = GP,
    HOS = HOS,
    ICU = ICU,
    DEA = DEA,

    Npop = c(8999200, 33933500, 8159400),  ## Population size by age group
    
    pA   = c(1,40,1,1,1,1),  ## a & b parameters for beta distributions
    pB   = c(1,60,9,1,1,1),  ## ~40% (30-50%) for S|I and ~10% (0-33%) for G|S, flat for rest
             
    counts = c(GP,HOS,DEA),  ## levels at which absolute number of infections observed: G,H,D
    probs  = c(ICU,DEA),    ## conditional probabilities observed binomially: ICU|H, D|ICU

    ## observed number of infections
    y = t(structure(.Data = c(
                        ##  G,   H,  D
                        101039, 1247,  94, ## <15
                        331773, 2256, 576, ## 15-64
                          6982,  210, 144  ## 65+
                        ), .Dim = c(3,3))),

    ## Uniform(odL, odU) prior for over-dispersion parameter r, on inverse log-scale
    odL = 0.05,
    odU = 0.5,

    ## for observed hospitalisation where outcome (ICU, D, discharge) is known
    ## informs conditional probabilities: ICU|H, D|ICU
    yp = t(structure(.Data = c(
                         ## ICU,   D
                             41,  41, ## <15
                            350, 374, ## 15-64
                             22,  13  ## 65+
                         ), .Dim = c(2,3))),
    np = t(structure(.Data = c(
                         ##   H,  ICU
                            883,  210,  ## <15
                           1722, 1676,  ## 15-64
                            168,   52   ## 65+
                         ), .Dim = c(2,3))),

    ## observed prevalences pre- and post-pandemic
    ypi = t(structure(.Data = c(
                          ## pre, post,
                              10,  236,  ## <15
                              66,  379,  ## 15-64
                             128,  168   ## 65+
                          ), .Dim = c(2,3))),
    npi = t(structure(.Data = c(
                          ## pre, post,
                             359,  423,  ## <15
                             495,  978,  ## 15-64
                             549,  445   ## 65+
                          ), .Dim = c(2,3)))
)


## List of lists of initial values, one list for each chain
fluInits <- list(
    list(
        p = t(structure(.Data = c(
                             0.3, 0.4, 0.1, 0.1, 0.5, 0.5,
                             0.3, 0.4, 0.1, 0.1, 0.5, 0.5,
                             0.3, 0.4, 0.1, 0.1, 0.5, 0.5
                             ), .Dim = c(6,3))),
        pi = t(structure(.Data = c(
                              0.1,NA,
                              0.1,NA,
                              0.1,NA
                          ), .Dim = c(2,3))),
        invlogr = c(0.2,0.4,0.3)
    ),
    
    list(
        
        p = t(structure(.Data = c(
                             0.5, 0.2, 0.2, 0.1, 0.5, 0.5,
                             0.5, 0.2, 0.2, 0.1, 0.5, 0.5,
                             0.5, 0.2, 0.2, 0.1, 0.5, 0.5
                             ), .Dim = c(6,3))),
        pi = t(structure(.Data = c(
                              0.5,NA,
                              0.5,NA,
                              0.5,NA
                          ), .Dim = c(2,3))),
        invlogr = c(0.1,0.1,0.1)
    )
)


## Parameters to monitor
## Key parameters N (number at each severity level); c (conditional probabilities)
## Plus nodes monitored to use in calculating deviance summaries
params <- c("N","p","pi","invlogr","logr","r","psi","lgtpi","lgtpsi","lgtp","yhatpi",
            "devpi","devlgtpi","yhatpsi","devpsi","devlgtpsi","yhatp","devp","devlgtp")


## Run in R2OpenBUGS
print(ptm <- proc.time())
outFlu <- bugs(model = "fluModel.txt",
               data = fluDat,
               inits = fluInits,
               n.chains = 2,
               parameters.to.save = params,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE)
print(proc.time() - ptm)


## pD and DIC based on unstandardised deviance (OpenBUGS default)
outFlu$pD
## [1] 11.24
outFlu$DIC
## [1] 281.1

## mcmc.list object to use in summary/plotting functions
outFluMC <- as.mcmc.list(outFlu)


## Check traces for convergence - the string in the second argument is one that can be passed to the grep command to extract particular parameters from the mcmc.list
chkTraces(outFluMC, "^pi")
chkTraces(outFluMC, "^N")
chkTraces(outFluMC, "^p\\[")
chkTraces(outFluMC, "^r")
chkTraces(outFluMC, "invlogr")
chkTraces(outFluMC, "^psi")

## Summary statistics
summary(outFluMC)


## Calculate deviance summaries, using posterior mean for plug-in
devsFluMean <- getDevs(post = outFluMC,
                       daty = list(fluDat$yp, fluDat$ypi, fluDat$y),
                       datn = list(fluDat$np, fluDat$npi, rep(NA, length(fluDat$y))),
                       pNames = c("^p\\[.,[5|6]\\]","^pi","^psi\\["),
                       likTypes = c("binomial","binomial","negbin"),
                       aux = c(NA, NA, "^r\\["),
                       devstem = "dev[^i|^lgt]",
                       yhatstem = "yhat",
                       plugin = "mean")
round(devsFluMean, digits = 2)
## Notice that, as the posterior distributions of some parameters are skewed, the posterior mean is not a good summary of the posterior, hence the strange large values for the plug-in deviance Dhat and therefore the large negative values of pD and DIC.
##                  y    n pobs   yexpbar thetabar  Dbar        Dhat           pD          DIC
## devp[1,1]       41  883 0.05     44.10     0.05  1.12        0.24         0.89         2.01
## devp[1,2]       41  210 0.20     41.64     0.20  1.00        0.01         0.99         1.99
## devp[2,1]      350 1722 0.20    352.22     0.20  1.00        0.02         0.98         1.98
## devp[2,2]      374 1676 0.22    374.53     0.22  1.01        0.00         1.01         2.01
## devp[3,1]       22  168 0.13     27.48     0.16  2.16        1.39         0.78         2.94
## devp[3,2]       13   52 0.25     13.47     0.26  0.99        0.02         0.96         1.95
## devpi[1,1]      10  359 0.03     11.85     0.03  1.25        0.32         0.94         2.19
## devpi[1,2]     236  423 0.56    228.54     0.54  1.54        0.53         1.01         2.55
## devpi[2,1]      66  495 0.13     72.16     0.15  1.62        0.63         0.99         2.62
## devpi[2,2]     379  978 0.39    369.45     0.38  1.42        0.40         1.02         2.44
## devpi[3,1]     128  549 0.23    146.30     0.27  3.94        3.21         0.73         4.66
## devpi[3,2]     168  445 0.38    147.32     0.33  4.89        4.25         0.65         5.54
## devpsi[1,1] 101039   NA   NA 808959.18     0.00 18.82       18.57         0.25        19.07
## devpsi[1,2]   1247   NA   NA   2247.87     0.00  2.53        1.79         0.74         3.27
## devpsi[1,3]     94   NA   NA     84.53     0.33  1.09  3954266.98  -3954265.89  -3954264.80
## devpsi[2,1] 331773   NA   NA 680025.48     0.00  3.23        2.85         0.38         3.62
## devpsi[2,2]   2256   NA   NA   3590.11     0.00  1.75        1.11         0.64         2.39
## devpsi[2,3]    576   NA   NA    544.61     0.17  0.92  6268936.72  -6268935.80  -6268934.87
## devpsi[3,1]   6982   NA   NA  13518.17     0.00  2.91        1.38         1.53         4.43
## devpsi[3,2]    210   NA   NA    803.18     0.01  9.73        8.80         0.93        10.65
## devpsi[3,3]    144   NA   NA    110.48     0.30  2.41  4329859.98  -4329857.56  -4329855.15
## TOTAL           NA   NA   NA        NA       NA 65.34 14553109.18 -14553043.85 -14552978.51



## Calculate deviance summaries, using posterior median for plug-in
devsFluMedian <- getDevs(post = outFluMC,
                         daty = list(fluDat$yp, fluDat$ypi, fluDat$y),
                         datn = list(fluDat$np, fluDat$npi, rep(NA, length(fluDat$y))),
                         pNames = c("^p\\[.,[5|6]\\]","^pi","^psi\\["),
                         likTypes = c("binomial","binomial","negbin"),
                         aux = c(NA, NA, "^r\\["),
                         devstem = "dev[^i|^lgt]",
                         yhatstem = "yhat",
                         plugin = "median")
round(devsFluMedian, digits = 2)
## The posterior median is a better summary, hence we obtain a much more reasonable estimate of Dhat. Note however that we still obtain negative values of the effective number of parameters, pD, which can occur when posterior is skewed or if sampling distribution is not log-concave, or if there is strong prior data conflict (see next lecture).
##                  y    n pobs   yexpbar thetabar  Dbar  Dhat    pD   DIC
## devp[1,1]       41  883 0.05     43.81     0.05  1.12  0.19  0.93  2.05
## devp[1,2]       41  210 0.20     41.43     0.20  1.00  0.01  0.99  1.99
## devp[2,1]      350 1722 0.20    352.10     0.20  1.00  0.02  0.98  1.98
## devp[2,2]      374 1676 0.22    374.50     0.22  1.01  0.00  1.01  2.01
## devp[3,1]       22  168 0.13     27.27     0.16  2.16  1.28  0.88  3.04
## devp[3,2]       13   52 0.25     13.31     0.26  0.99  0.01  0.98  1.96
## devpi[1,1]      10  359 0.03     11.52     0.03  1.25  0.22  1.03  2.29
## devpi[1,2]     236  423 0.56    228.60     0.54  1.54  0.52  1.02  2.56
## devpi[2,1]      66  495 0.13     72.12     0.15  1.62  0.62  1.00  2.63
## devpi[2,2]     379  978 0.39    369.00     0.38  1.42  0.43  0.98  2.40
## devpi[3,1]     128  549 0.23    145.90     0.27  3.94  3.08  0.85  4.79
## devpi[3,2]     168  445 0.38    147.20     0.33  4.89  4.30  0.60  5.49
## devpsi[1,1] 101039   NA   NA 799300.00     0.00 18.82 19.41 -0.59 18.23
## devpsi[1,2]   1247   NA   NA   2113.00     0.00  2.53  1.52  1.01  3.53
## devpsi[1,3]     94   NA   NA     84.00     0.17  1.09  0.60  0.50  1.59
## devpsi[2,1] 331773   NA   NA 665700.00     0.00  3.23  3.27 -0.04  3.20
## devpsi[2,2]   2256   NA   NA   3346.00     0.00  1.75  0.83  0.92  2.67
## devpsi[2,3]    576   NA   NA    548.00     0.03  0.92  0.27  0.65  1.57
## devpsi[3,1]   6982   NA   NA  11790.00     0.00  2.91  1.94  0.97  3.88
## devpsi[3,2]    210   NA   NA    772.00     0.01  9.73  8.43  1.30 11.03
## devpsi[3,3]    144   NA   NA    110.00     0.14  2.41  2.42 -0.01  2.40
## TOTAL           NA   NA   NA        NA       NA 65.34 49.39 15.95 81.29


## Calculate deviance summaries on log/lgt scale, using posterior mean for plug-in
devsFluLgt <- getDevs(post = outFluMC,
                      daty = list(fluDat$yp, fluDat$ypi, fluDat$y),
                      datn = list(fluDat$np, fluDat$npi, rep(NA, length(fluDat$y))),
                      pNames = c("^lgtp\\[.,[1|2]\\]","^lgtpi","^lgtpsi\\["),
                      likTypes = c("binomial.lgt","binomial.lgt","negbin.lgt"),
                      aux = c(NA, NA, "^logr\\["),
                      devstem = "devlgt",
                      yhatstem = "yhat",
                      plugin = "mean")
round(devsFluLgt, digits = 2)
## Calculating the deviance on a scale for the parameters which is more likely to approach Gaussianity (e.g. logit scale for proportions, log scale for counts), helps to obtain values of pD that are positive. Note we still have one negative value for pD.
##                     y    n pobs   yexpbar thetabar  Dbar  Dhat    pD   DIC
## devlgtp[1,1]       41  883 0.05     44.10     0.05  1.12  0.18  0.95  2.07
## devlgtp[1,2]       41  210 0.20     41.64     0.20  1.00  0.00  1.00  2.00
## devlgtp[2,1]      350 1722 0.20    352.22     0.20  1.00  0.01  0.98  1.98
## devlgtp[2,2]      374 1676 0.22    374.53     0.22  1.01  0.00  1.01  2.02
## devlgtp[3,1]       22  168 0.13     27.48     0.16  2.16  1.24  0.92  3.08
## devlgtp[3,2]       13   52 0.25     13.47     0.26  0.99  0.01  0.98  1.97
## devlgtpi[1,1]      10  359 0.03     11.85     0.03  1.25  0.17  1.08  2.33
## devlgtpi[1,2]     236  423 0.56    228.54     0.54  1.54  0.53  1.02  2.56
## devlgtpi[2,1]      66  495 0.13     72.16     0.15  1.62  0.56  1.07  2.69
## devlgtpi[2,2]     379  978 0.39    369.45     0.38  1.42  0.41  1.01  2.43
## devlgtpi[3,1]     128  549 0.23    146.30     0.27  3.94  3.15  0.79  4.72
## devlgtpi[3,2]     168  445 0.38    147.32     0.33  4.89  4.29  0.60  5.49
## devlgtpsi[1,1] 101039   NA   NA 808959.18     0.00 19.32 19.25  0.07 19.39
## devlgtpsi[1,2]   1247   NA   NA   2247.87     0.00  2.59  2.18  0.41  3.00
## devlgtpsi[1,3]     94   NA   NA     84.53     0.33  1.10  0.62  0.47  1.57
## devlgtpsi[2,1] 331773   NA   NA 680025.48     0.00  3.42  3.22  0.20  3.62
## devlgtpsi[2,2]   2256   NA   NA   3590.11     0.00  1.80  1.39  0.41  2.22
## devlgtpsi[2,3]    576   NA   NA    544.61     0.17  0.93  0.32  0.61  1.53
## devlgtpsi[3,1]   6982   NA   NA  13518.17     0.00  3.07  2.15  0.92  3.99
## devlgtpsi[3,2]    210   NA   NA    803.18     0.01  9.86  9.69  0.17 10.02
## devlgtpsi[3,3]    144   NA   NA    110.48     0.30  2.42  3.80 -1.38  1.04
## TOTAL              NA   NA   NA        NA       NA 66.44 53.15 13.29 79.73




##################################################
## QUESTION 3 - 'Flu severity model             ##
##            - Calculate case-severity risks   ##
##################################################



## Calculate case-severity risks from conditional probabilities. Since we already have posterior samples of the conditional probabilities p[a,l], we can calculate the functional parameters directly:
ps <- getNode(outFluMC, "^p\\[")
colnames(ps)

## symptomatic case-severity risks are the probabilities of a severe event conditional on a symptomatic infection (level 2)
sCHR <- ps[,grep(",4",colnames(ps))] * ps[,grep(",3",colnames(ps))]
sCIR <- ps[,grep(",5",colnames(ps))] * sCHR
sCFR <- ps[,grep(",6",colnames(ps))] * sCIR

## case-severity risks are the probabilities of a severe event conditional on an infection (level 1)
## and can be calculated as a symptomatic case-severity risk multiplied by Pr{Symptomatic | Infected}
CHR <-  sCHR * ps[,grep(",2",colnames(ps))]
CIR <-  sCIR * ps[,grep(",2",colnames(ps))]
CFR <-  sCFR * ps[,grep(",2",colnames(ps))]


## look at estimates
sumStats(sCHR)
sumStats(sCIR)
sumStats(sCFR)
sumStats(CHR)
sumStats(CIR)
sumStats(CFR)


## Alternatively, amend the model code to calculate the case-severity risks within OpenBUGS

## Model code
fluModelCSR <- function()
{
    ## For each age group
    for(a in 1:A)
    {
        ## Number of infections at each severity level l
        ## Population size is a constant
        for(l in INF:INF)
        {
            Ncont[a,l] <- p[a,l] * Npop[a]
            N[a,l] <- round(Ncont[a,l])
        }
        ## Each number of infections is (mean) pr{l|l-1}*number at level below
        ## SYM, GP, HOS, ICU, DEA
        for(l in SYM:DEA)
        {
            Ncont[a,l] <- Ncont[a,l-1] * p[a,l-1]
            N[a,l] <- round(Ncont[a,l])
        }

        ## Four conditional probabilities with flat priors
        ## (I|P, H|G, ICU|H, D|ICU)
        ## Two conditional probabilities with informative priors
        ## S|I, G|S
        for(l in INF:DEA)
        {
            p[a,l] ~ dbeta(pA[l], pB[l])
        }

        ## Observed prevalences indirectly informing incidence
        ## pre- and post-pandemic
        for(t in 1:2)
        {
            ypi[a,t] ~ dbin(pi[a,t], npi[a,t])

            ## for model criticism - on original scale
            yhatpi[a,t] <- pi[a,t] * npi[a,t]
            devpi[a,t] <- 2 * (ypi[a,t] * (log(ypi[a,t])-log(yhatpi[a,t])) + (npi[a,t]-ypi[a,t]) * (log(npi[a,t]-ypi[a,t]) - log(npi[a,t]-yhatpi[a,t])))

            ## on logit scale
            lgtpi[a,t] <- logit(pi[a,t])
            devlgtpi[a,t] <- -2 * ((ypi[a,t] * (lgtpi[a,t] - log(ypi[a,t]/npi[a,t]) + log((npi[a,t]-ypi[a,t])/npi[a,t]))) - (npi[a,t] * (log(1 + exp(lgtpi[a,t])) + log((npi[a,t]-ypi[a,t])/npi[a,t]))))
        }

        ## Flat prior for pre-pandemic prevalence,
        ## post-pandemic prevalence defined in terms of pre-pandemic and incidence I|P
        pi[a,1] ~ dbeta(1,1)
        pi[a,2] <- pi[a,1] + p[a,INF]

        ## Observed number at 3 levels, NegBin likelihood, allowing for over-dispersion
        ## G, H, D
        for(l in 1:3)
        {
            y[a,l] ~ dnegbin(psi[a,l],r[l])
            psi[a,l] <- r[l] / (N[a,counts[l]] + r[l])

            ## for model criticism - on original scale
            yhatpsi[a,l] <- r[l] * (1 - psi[a,l]) / psi[a,l]
            devpsi[a,l] <- 2 * (r[l]*(log(r[l]) - log(psi[a,l] * (r[l]+y[a,l]))) + y[a,l]*(log(y[a,l]) - log((1 - psi[a,l])*(r[l]+y[a,l]))))

            ## on logit/log scale
            lgtpsi[a,l] <- logit(psi[a,l])
            devlgtpsi[a,l] <- 2 * ((exp(logr[l])*(logr[l] - lgtpsi[a,l])) - ((exp(logr[l]) + y[a,l]) * (log(exp(logr[l]) + y[a,l]) - log(1 + exp(lgtpsi[a,l])))) + (y[a,l] * log(y[a,l])))
        }

        ## Observations informing 2 conditional probabilities
        ## ICU|HOS, DEA|ICU
        for(l in 1:2)
        {
            yp[a,l] ~ dbin(p[a,probs[l]], np[a,l])

            ## for model criticism - on original scale
            yhatp[a,l] <- p[a,probs[l]] * np[a,l]
            devp[a,l] <- 2 * (yp[a,l] * (log(yp[a,l])-log(yhatp[a,l])) +  (np[a,l]-yp[a,l]) * (log(np[a,l]-yp[a,l]) - log(np[a,l]-yhatp[a,l])))

            ## on logit/log scale
            lgtp[a,l] <- logit(p[a,probs[l]])
            devlgtp[a,l] <- -2 * ((yp[a,l] * (lgtp[a,l] - log(yp[a,l]/np[a,l]) + log((np[a,l]-yp[a,l])/np[a,l]))) - (np[a,l] * (log(1 + exp(lgtp[a,l])) + log((np[a,l]-yp[a,l])/np[a,l]))))
        }
    }

    ## over-dispersion parameter for negative binomial counts of severe outcomes
    for(l in 1:3)
    {
        r[l] <- round(exp(logr[l]))
        logr[l] <- 1 / invlogr[l]
        invlogr[l] ~ dunif(odL,odU)
    }


    ## Calculate functional parameters the case-severity risks
    for(a in 1:A)
    {
        ## Pr(HOS | SYM) = Pr(HOS | GP) * Pr(GP | SYM)
        sCHR[a] <- p[a,HOS] * p[a,GP]
        ## Pr(ICU | SYM) = Pr(ICU | HOS) * Pr(HOS | SYM)
        sCIR[a] <- p[a,ICU] * sCHR[a]
        ## Pr(DEA | SYM) = Pr(DEA | ICU) * Pr(ICU | SYM)
        sCFR[a] <- p[a,DEA] * sCIR[a]

        ## Corresponding infection-severity risks are symptomatic case-severity risks multiplied
        ## by Pr(SYM | INF)
        CHR[a] <- sCHR[a] * p[a,SYM]
        CIR[a] <- sCIR[a] * p[a,SYM]
        CFR[a] <- sCFR[a] * p[a,SYM]
    }
}

## Write model code to a text file
write.model(fluModelCSR, con = "fluModelCSR.txt")


## Monitor case-severity risks
paramsCSR <- c("sCHR","sCIR","sCFR","CHR","CIR","CFR")


## Run in R2OpenBUGS
print(ptm <- proc.time())
outFluCSR <- bugs(model = "fluModelCSR.txt",
               data = fluDat,
               inits = fluInits,
               n.chains = 2,
               parameters.to.save = paramsCSR,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE)
print(proc.time() - ptm)


## mcmc.list
outFluCSRmc <- as.mcmc.list(outFluCSR)

## Check convergence
chkTraces(outFluCSRmc, "C")

## Posterior summaries
CSRs <- getNode(outFluCSRmc, "C")
sumStats(CSRs)





##################################################
## QUESTION 4 - 'Flu severity model             ##
##            - Explore over-dispersion         ##
##################################################



## Look at posterior distribution of the over-dispersion parameter r, on original scale.
odR <- getNode(outFluMC, "^r\\[")
sumStats(odR)
chkTraces(outFluMC, "^r\\[")

## Look at posterior distribution of 1 / log(r) and compare to its Uniform(0.05,0.5) prior.
odinvlogr <- getNode(outFluMC, "invlogr")
sumStats(odinvlogr)
chkTraces(outFluMC, "invlogr")

## Note that a smaller value of r implies more over-dispersion, whereas a larger value of r imples less over-dispersion. So conversely, a smaller value of 1 / log(r) implies less over-dispersion and a larger values implies more. Note that for the first two severity levels (GP, HOS), the posterior of 1 / log(r) is concentrated towards the upper bound (0.5) of its prior distribution, suggesting we should allow more over-dispersion for these levels than the prior allows. By contrast, the DEA level's posterior is more evenly spread over the prior space, although still with a mode towards the upper end of the prior range.



## Allow for more over-dispersion, by increasing the upper bound of the Uniform prior for 1 / log(r)
fluDatMoreOD <- fluDat
fluDatMoreOD$odU <- 3


## Run in R2OpenBUGS
print(ptm <- proc.time())
outFluMoreOD <- bugs(model = "fluModel.txt",
               data = fluDatMoreOD,
               inits = fluInits,
               n.chains = 2,
               parameters.to.save = params,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE)
print(proc.time() - ptm)

outFluMoreODmc <- as.mcmc.list(outFluMoreOD)

## Look at posterior distribution of the over-dispersion parameter r, on original scale.
odRmoreOD <- getNode(outFluMoreODmc, "^r\\[")
sumStats(odRmoreOD)
chkTraces(outFluMoreODmc, "^r\\[")

## Look at posterior distribution of 1 / log(r) and compare to its Uniform(0.05,3) prior.
odinvlogrmoreOD <- getNode(outFluMoreODmc, "invlogr")
sumStats(odinvlogrmoreOD)
chkTraces(outFluMoreODmc, "invlogr")


## Look at deviance summaries for this model on logit/log scale
devsFluMoreODLgt <- getDevs(post = outFluMoreODmc,
                            daty = list(fluDatMoreOD$yp, fluDatMoreOD$ypi, fluDatMoreOD$y),
                            datn = list(fluDatMoreOD$np, fluDatMoreOD$npi, rep(NA, length(fluDatMoreOD$y))),
                            pNames = c("^lgtp\\[.,[1|2]\\]","^lgtpi","^lgtpsi\\["),
                            likTypes = c("binomial.lgt","binomial.lgt","negbin.lgt"),
                            aux = c(NA, NA, "^logr\\["),
                            devstem = "devlgt",
                            yhatstem = "yhat",
                            plugin = "mean")
round(devsFluMoreODLgt, digits = 2)
##                     y    n pobs   yexpbar thetabar  Dbar  Dhat    pD   DIC
## devlgtp[1,1]       41  883 0.05     41.93     0.05  1.00  0.01  0.99  1.99
## devlgtp[1,2]       41  210 0.20     41.54     0.20  0.98  0.00  0.98  1.96
## devlgtp[2,1]      350 1722 0.20    350.80     0.20  1.00  0.00  1.00  1.99
## devlgtp[2,2]      374 1676 0.22    374.68     0.22  1.01  0.00  1.01  2.03
## devlgtp[3,1]       22  168 0.13     23.36     0.14  1.03  0.05  0.98  2.01
## devlgtp[3,2]       13   52 0.25     13.47     0.26  1.00  0.01  1.00  2.00
## devlgtpi[1,1]      10  359 0.03     11.17     0.03  1.07  0.05  1.02  2.10
## devlgtpi[1,2]     236  423 0.56    233.58     0.55  1.01  0.05  0.96  1.98
## devlgtpi[2,1]      66  495 0.13     68.33     0.14  1.09  0.07  1.02  2.11
## devlgtpi[2,2]     379  978 0.39    375.65     0.38  1.04  0.05  0.99  2.03
## devlgtpi[3,1]     128  549 0.23    133.54     0.24  1.32  0.28  1.04  2.37
## devlgtpi[3,2]     168  445 0.38    161.82     0.36  1.40  0.39  1.01  2.40
## devlgtpsi[1,1] 101039   NA   NA 945897.89     0.00  4.80  4.78  0.02  4.82
## devlgtpsi[1,2]   1247   NA   NA   9957.57     0.00  3.02  2.82  0.20  3.23
## devlgtpsi[1,3]     94   NA   NA    148.26     0.12  0.96  0.14  0.82  1.78
## devlgtpsi[2,1] 331773   NA   NA 806346.46     0.00  1.31  1.27  0.04  1.36
## devlgtpsi[2,2]   2256   NA   NA  12249.80     0.00  2.36  2.10  0.26  2.62
## devlgtpsi[2,3]    576   NA   NA    908.32     0.06  0.93  0.17  0.76  1.70
## devlgtpsi[3,1]   6982   NA   NA  48945.22     0.00  3.69  3.62  0.06  3.75
## devlgtpsi[3,2]    210   NA   NA   2118.55     0.00  4.33  4.29  0.04  4.38
## devlgtpsi[3,3]    144   NA   NA    152.86     0.10  0.94  0.06  0.88  1.81
## TOTAL              NA   NA   NA        NA       NA 35.30 20.20 15.10 50.40

## Compared to the original model, DIC is around 30 smaller, no more negative pD, posterior mean deviances Dbar are closer to 1 for each data point. However the posteriors for 1 / log(r) for the GP and HOS levels are still concentrated towards the top end of the prior range, so we will try increasing the upper bound even further, to allow for even more over-dispersion:
fluDatMoreOD$odU <- 10


## Run in R2OpenBUGS
print(ptm <- proc.time())
outFluMoreOD <- bugs(model = "fluModel.txt",
               data = fluDatMoreOD,
               inits = fluInits,
               n.chains = 2,
               parameters.to.save = params,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE)
print(proc.time() - ptm)

outFluMoreODmc <- as.mcmc.list(outFluMoreOD)

## Look at posterior distribution of the over-dispersion parameter r, on original scale.
odRmoreOD <- getNode(outFluMoreODmc, "^r\\[")
sumStats(odRmoreOD)
chkTraces(outFluMoreODmc, "^r\\[")

## Look at posterior distribution of 1 / log(r) and compare to its Uniform(0.05,3) prior.
odinvlogrmoreOD <- getNode(outFluMoreODmc, "invlogr")
sumStats(odinvlogrmoreOD)
chkTraces(outFluMoreODmc, "invlogr")



## Look at deviance summaries for this model on logit/log scale
devsFluMoreODLgt <- getDevs(post = outFluMoreODmc,
                            daty = list(fluDatMoreOD$yp, fluDatMoreOD$ypi, fluDatMoreOD$y),
                            datn = list(fluDatMoreOD$np, fluDatMoreOD$npi, rep(NA, length(fluDatMoreOD$y))),
                            pNames = c("^lgtp\\[.,[1|2]\\]","^lgtpi","^lgtpsi\\["),
                            likTypes = c("binomial.lgt","binomial.lgt","negbin.lgt"),
                            aux = c(NA, NA, "^logr\\["),
                            devstem = "devlgt",
                            yhatstem = "yhat",
                            plugin = "mean")
round(devsFluMoreODLgt, digits = 2)
##                     y    n pobs   yexpbar thetabar  Dbar  Dhat    pD   DIC
## devlgtp[1,1]       41  883 0.05     41.84     0.05  1.01  0.00  1.01  2.01
## devlgtp[1,2]       41  210 0.20     41.62     0.20  1.00  0.00  1.00  2.00
## devlgtp[2,1]      350 1722 0.20    350.40     0.20  0.99  0.00  0.99  1.99
## devlgtp[2,2]      374 1676 0.22    374.30     0.22  0.99  0.00  0.99  1.98
## devlgtp[3,1]       22  168 0.13     23.12     0.14  0.98  0.03  0.95  1.93
## devlgtp[3,2]       13   52 0.25     13.49     0.26  0.98  0.01  0.97  1.95
## devlgtpi[1,1]      10  359 0.03     11.17     0.03  1.03  0.05  0.97  2.00
## devlgtpi[1,2]     236  423 0.56    233.81     0.55  1.06  0.04  1.02  2.08
## devlgtpi[2,1]      66  495 0.13     68.30     0.14  1.10  0.06  1.04  2.14
## devlgtpi[2,2]     379  978 0.39    375.70     0.38  1.04  0.05  0.99  2.03
## devlgtpi[3,1]     128  549 0.23    132.37     0.24  1.18  0.17  1.01  2.20
## devlgtpi[3,2]     168  445 0.38    163.59     0.37  1.16  0.20  0.96  2.12
## devlgtpsi[1,1] 101039   NA   NA 954827.47     0.00  3.70  3.67  0.03  3.74
## devlgtpsi[1,2]   1247   NA   NA  14154.38     0.00  2.95  2.73  0.22  3.17
## devlgtpsi[1,3]     94   NA   NA    235.73     0.04  1.05  0.43  0.62  1.67
## devlgtpsi[2,1] 331773   NA   NA 810503.06     0.00  1.01  0.97  0.04  1.04
## devlgtpsi[2,2]   2256   NA   NA  18316.16     0.00  2.47  2.23  0.24  2.71
## devlgtpsi[2,3]    576   NA   NA   1358.01     0.02  1.01  0.42  0.59  1.61
## devlgtpsi[3,1]   6982   NA   NA  53706.87     0.00  3.06  3.00  0.07  3.13
## devlgtpsi[3,2]    210   NA   NA   2674.69     0.00  3.69  3.60  0.09  3.79
## devlgtpsi[3,3]    144   NA   NA    183.66     0.03  0.75  0.01  0.74  1.49
## TOTAL              NA   NA   NA        NA       NA 32.23 17.68 14.55 46.79


## Now the prior (Uniform(0.05,10)) and posterior distributions of 1 / log(r) are more similar. DIC has not improved much over the Uniform(0.05,3) model.



## Finally, we look at how important allowing for over-dispersion in this example is, by reducing the upper bound of the 1 / log(r) prior to a small value, so that the likelihood approaches Poisson rather than negative binomial (alternatively, you could change the model code to use a Poisson likelihood rather than the negative binomial).
fluDatLessOD <- fluDat
fluDatLessOD$odU <- 0.1



## Need to amend initial values so lie in prior range
fluInitsLessOD <- list(
    list(
        p = t(structure(.Data = c(
                             0.3, 0.4, 0.1, 0.1, 0.5, 0.5,
                             0.3, 0.4, 0.1, 0.1, 0.5, 0.5,
                             0.3, 0.4, 0.1, 0.1, 0.5, 0.5
                             ), .Dim = c(6,3))),
        pi = t(structure(.Data = c(
                              0.1,NA,
                              0.1,NA,
                              0.1,NA
                          ), .Dim = c(2,3))),
        invlogr = c(0.05,0.07,0.09)
    ),
    
    list(
        
        p = t(structure(.Data = c(
                             0.5, 0.2, 0.2, 0.1, 0.5, 0.5,
                             0.5, 0.2, 0.2, 0.1, 0.5, 0.5,
                             0.5, 0.2, 0.2, 0.1, 0.5, 0.5
                             ), .Dim = c(6,3))),
        pi = t(structure(.Data = c(
                              0.5,NA,
                              0.5,NA,
                              0.5,NA
                          ), .Dim = c(2,3))),
        invlogr = c(0.075,0.075,0.075)
    )
)



## Run in R2OpenBUGS
print(ptm <- proc.time())
outFluLessOD <- bugs(model = "fluModel.txt",
               data = fluDatLessOD,
               inits = fluInitsLessOD,
               n.chains = 2,
               parameters.to.save = params,
               n.burnin = 1000,
               n.iter = 9000,
               n.thin = 1,
               DIC = TRUE)
print(proc.time() - ptm)

outFluLessODmc <- as.mcmc.list(outFluLessOD)

## Note the poor convergence/mixing for some of the numbers of infections at different severity levels, N:
chkTraces(outFluLessODmc, "^N\\[")

