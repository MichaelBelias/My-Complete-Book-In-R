library(rjags)

## Set of functions for processing posterior samples in R
## including functions for calculating deviance summaries
source("./fns.R")


#############################################
## QUESTION 1 - HIV model                  ##
##            - explore evidence synthesis ##
#############################################

## Model A: using single datum (y1), flat priors
hivModelA <- "
model
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
"

## DATA - datum implies maximum likelihood estimate of 0.05 for pi
hivDatA <- list(
  y = c(  5),
  n = c(100),
  a.pi = 1, b.pi = 1,
  a.delta = 1, b.delta = 1
)

## INITIAL VALUES FOR 2 CHAINS
hivInits <- list(
  ## chain 1
  list(
    pi = 0.1,
    delta = 0.9,
    .RNG.name = c("base::Mersenne-Twister"),
    .RNG.seed = c(7195)
  ),
  ## chain 2
  list(
    pi = 0.2,
    delta = 0.2,
    .RNG.name = c("base::Mersenne-Twister"),
    .RNG.seed = c(168422)
  )
)
  


## Parameters to monitor
## Basic & functional parameters p
hivParams <- c("p")


## Initialise model
hiv.jm <- jags.model(textConnection(hivModelA),
                     data = hivDatA,
                     inits = hivInits,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(hiv.jm, n.iter = 1000)
## samples to keep
outHIVA1 <- coda.samples(hiv.jm,
                       variable.names = hivParams,
                       n.iter = 9000,
                       n.thin = 1)
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
# Iterations = 1001:10000
# Thinning interval = 1 
# Number of chains = 2 
# Sample size per chain = 9000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean      SD  Naive SE Time-series SE
# p[1] 0.05888 0.02331 0.0001738      0.0001754
# p[2] 0.02948 0.02176 0.0001622      0.0001622
# p[3] 0.49946 0.28844 0.0021499      0.0021689
# 
# 2. Quantiles for each variable:
#   
#   2.5%     25%     50%     75%   97.5%
# p[1] 0.021828 0.04178 0.05601 0.07274 0.11237
# p[2] 0.001222 0.01246 0.02529 0.04209 0.08091
# p[3] 0.025687 0.24799 0.50207 0.74627 0.97500

## Summary statistics for a particular parameter only, using functions provided
## in fns.R - second argument of getNode() function is a regular expression that can 
## be passed to the grep command, e.g. to summarise pi = p[1]:
sumStats(as.matrix(getNode(outHIVA1, node = "p\\[1\\]"), ncol = 1))
# [,1]
# Mean     0.05888075
# SD       0.02331119
# Median   0.05601417
# 2.5%ile  0.02182798
# 97.5%ile 0.11236542


## Add in informative prior for delta (mean 0.75) and re-run
hivDatA$a.delta <- 75
hivDatA$b.delta <- 25

## Initialise model
hiv.jm <- jags.model(textConnection(hivModelA),
                     data = hivDatA,
                     inits = hivInits,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(hiv.jm, n.iter = 1000)
## samples to keep
outHIVA2 <- coda.samples(hiv.jm,
                       variable.names = hivParams,
                       n.iter = 9000,
                       n.thin = 1)
print(proc.time() - ptm)

## Summary statistics - note that all parameters (basic + functional) are now identified
sumStats(as.matrix(getNode(outHIVA2, node = "p"), ncol = 1))
#                p[1]        p[2]       p[3]
# Mean     0.05887241 0.014722985 0.74983274
# SD       0.02319329 0.006397553 0.04313006
# Median   0.05588192 0.013699807 0.75177879
# 2.5%ile  0.02241137 0.005205073 0.66038731
# 97.5%ile 0.11177326 0.029685127 0.82999782



## Add in informative prior for pi (mean 0.15) and re-run
hivDatA$a.pi <- 15
hivDatA$b.pi <- 85

## Initialise model
hiv.jm <- jags.model(textConnection(hivModelA),
                     data = hivDatA,
                     inits = hivInits,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(hiv.jm, n.iter = 1000)
## samples to keep
outHIVA3 <- coda.samples(hiv.jm,
                        variable.names = hivParams,
                        n.iter = 9000,
                        n.thin = 1)
print(proc.time() - ptm)

## Summary statistics - the informative prior for pi implies a different region of 
## support to the data, so note the posterior for pi = p[1] is a compromise between
## prior and likelihood
sumStats(as.matrix(getNode(outHIVA3, node = "p"), ncol = 1))
#                p[1]        p[2]       p[3]
# Mean     0.10015104 0.025023359 0.75005800
# SD       0.02114698 0.006869288 0.04320262
# Median   0.09864083 0.024255719 0.75174263
# 2.5%ile  0.06261560 0.013743491 0.65981874
# 97.5%ile 0.14497874 0.040277683 0.82986808


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
hivModelB <- "
model
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
"

## DATA
hivDatB <- list(
  y = c(  5,   3000),
  n = c(100, 100000),
  a.pi = 1, b.pi = 1,
  a.delta = 1, b.delta = 1
)

## Initialise model
hiv.jm <- jags.model(textConnection(hivModelB),
                     data = hivDatB,
                     inits = hivInits,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(hiv.jm, n.iter = 1000)
## samples to keep
outHIVB1 <- coda.samples(hiv.jm,
                         variable.names = hivParams,
                         n.iter = 9000,
                         n.thin = 1)
print(proc.time() - ptm)

## Summary statistics - all parameters are now identifiable even with flat priors,
## although the information on delta = p[3] is quite uncertain, since it is indirect
sumStats(as.matrix(getNode(outHIVB1, node = "p"), ncol = 1))
#                p[1]         p[2]       p[3]
# Mean     0.05521550 0.0300026360 0.41135028
# SD       0.01609751 0.0005376512 0.16384056
# Median   0.05219901 0.0299955995 0.42516734
# 2.5%ile  0.03197173 0.0289625289 0.06009068
# 97.5%ile 0.09294816 0.0310736516 0.67776440



## Check traces - note that pi = p[1] and delta = p[3] are now not mixing well. This is because the introduction of the data y2 informing pi(1-delta) induces correlation between pi and delta.
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
hivModelC <- "
model
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
"

## DATA
hivDatC <- list(
  y = c(  5,   3000,  90),
  n = c(100, 100000, 100),
  a.pi = 1, b.pi = 1,
  a.delta = 1, b.delta = 1
)

## Initialise model
hiv.jm <- jags.model(textConnection(hivModelC),
                     data = hivDatC,
                     inits = hivInits,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(hiv.jm, n.iter = 1000)
## samples to keep
outHIVC1 <- coda.samples(hiv.jm,
                         variable.names = hivParams,
                         n.iter = 9000,
                         n.thin = 1)
print(proc.time() - ptm)

## Summary statistics - the three data points are not consistent, hence
## the estimate of pi = p[1] is higher and more uncertain than the 
## corresponding estimate in model B1. Note that since the data on
## prevalence has the smallest sample size, it has the least influence
## and hence the estimate of pi is closer to that implied by y2, y3
## than that implied by y1
sumStats(as.matrix(getNode(outHIVC1, node = "p"), ncol = 1))
#                p[1]         p[2]       p[3]
# Mean     0.15431331 0.0298929180 0.80211580
# SD       0.02272932 0.0005375584 0.02905723
# Median   0.15230689 0.0298837451 0.80379154
# 2.5%ile  0.11383399 0.0288607067 0.73712022
# 97.5%ile 0.20431686 0.0309659655 0.85359963


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




########################################################################################################################################################



###################################################
## QUESTION 2 - 'Flu severity model              ##
##            - exploration & deviance summaries ##
###################################################

## Model code
fluModel <- "
model
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
"

## Set up data list
INF <- 1; SYM <- 2; GP <- 3; HOS <- 4; ICU <- 5; DEA <- 6;

fluDat <- list(
    A   = 3,  ## Number of age groups
    INF = INF,  ## Severity levels
    SYM = SYM,
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
                              0.2,NA,
                              0.2,NA,
                              0.2,NA
                          ), .Dim = c(2,3))),
        invlogr = c(0.1,0.1,0.1)
    )
)


## Parameters to monitor
## Key parameters N (number at each severity level); c (conditional probabilities)
## Plus nodes monitored to use in calculating deviance summaries
params <- c("N","p","pi","invlogr","logr","r","psi","lgtpi","lgtpsi","lgtp","yhatpi",
            "devpi","devlgtpi","yhatpsi","devpsi","devlgtpsi","yhatp","devp","devlgtp")


## Initialise model
flu.jm <- jags.model(textConnection(fluModel),
                     data = fluDat,
                     inits = fluInits,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(flu.jm, n.iter = 1000)
## samples to keep
outFlu <- coda.samples(flu.jm,
               variable.names = params,
               n.iter = 9000,
               n.thin = 1)
## in-built calculation of DIC, using Plummer's definition of pD
dicFlu_pD <- dic.samples(flu.jm, n.iter = 9000, thin = 1, type = "pD")
## in-built calculation of penalised expected deviance = Dbar + popt, where popt
## is the optimism (Plummer 2008)
dicFlu_popt <- dic.samples(flu.jm, n.iter = 9000, thin = 1, type = "popt")
print(proc.time() - ptm)


## pD and DIC based on Plummer (2008)
dicFlu_pD
# Mean deviance:  269.7 
# penalty 25.09 
# Penalized deviance: 294.8

## penalised expected deviance based on Plummer (2008)
dicFlu_popt
# Mean deviance:  269.7 
# penalty 137.4 
# Penalized deviance: 407.1


## Check traces for convergence - the string in the second argument is one that can be passed to the grep command to extract particular parameters from the mcmc.list
chkTraces(outFlu, "^pi")
dev.new()
chkTraces(outFlu, "^N")
chkTraces(outFlu, "^p\\[")
chkTraces(outFlu, "^r")
chkTraces(outFlu, "invlogr")
chkTraces(outFlu, "^psi")

## Summary statistics
summary(outFlu)


## Calculate deviance summaries, using posterior mean for plug-in
devsFluMean <- getDevs(post = outFlu,
                       daty = list(t(fluDat$yp), t(fluDat$ypi), t(fluDat$y)),
                       datn = list(t(fluDat$np), t(fluDat$npi), t(rep(NA, length(fluDat$y)))),
                       pNames = c("^p\\[.,[5|6]\\]","^pi","^psi\\["),
                       likTypes = c("binomial","binomial","negbin"),
                       aux = c(NA, NA, "^r\\["),
                       devstem = "dev[^i|^lgt]",
                       yhatstem = "yhat",
                       plugin = "mean")
round(devsFluMean, digits = 2)
## Notice that, as the posterior distributions of some parameters are skewed, the posterior mean is not a good summary of the posterior, hence the strange large values for the plug-in deviance Dhat and therefore the large negative values of pD and DIC.
##                  y    n pobs   yexpbar thetabar  Dbar        Dhat           pD          DIC
## devp[1,1]       41  883 0.05     44.16     0.05  1.11        0.24         0.87         1.98
## devp[1,2]       41  210 0.20     41.60     0.20  0.99        0.01         0.98         1.97
## devp[2,1]      350 1722 0.20    351.83     0.20  1.00        0.01         0.98         1.98
## devp[2,2]      374 1676 0.22    374.43     0.22  1.00        0.00         1.00         2.00
## devp[3,1]       22  168 0.13     27.47     0.16  2.15        1.38         0.78         2.93
## devp[3,2]       13   52 0.25     13.51     0.26  0.96        0.03         0.94         1.90
## devpi[1,1]      10  359 0.03     12.06     0.03  1.34        0.39         0.95         2.29
## devpi[1,2]     236  423 0.56    228.02     0.54  1.64        0.61         1.03         2.67
## devpi[2,1]      66  495 0.13     71.82     0.15  1.51        0.56         0.95         2.46
## devpi[2,2]     379  978 0.39    369.69     0.38  1.36        0.38         0.98         2.34
## devpi[3,1]     128  549 0.23    146.48     0.27  3.92        3.27         0.65         4.58
## devpi[3,2]     168  445 0.38    147.82     0.33  4.64        4.04         0.60         5.25
## devpsi[1,1] 101039   NA   NA 803978.23     0.00 18.76       18.51         0.25        19.00
## devpsi[1,2]   1247   NA   NA   2225.76     0.00  2.46        1.73         0.74         3.20
## devpsi[1,3]     94   NA   NA     84.27     0.32  1.09  4010155.85  -4010154.75  -4010153.66
## devpsi[2,1] 331773   NA   NA 686429.85     0.00  3.31        2.93         0.38         3.69
## devpsi[2,2]   2256   NA   NA   3579.29     0.00  1.74        1.10         0.64         2.38
## devpsi[2,3]    576   NA   NA    540.55     0.17  0.93  6328082.34  -6328081.41  -6328080.48
## devpsi[3,1]   6982   NA   NA  13702.58     0.00  2.99        1.54         1.45         4.44
## devpsi[3,2]    210   NA   NA    801.65     0.01  9.71        8.74         0.96        10.67
## devpsi[3,3]    144   NA   NA    109.72     0.29  2.47  4381827.51  -4381825.04  -4381822.58
## TOTAL           NA   NA   NA        NA       NA 65.09 14720111.17 -14720046.08 -14719981.00



## Calculate deviance summaries, using posterior median for plug-in
devsFluMedian <- getDevs(post = outFlu,
                         daty = list(t(fluDat$yp), t(fluDat$ypi), t(fluDat$y)),
                         datn = list(t(fluDat$np), t(fluDat$npi), t(rep(NA, length(fluDat$y)))),
                         pNames = c("^p\\[.,[5|6]\\]","^pi","^psi\\["),
                         likTypes = c("binomial","binomial","negbin"),
                         aux = c(NA, NA, "^r\\["),
                         devstem = "dev[^i|^lgt]",
                         yhatstem = "yhat",
                         plugin = "median")
round(devsFluMedian, digits = 2)
## The posterior median is a better summary, hence we obtain a much more reasonable estimate of Dhat. Note however that we still obtain negative values of the effective number of parameters, pD, which can occur when posterior is skewed or if sampling distribution is not log-concave, or if there is strong prior data conflict (see next lecture).
##                  y    n pobs   yexpbar thetabar  Dbar  Dhat    pD   DIC
## devp[1,1]       41  883 0.05     43.86     0.05  1.11  0.20  0.91  2.02
## devp[1,2]       41  210 0.20     41.40     0.20  0.99  0.00  0.98  1.97
## devp[2,1]      350 1722 0.20    351.53     0.20  1.00  0.01  0.99  1.98
## devp[2,2]      374 1676 0.22    374.12     0.22  1.00  0.00  1.00  2.00
## devp[3,1]       22  168 0.13     27.24     0.16  2.15  1.27  0.88  3.04
## devp[3,2]       13   52 0.25     13.35     0.26  0.96  0.01  0.95  1.92
## devpi[1,1]      10  359 0.03     11.71     0.03  1.34  0.27  1.07  2.41
## devpi[1,2]     236  423 0.56    228.21     0.54  1.64  0.58  1.06  2.70
## devpi[2,1]      66  495 0.13     71.64     0.14  1.51  0.53  0.98  2.49
## devpi[2,2]     379  978 0.39    369.64     0.38  1.36  0.38  0.98  2.33
## devpi[3,1]     128  549 0.23    146.47     0.27  3.92  3.27  0.65  4.58
## devpi[3,2]     168  445 0.38    147.64     0.33  4.64  4.11  0.53  5.18
## devpsi[1,1] 101039   NA   NA 796283.50     0.00 18.76 19.31 -0.55 18.21
## devpsi[1,2]   1247   NA   NA   2095.50     0.00  2.46  1.47  0.99  3.46
## devpsi[1,3]     94   NA   NA     84.00     0.17  1.09  0.70  0.40  1.49
## devpsi[2,1] 331773   NA   NA 675788.00     0.00  3.31  3.39 -0.08  3.23
## devpsi[2,2]   2256   NA   NA   3351.00     0.00  1.74  0.83  0.90  2.64
## devpsi[2,3]    576   NA   NA    544.00     0.03  0.93  0.42  0.51  1.45
## devpsi[3,1]   6982   NA   NA  12112.50     0.00  2.99  2.12  0.88  3.87
## devpsi[3,2]    210   NA   NA    768.00     0.01  9.71  8.37  1.34 11.04
## devpsi[3,3]    144   NA   NA    110.00     0.14  2.47  2.82 -0.36  2.11
## TOTAL           NA   NA   NA        NA       NA 65.09 50.07 15.02 80.10


## Calculate deviance summaries on log/lgt scale, using posterior mean for plug-in
devsFluLgt <- getDevs(post = outFlu,
                      daty = list(t(fluDat$yp), t(fluDat$ypi), t(fluDat$y)),
                      datn = list(t(fluDat$np), t(fluDat$npi), t(rep(NA, length(fluDat$y)))),
                      pNames = c("^lgtp\\[.,[1|2]\\]","^lgtpi","^lgtpsi\\["),
                      likTypes = c("binomial.lgt","binomial.lgt","negbin.lgt"),
                      aux = c(NA, NA, "^logr\\["),
                      devstem = "devlgt",
                      yhatstem = "yhat",
                      plugin = "mean")
round(devsFluLgt, digits = 2)
## Calculating the deviance on a scale for the parameters which is more likely to approach Gaussianity (e.g. logit scale for proportions, log scale for counts), helps to obtain values of pD that are positive. Note we still have one negative value for pD.
##                     y    n pobs   yexpbar thetabar  Dbar  Dhat    pD   DIC
## devlgtp[1,1]       41  883 0.05     44.16     0.05  1.11  0.19  0.92  2.04
## devlgtp[1,2]       41  210 0.20     41.60     0.20  0.99  0.00  0.99  1.98
## devlgtp[2,1]      350 1722 0.20    351.83     0.20  1.00  0.01  0.99  1.98
## devlgtp[2,2]      374 1676 0.22    374.43     0.22  1.00  0.00  1.00  2.00
## devlgtp[3,1]       22  168 0.13     27.47     0.16  2.15  1.23  0.92  3.07
## devlgtp[3,2]       13   52 0.25     13.51     0.26  0.96  0.01  0.96  1.92
## devlgtpi[1,1]      10  359 0.03     12.06     0.03  1.34  0.22  1.12  2.46
## devlgtpi[1,2]     236  423 0.56    228.02     0.54  1.64  0.60  1.04  2.68
## devlgtpi[2,1]      66  495 0.13     71.82     0.15  1.51  0.50  1.01  2.52
## devlgtpi[2,2]     379  978 0.39    369.69     0.38  1.36  0.39  0.97  2.33
## devlgtpi[3,1]     128  549 0.23    146.48     0.27  3.92  3.21  0.71  4.63
## devlgtpi[3,2]     168  445 0.38    147.82     0.33  4.64  4.08  0.56  5.21
## devlgtpsi[1,1] 101039   NA   NA 803978.23     0.00 19.24 19.18  0.07 19.31
## devlgtpsi[1,2]   1247   NA   NA   2225.76     0.00  2.53  2.13  0.40  2.93
## devlgtpsi[1,3]     94   NA   NA     84.27     0.32  1.10  0.63  0.47  1.56
## devlgtpsi[2,1] 331773   NA   NA 686429.85     0.00  3.49  3.30  0.19  3.68
## devlgtpsi[2,2]   2256   NA   NA   3579.29     0.00  1.79  1.38  0.41  2.20
## devlgtpsi[2,3]    576   NA   NA    540.55     0.17  0.94  0.36  0.57  1.51
## devlgtpsi[3,1]   6982   NA   NA  13702.58     0.00  3.15  2.30  0.85  4.01
## devlgtpsi[3,2]    210   NA   NA    801.65     0.01  9.83  9.67  0.17 10.00
## devlgtpsi[3,3]    144   NA   NA    109.72     0.29  2.47  3.80 -1.33  1.13
## TOTAL              NA   NA   NA        NA       NA 66.17 53.18 12.99 79.16



##################################################
## QUESTION 3 - 'Flu severity model             ##
##            - Calculate case-severity risks   ##
##################################################



## Calculate case-severity risks from conditional probabilities. Since we already have posterior samples of the conditional probabilities p[a,l], we can calculate the functional parameters directly:
ps <- getNode(outFlu, "^p\\[")
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
fluModelCSR <- "
model
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
"

## Monitor case-severity risks
paramsCSR <- c("sCHR","sCIR","sCFR","CHR","CIR","CFR")


## GP, HOS and ICU constants now needed, so set in data file
fluDatCSR <- fluDat
fluDatCSR$GP <- GP
fluDatCSR$HOS <- HOS
fluDatCSR$ICU <- ICU


## Initialise model
flu.jm.CSR <- jags.model(textConnection(fluModelCSR),
                     data = fluDatCSR,
                     inits = fluInits,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(flu.jm.CSR, n.iter = 1000)
## samples to keep
outFluCSR <- coda.samples(flu.jm.CSR,
               variable.names = paramsCSR,
               n.iter = 9000,
               n.thin = 1)
print(proc.time() - ptm)

## Check convergence
chkTraces(outFluCSR, "C")

## Posterior summaries
CSRs <- getNode(outFluCSR, "C")
sumStats(CSRs)





##################################################
## QUESTION 4 - 'Flu severity model             ##
##            - Explore over-dispersion         ##
##################################################



## Look at posterior distribution of the over-dispersion parameter r, on original scale.
odR <- getNode(outFlu, "^r\\[")
sumStats(odR)
chkTraces(outFlu, "^r\\[")

## Look at posterior distribution of 1 / log(r) and compare to its Uniform(0.05,0.5) prior.
odinvlogr <- getNode(outFlu, "invlogr")
sumStats(odinvlogr)
chkTraces(outFlu, "invlogr")

## Note that a smaller value of r implies more over-dispersion, whereas a larger value of r imples less over-dispersion. So conversely, a smaller value of 1 / log(r) implies less over-dispersion and a larger values implies more. Note that for the first two severity levels (GP, HOS), the posterior of 1 / log(r) is concentrated towards the upper bound (0.5) of its prior distribution, suggesting we should allow more over-dispersion for these levels than the prior allows. By contrast, the DEA level's posterior is more evenly spread over the prior space, although still with a mode towards the upper end of the prior range.



## Allow for more over-dispersion, by increasing the upper bound of the Uniform prior for 1 / log(r)
fluDatMoreOD <- fluDat
fluDatMoreOD$odU <- 3



## Initialise model
flu.jm.moreOD <- jags.model(textConnection(fluModel),
                     data = fluDatMoreOD,
                     inits = fluInits,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(flu.jm.moreOD, n.iter = 1000)
## samples to keep
outFluMoreOD <- coda.samples(flu.jm.moreOD,
               variable.names = params,
               n.iter = 9000,
               n.thin = 1)
print(proc.time() - ptm)


## Look at posterior distribution of the over-dispersion parameter r, on original scale.
odRmoreOD <- getNode(outFluMoreOD, "^r\\[")
sumStats(odRmoreOD)
chkTraces(outFluMoreOD, "^r\\[")

## Look at posterior distribution of 1 / log(r) and compare to its Uniform(0.05,3) prior.
odinvlogrmoreOD <- getNode(outFluMoreOD, "invlogr")
sumStats(odinvlogrmoreOD)
chkTraces(outFluMoreOD, "invlogr")


## Look at deviance summaries for this model on logit/log scale
devsFluMoreODLgt <- getDevs(post = outFluMoreOD,
                            daty = list(t(fluDatMoreOD$yp), t(fluDatMoreOD$ypi), t(fluDatMoreOD$y)),
                            datn = list(t(fluDatMoreOD$np), t(fluDatMoreOD$npi), t(rep(NA, length(fluDatMoreOD$y)))),
                            pNames = c("^lgtp\\[.,[1|2]\\]","^lgtpi","^lgtpsi\\["),
                            likTypes = c("binomial.lgt","binomial.lgt","negbin.lgt"),
                            aux = c(NA, NA, "^logr\\["),
                            devstem = "devlgt",
                            yhatstem = "yhat",
                            plugin = "mean")
round(devsFluMoreODLgt, digits = 2)
##                     y    n pobs   yexpbar thetabar  Dbar  Dhat    pD   DIC
## devlgtp[1,1]       41  883 0.05     41.83     0.05  1.01  0.00  1.00  2.01
## devlgtp[1,2]       41  210 0.20     41.52     0.20  0.99  0.00  0.99  1.97
## devlgtp[2,1]      350 1722 0.20    350.48     0.20  0.97  0.00  0.97  1.93
## devlgtp[2,2]      374 1676 0.22    374.73     0.22  0.99  0.00  0.99  1.99
## devlgtp[3,1]       22  168 0.13     23.39     0.14  1.03  0.05  0.98  2.02
## devlgtp[3,2]       13   52 0.25     13.50     0.26  0.97  0.01  0.97  1.94
## devlgtpi[1,1]      10  359 0.03     11.15     0.03  1.06  0.04  1.02  2.08
## devlgtpi[1,2]     236  423 0.56    233.86     0.55  1.06  0.04  1.02  2.08
## devlgtpi[2,1]      66  495 0.13     68.36     0.14  1.07  0.07  1.00  2.07
## devlgtpi[2,2]     379  978 0.39    375.88     0.38  1.05  0.05  1.01  2.06
## devlgtpi[3,1]     128  549 0.23    133.83     0.24  1.35  0.31  1.04  2.39
## devlgtpi[3,2]     168  445 0.38    161.62     0.36  1.42  0.41  1.01  2.42
## devlgtpsi[1,1] 101039   NA   NA 951292.62     0.00  4.80  4.78  0.02  4.82
## devlgtpsi[1,2]   1247   NA   NA   8672.32     0.00  2.95  2.77  0.18  3.13
## devlgtpsi[1,3]     94   NA   NA    144.40     0.10  0.93  0.11  0.82  1.75
## devlgtpsi[2,1] 331773   NA   NA 810753.87     0.00  1.32  1.28  0.04  1.36
## devlgtpsi[2,2]   2256   NA   NA  12122.74     0.00  2.44  2.20  0.24  2.68
## devlgtpsi[2,3]    576   NA   NA    900.15     0.05  0.94  0.15  0.79  1.73
## devlgtpsi[3,1]   6982   NA   NA  48153.07     0.00  3.62  3.55  0.06  3.68
## devlgtpsi[3,2]    210   NA   NA   2101.66     0.00  4.32  4.32  0.00  4.33
## devlgtpsi[3,3]    144   NA   NA    151.73     0.09  0.94  0.06  0.88  1.82
## TOTAL              NA   NA   NA        NA       NA 35.25 20.21 15.04 50.28

## Compared to the original model, DIC is around 30 smaller, no more negative pD, posterior mean deviances Dbar are closer to 1 for each data point. However the posteriors for 1 / log(r) for the GP and HOS levels are still concentrated towards the top end of the prior range, so we will try increasing the upper bound even further, to allow for even more over-dispersion:
fluDatMoreOD$odU <- 10



## Initialise model
flu.jm.moreOD <- jags.model(textConnection(fluModel),
                     data = fluDatMoreOD,
                     inits = fluInits,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(flu.jm.moreOD, n.iter = 1000)
## samples to keep
outFluMoreOD <- coda.samples(flu.jm.moreOD,
               variable.names = params,
               n.iter = 9000,
               n.thin = 1)
print(proc.time() - ptm)


## Look at posterior distribution of the over-dispersion parameter r, on original scale.
odRmoreOD <- getNode(outFluMoreOD, "^r\\[")
sumStats(odRmoreOD)
chkTraces(outFluMoreOD, "^r\\[")

## Look at posterior distribution of 1 / log(r) and compare to its Uniform(0.05,10) prior.
odinvlogrmoreOD <- getNode(outFluMoreOD, "invlogr")
sumStats(odinvlogrmoreOD)
chkTraces(outFluMoreOD, "invlogr")



## Look at deviance summaries for this model on logit/log scale
devsFluMoreODLgt <- getDevs(post = outFluMoreOD,
                            daty = list(t(fluDatMoreOD$yp), t(fluDatMoreOD$ypi), t(fluDatMoreOD$y)),
                            datn = list(t(fluDatMoreOD$np), t(fluDatMoreOD$npi), t(rep(NA, length(fluDatMoreOD$y)))),
                            pNames = c("^lgtp\\[.,[1|2]\\]","^lgtpi","^lgtpsi\\["),
                            likTypes = c("binomial.lgt","binomial.lgt","negbin.lgt"),
                            aux = c(NA, NA, "^logr\\["),
                            devstem = "devlgt",
                            yhatstem = "yhat",
                            plugin = "mean")
round(devsFluMoreODLgt, digits = 2)
##                     y    n pobs   yexpbar thetabar  Dbar  Dhat    pD   DIC
## devlgtp[1,1]       41  883 0.05     41.77     0.05  0.98  0.00  0.98  1.97
## devlgtp[1,2]       41  210 0.20     41.61     0.20  0.98  0.00  0.98  1.96
## devlgtp[2,1]      350 1722 0.20    350.45     0.20  0.99  0.00  0.99  1.98
## devlgtp[2,2]      374 1676 0.22    374.32     0.22  0.97  0.00  0.97  1.93
## devlgtp[3,1]       22  168 0.13     23.06     0.14  1.03  0.02  1.01  2.04
## devlgtp[3,2]       13   52 0.25     13.48     0.26  0.96  0.01  0.95  1.91
## devlgtpi[1,1]      10  359 0.03     11.14     0.03  1.04  0.04  1.00  2.04
## devlgtpi[1,2]     236  423 0.56    234.24     0.55  1.04  0.03  1.02  2.06
## devlgtpi[2,1]      66  495 0.13     68.02     0.14  1.04  0.05  0.99  2.03
## devlgtpi[2,2]     379  978 0.39    376.66     0.39  1.02  0.03  0.99  2.01
## devlgtpi[3,1]     128  549 0.23    132.53     0.24  1.22  0.18  1.04  2.26
## devlgtpi[3,2]     168  445 0.38    163.31     0.37  1.23  0.22  1.00  2.23
## devlgtpsi[1,1] 101039   NA   NA 958367.51     0.00  3.71  3.68  0.03  3.75
## devlgtpsi[1,2]   1247   NA   NA  13687.53     0.00  2.99  2.79  0.21  3.20
## devlgtpsi[1,3]     94   NA   NA    233.26     0.04  1.05  0.45  0.60  1.66
## devlgtpsi[2,1] 331773   NA   NA 820274.59     0.00  1.02  0.99  0.04  1.06
## devlgtpsi[2,2]   2256   NA   NA  17758.92     0.00  2.44  2.20  0.24  2.68
## devlgtpsi[2,3]    576   NA   NA   1328.19     0.02  1.00  0.43  0.57  1.57
## devlgtpsi[3,1]   6982   NA   NA  53218.55     0.00  3.03  2.96  0.07  3.10
## devlgtpsi[3,2]    210   NA   NA   2616.76     0.00  3.66  3.57  0.09  3.75
## devlgtpsi[3,3]    144   NA   NA    177.18     0.04  0.75  0.00  0.75  1.50
## TOTAL              NA   NA   NA        NA       NA 32.16 17.65 14.51 46.67


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
                              0.2,NA,
                              0.2,NA,
                              0.2,NA
                          ), .Dim = c(2,3))),
        invlogr = c(0.075,0.075,0.075)
    )
)


## Initialise model
flu.jm.lessOD <- jags.model(textConnection(fluModel),
                     data = fluDatLessOD,
                     inits = fluInitsLessOD,
                     n.chains = 2)

## Run in JAGS
print(ptm <- proc.time())
## burn-in
update(flu.jm.lessOD, n.iter = 1000)
## samples to keep
outFluLessOD <- coda.samples(flu.jm.lessOD,
               variable.names = params,
               n.iter = 9000,
               n.thin = 1)
print(proc.time() - ptm)


## Note the poor convergence/mixing for some of the numbers of infections at different severity levels, N:
chkTraces(outFluLessOD, "^N\\[")

