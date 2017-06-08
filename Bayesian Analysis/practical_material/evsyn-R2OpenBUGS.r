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
## Model as in slide 149 - FILL IN FUNCTIONAL PARAMETERS IN PLACE OF QUESTION MARKS
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
  p[1] <- ?
  p[2] <- ?
  p[3] <- ?
}



## Write model code to a text file
write.model(hivModelA, con = "hivModelA.txt")


## DATA - Replace question marks with data. You can either hard code
## the parameters of the Beta priors in the model code, or set them in
## the data (a.pi, b.pi, a.delta, b.delta).
hivDatA <- list(
  y = c(  ?, NA, NA),
  n = c(  ?, NA, NA),
  a.pi = ?, b.pi = ?,
  a.delta = ?, b.delta = ?
)

## INITIAL VALUES FOR 2 CHAINS - replace question marks with initial values
hivInits <- list(
  ## chain 1
  list(
    pi = ?,
    delta = ?
  ),
  ## chain 2
  list(
    pi = ?,
    delta = ?
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
chkTraces(outHIVA1, node = "p")

## Summary statistics - using rjags provided command to look at summaries of all
## monitored parameters
summary(outHIVA1)


## Summary statistics for a particular parameter only, using functions provided
## in fns.R - second argument of getNode() function is a regular expression that can 
## be passed to the grep command, e.g. to summarise pi = p[1]:
sumStats(as.matrix(getNode(outHIVA1, node = "p\\[1\\]"), ncol = 1))



## Add in informative prior for delta (mean 0.75) and re-run


## Add in informative prior for pi (mean 0.15) and re-run



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
?
}

## Write model code to a text file
write.model(hivModelB, con = "hivModelB.txt")

## DATA
hivDatB <- list(
  ?
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

## Summary statistics


## Check traces

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
  ?
}

## Write model code to a text file
write.model(hivModelC, con = "hivModelC.txt")


## DATA
hivDatC <- list(
  ?
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

## Summary statistics


## Check traces


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
outFlu$DIC

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




##################################################
## QUESTION 3 - 'Flu severity model             ##
##            - Calculate case-severity risks   ##
##################################################



## Calculate case-severity risks from conditional probabilities. Since we already have posterior samples of the conditional probabilities p[a,l], we can calculate the functional parameters directly:
ps <- getNode(outFluMC, "^p\\[")
colnames(ps)

## symptomatic case-severity risks are the probabilities of a severe event conditional on a symptomatic infection (level 2) - REPLACE question marks with appropriate functions of ps
sCHR <- ?
sCIR <- ?
sCFR <- ?

## case-severity risks are the probabilities of a severe event conditional on an infection (level 1)
## and can be calculated as a symptomatic case-severity risk multiplied by Pr{Symptomatic | Infected} - REPLACE question marks with appropriate functions of ps
CHR <-  ?
CIR <-  ?
CFR <-  ?


## look at estimates
sumStats(sCHR)
sumStats(sCIR)
sumStats(sCFR)
sumStats(CHR)
sumStats(CIR)
sumStats(CFR)


## Alternatively, amend the model code to calculate the case-severity risks within OpenBUGS - REPLACE question marks with appropriate functions

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
        sCHR[a] <- ?
        ## Pr(ICU | SYM) = Pr(ICU | HOS) * Pr(HOS | SYM)
        sCIR[a] <- ?
        ## Pr(DEA | SYM) = Pr(DEA | ICU) * Pr(ICU | SYM)
        sCFR[a] <- ?

        ## Corresponding infection-severity risks are symptomatic case-severity risks multiplied
        ## by Pr(SYM | INF)
        CHR[a] <- ?
        CIR[a] <- ?
        CFR[a] <- ?
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

## Posterior summaries





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


## Allow for more over-dispersion, by increasing the upper bound of the Uniform prior for 1 / log(r) - REPLACE question mark
fluDatMoreOD <- fluDat
fluDatMoreOD$odU <- ?


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




## Finally, we look at how important allowing for over-dispersion in this example is, by reducing the upper bound of the 1 / log(r) prior to a small value, so that the likelihood approaches Poisson rather than negative binomial (alternatively, you could change the model code to use a Poisson likelihood rather than the negative binomial) - REPLACE question mark.
fluDatLessOD <- fluDat
fluDatLessOD$odU <- ?



## NOTE depending on what value you chose for odU, you may need to amend initial values so lie in prior range - REPLACE question marks
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
        invlogr = c(?,?,?)
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
        invlogr = c(?,?,?)
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

## Check traces for N
chkTraces(outFluLessODmc, "^N\\[")

