## Solutions for practical1

## ---- Question A1 ----
#  file.show("coins-model.txt")

## ---- Question A2 ----
coins.sim <- bugs(data="nodata.txt",
                  inits="nodata.txt",
                  model="coins-model.txt",
                  n.chains=1,
                  n.burnin=0,
                  n.iter=10000,
                  n.thin=1,
                  parameters.to.save=c("Y","P8"),
                  DIC=FALSE)

## ---- Question A2b ----
coins.sim

## ---- Question A3 ----
coins.coda <- as.mcmc.list(coins.sim)

## ---- Question A5 ----
summary(coins.coda)
par(mfrow=c(2,2))
traceplot(coins.coda)

# densplot.discrete() provided in plot-functions.R
# CODA's default densplot() is misleading for discrete parameters
densplot.discrete(coins.coda)
model <- normalizePath("solutions/practical1/coins-model-A5.txt")
coins.A5.sim <- bugs(data="nodata.txt",
                     inits="nodata.txt",
                     model=model,
                     n.chains=1,
                     n.burnin=0,
                     n.iter=10000,
                     n.thin=1,
                     parameters.to.save=c("Y","P15"),
                     DIC=FALSE)
coins.A5.coda <- as.mcmc.list(coins.A5.sim)
summary(coins.A5.coda)
par(mfrow=c(2, 2))
traceplot(coins.A5.coda)
densplot.discrete(coins.A5.coda)

## ---- Question B1 ----
drug.sim <- bugs(model="drug-MC-model.txt",
                 data="nodata.txt",
                 inits="nodata.txt",
                 n.chains=1,
                 n.burnin=0,
                 n.iter=30000,
                 n.thin=1,
                 parameters.to.save=c("theta","r","P.crit"),
                 DIC=FALSE)
drug.sim
drug.coda <- as.mcmc.list(drug.sim)
drug.coda <- drug.coda[,c("r","theta")]
par(mfrow=c(3,2))
traceplot(drug.coda)
densplot.discrete(drug.coda[, "r"]) # r is a discrete parameter
densplot(drug.coda[, "theta"]) # theta is a continuous parameter
autocorr.plot(drug.coda, auto.layout=FALSE)

## ---- Question B2 ----
model <- normalizePath("solutions/practical1/drug-MC-model-B2.txt")
drug.B2.sim <- bugs(data="nodata.txt",
                    inits="nodata.txt",
                    model=model,
                    n.chains=1,
                    n.burnin=0,
                    n.iter=30000,
                    n.thin=1,
                    parameters.to.save=c("theta","r","P.crit"),
                    DIC=FALSE)
drug.B2.coda <- as.mcmc.list(drug.B2.sim)
drug.B2.coda.r <- drug.B2.coda[[1]][,"r"]
densplot.discrete(drug.B2.coda.r)
summary(drug.B2.coda)
drug.B2.coda.theta <- drug.B2.coda[,c("theta")]
densplot(drug.B2.coda.theta, ylim=c(0, 1.2), show.obs=FALSE)

## ---- Question C1 ----
predpower.sim <- bugs(model="predpower-model.txt",
                      data="nodata.txt",
                      inits="nodata.txt",
                      n.chains=1,
                      n.burnin=0,
                      n.iter=50000,
                      n.thin=1,
                      parameters.to.save=c("n","sigma","theta","power"),
                      DIC=FALSE)
predpower.sim
predpower.coda <- as.mcmc.list(predpower.sim)
summary(predpower.coda)
par(mfrow=c(2, 2))
densplot(predpower.coda)

## ---- Question C2 ----
model <- normalizePath("solutions/practical1/predpower-model-C2.txt")
predpower.C2.sim <- bugs(data="nodata.txt",
                     inits="nodata.txt",
                     model=model,
                     n.chains=1,
                     n.burnin=0,
                     n.iter=50000,
                     n.thin=1,
                     parameters.to.save=c("n","sigma","theta","power"),
                     DIC=FALSE)
predpower.C2.coda <- as.mcmc.list(predpower.C2.sim)
par(mfrow=c(2, 2))
densplot(predpower.C2.coda)

## ---- Question D1 ----
model <- normalizePath("solutions/practical1/t-model-D1.txt")
t.D1.sim <- bugs(data="nodata.txt",
                 inits="nodata.txt",
                 model=model,
                 n.chains=1,
                 n.burnin=0,
                 n.iter=10000,
                 n.thin=1,
                 parameters.to.save="y",
                 DIC=FALSE)
t.D1.coda <- as.mcmc.list(t.D1.sim)
summary(t.D1.coda)
densplot(t.D1.coda)

## ---- Question D2 ----
model <- normalizePath("solutions/practical1/cubed-model-D2.txt")
cubed.D2.sim <- bugs(data="nodata.txt",
                     inits="nodata.txt",
                     model=model,
                     n.chains=1,
                     n.burnin=0,
                     n.iter=10000,
                     n.thin=1,
                     parameters.to.save=c("y","ycubed"),
                     DIC=FALSE)
cubed.D2.coda <- as.mcmc.list(cubed.D2.sim)
summary(cubed.D2.coda)
par(mfrow=c(1,2))
densplot(cubed.D2.coda)

