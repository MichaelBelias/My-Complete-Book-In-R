## Solutions for practical1

## ---- Question A1 ----
#  file.show("coins-model.txt")

## ---- Question A3 ----
coins.jag <- jags.model("coins-model.txt")

## ---- Question A4 ----
coins.coda <- coda.samples(coins.jag, c("Y","P8"), n.iter=10000)

## ---- Question A5 ----
summary(coins.coda)
par(mfrow=c(2,2))
traceplot(coins.coda)

# densplot.discrete() provided in plot-functions.R
# CODA's default densplot() is misleading for discrete parameters
densplot.discrete(coins.coda)
coins.A5.jag <- jags.model("solutions/practical1/coins-model-A5.txt",
                           quiet=TRUE)
coins.A5.coda <- coda.samples(coins.A5.jag, c("Y", "P15"), n.iter=10000)
summary(coins.A5.coda)
par(mfrow=c(2, 2))
traceplot(coins.A5.coda)
densplot.discrete(coins.A5.coda)

## ---- Question B1 ----
drug.jag <- jags.model("drug-MC-model.txt", quiet=TRUE)
drug.coda <- coda.samples(drug.jag,
                          variable.names = c("theta","r","P.crit"),
                          n.iter=30000)
summary(drug.coda)
drug.coda <- drug.coda[,c("r","theta")]
par(mfrow=c(3,2))
traceplot(drug.coda)
densplot.discrete(drug.coda[, "r"]) # r is a discrete parameter
densplot(drug.coda[, "theta"]) # theta is a continuous parameter
autocorr.plot(drug.coda, auto.layout=FALSE)

## ---- Question B2 ----
drug.B2.jag <- jags.model("solutions/practical1/drug-MC-model-B2.txt",
                          quiet=TRUE)
drug.B2.coda <- coda.samples(drug.B2.jag,
                             variable.names=c("theta","r","P.crit"),
                             n.iter=30000)
drug.B2.coda.r <- drug.B2.coda[[1]][,"r"]
densplot.discrete(drug.B2.coda.r)
summary(drug.B2.coda)
drug.B2.coda.theta <- drug.B2.coda[,c("theta")]
densplot(drug.B2.coda.theta, ylim=c(0, 1.2), show.obs=FALSE)

## ---- Question C1 ----
predpower.jag <- jags.model(file="predpower-model.txt",
                            quiet=TRUE)
variable.names <- c("n","sigma","theta","power")
predpower.coda <- coda.samples(predpower.jag,
                               variable.names=variable.names,
                               n.iter=50000)
summary(predpower.coda)
par(mfrow=c(2, 2))
densplot(predpower.coda)

## ---- Question C2 ----
file <- "solutions/practical1/predpower-model-C2.txt"
predpower.C2.jag <- jags.model(file=file, quiet=TRUE)
variable.names <- c("n","sigma","theta","power")
predpower.C2.coda <- coda.samples(predpower.C2.jag,
                                  variable.names=variable.names,
                                  n.iter=50000)
par(mfrow=c(2, 2))
densplot(predpower.C2.coda)

## ---- Question D1 ----
file <- "solutions/practical1/t-model-D1.txt"
t.D1.jag <- jags.model(file=file, quiet=TRUE)
t.D1.coda <- coda.samples(t.D1.jag,
                          variable.names="y",
                          n.iter=10000)
summary(t.D1.coda)
densplot(t.D1.coda)

## ---- Question D2 ----
file <- "solutions/practical1/cubed-model-D2.txt"
cubed.D2.jag <- jags.model(file=file, quiet=TRUE)
cubed.D2.coda <- coda.samples(cubed.D2.jag,
                              variable.names=c("y","ycubed"),
                              n.iter=10000)
summary(cubed.D2.coda)
par(mfrow=c(1,2))
densplot(cubed.D2.coda)

