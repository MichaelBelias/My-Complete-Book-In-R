## Solutions for practical4

## ---- Question A1 ----
# create a second set of initial values
dugongs_in2_A1 <- list(alpha = 10, beta = 5, tau = 2, gamma = 0.2)

parameters.to.save <- c("alpha","beta","gamma","mu","sigma","res","P.res",
                        "P.pred","Y.pred")
dugongs.A1.sim <- bugs(model="dugongs-check-model.txt",
                       data=dugongs,
                       inits=list(dugongs_in1, dugongs_in2_A1),
                       n.chains=2,
                       n.burnin=5000,
                       n.iter=55000,
                       n.thin=1,
                       parameters.to.save=parameters.to.save,
                       DIC=FALSE)
dugongs.A1.coda <- as.mcmc.list(dugongs.A1.sim)
summary(dugongs.A1.coda[,c("alpha","beta","gamma","sigma")])

## ---- Question A2 ----
bugs.boxplot(dugongs.A1.coda, "res")
density.strips(dugongs.A1.coda, "res")

## ---- Question A2b ----
catplot(dugongs.A1.coda, "P.res")
catplot(dugongs.A1.coda, "P.res", ordered=TRUE)

## ---- Question A2c ----
res <- get.coda.matrix(dugongs.A1.coda, "P.pred")
qs <- apply(res, 2, mean)
plot(sort(qs), length(qs):1, pch=20, yaxt="n",
     xlab="Predictive p-value", ylab="")

## ---- Question A2d ----
model.fit(dugongs$Y, dugongs$x, dugongs.A1.coda, "Y.pred",
          pch=20, lwd=2, xlab="Age (years)", ylab="Length (m)",
          ylim=c(1.5,3))

## ---- Question A2e ----
res <- get.coda.matrix(dugongs.A1.coda, "res")
plot(dugongs$x, colMeans(res), pch=20, ylim=c(-4,4), xlim=c(0,40),
     col="blue", xlab="Age (years)", ylab="Posterior mean residual")

## ---- Question A2f ----
res <- get.coda.matrix(dugongs.A1.coda, "res")
fit <- get.coda.matrix(dugongs.A1.coda, "mu")
plot(colMeans(fit), colMeans(res), pch=20, col="blue",
     ylim=c(-4, 4), xlim=c(1.75, 2.75),
     xlab="Posterior mean fitted values (mu)", ylab="Posterior mean residual")

## ---- Question A3 ----
# create a second set of initial values
dugongs_in2_A1 <- list(alpha = 10, beta = 5, tau = 2, gamma = 0.2)

# edit the data to make predictions at age 35 and 40
dugongs_A3 <- dugongs
dugongs_A3$x <- c(dugongs_A3$x, 35, 40)
dugongs_A3$Y <- c(dugongs_A3$Y, NA, NA)
dugongs_A3$N <- 29

parameters.to.save <- c("alpha","beta","gamma","mu","sigma","res","P.res",
                        "P.pred","Y.pred", "Y[28]","Y[29]")
dugongs.A3.sim <- bugs(model="dugongs-check-model.txt",
                       data=dugongs_A3,
                       inits=list(dugongs_in1, dugongs_in2_A1),
                       n.chains=2,
                       n.burnin=5000,
                       n.iter=55000,
                       n.thin=1,
                       parameters.to.save=parameters.to.save,
                       DIC=FALSE)
dugongs.A3.coda <- as.mcmc.list(dugongs.A3.sim)
summary(dugongs.A3.coda[,c("Y[28]","Y[29]")])

## ---- Question A3b ----
model.fit(dugongs_A3$Y, dugongs_A3$x, dugongs.A3.coda, "mu",
          pch=20, lwd=2, xlab="Age (years)", ylab="Length (m)",
          ylim=c(1.5,3))

## ---- Question A4 ----
# create a second set of initial values
dugongs_in2_A4 <- list(alpha = 10, beta = 5, tau = 2, gamma = 0.2)

parameters.to.save <- c("alpha","beta","gamma","mu","sigma","res","P.res",
                        "P.pred","Y.pred")
dugongs.A4.sim <- bugs(model="dugongs-check-model.txt",
                       data=dugongs_outl,
                       inits=list(dugongs_in1, dugongs_in2_A4),
                       n.chains=2,
                       n.burnin=5000,
                       n.iter=55000,
                       n.thin=1,
                       parameters.to.save=parameters.to.save,
                       DIC=FALSE)
dugongs.A4.coda <- as.mcmc.list(dugongs.A4.sim)
summary(dugongs.A4.coda[,c("alpha","beta","gamma","sigma")])

## ---- Question A4b ----
model.fit(dugongs_outl$Y, dugongs_outl$x, dugongs.A4.coda, "mu",
          pch=20, lwd=2, xlab="Age (years)", ylab="Length (m)",
          ylim=c(1.5,3))

## ---- Question A4c ----
model.fit(dugongs_outl$Y, dugongs_outl$x, dugongs.A4.coda, "Y.pred",
          pch=20, lwd=2, xlab="Age (years)", ylab="Length (m)",
          ylim=c(1.5,3))

## ---- Question A4d ----
bugs.boxplot(dugongs.A4.coda, "res")
density.strips(dugongs.A4.coda, "res")

## ---- Question A4e ----
res <- get.coda.matrix(dugongs.A4.coda, "res")
plot(dugongs_outl$x, colMeans(res), pch=20, col="blue",
     ylim=c(-4,4), xlim=c(0,40),
     xlab="Age (years)", ylab="Posterior mean residual")

## ---- Question A4f ----
catplot(dugongs.A4.coda, "P.res", ordered=TRUE)

## ---- Question A4g ----
P.res.names <- paste0("P.res[", 1:27,"]")
summary(dugongs.A4.coda[, P.res.names])

## ---- Question A5 ----
dugongs_model_A5_file <- normalizePath("solutions/practical4/dugongs-model-A5.txt")
dugongs_in2_A5 <- list(alpha = 10, beta = 5, tau = 2, gamma = 0.2)
parameters.to.save <- c("alpha","beta","gamma","mu","sigma","res","P.res",
                        "P.pred","Y.pred")

dugongs.A5.sim <- bugs(model=dugongs_model_A5_file,
                       data=dugongs_outl,
                       inits=list(dugongs_in1, dugongs_in2_A5),
                       n.chains=2,
                       n.burnin=5000,
                       n.iter=55000,
                       n.thin=1,
                       parameters.to.save=parameters.to.save,
                       DIC=FALSE)
dugongs.A5.coda <- as.mcmc.list(dugongs.A5.sim)
summary(dugongs.A5.coda[,c("P.res[27]")])

## ---- Question A5b ----
bugs.boxplot(dugongs.A5.coda, "res")
density.strips(dugongs.A5.coda, "res")

## ---- Question A5c ----
 model.fit(dugongs_outl$Y, dugongs_outl$x, dugongs.A5.coda, "Y.pred",
          pch=20, lwd=2, xlab="Age (years)", ylab="Length (m)",
          ylim=c(1.5,3))

## ---- Question A5d ----
 model.fit(dugongs_outl$Y, dugongs_outl$x, dugongs.A1.coda, "Y.pred",
          pch=20, lwd=2, xlab="Age (years)", ylab="Length (m)",
          ylim=c(1.5,3))

## ---- Question A6 ----
dugongs_model_A6_file <- normalizePath("solutions/practical4/dugongs-model-A6.txt")
dugongs_in2_A6 <- list(alpha = 10, beta = 5, tau = 2, gamma = 0.2)
parameters.to.save <- c("alpha","beta","gamma","mu","sigma","res","P.res",
                        "P.pred","Y.pred")

dugongs.A6.sim <- bugs(model=dugongs_model_A6_file,
                       data=dugongs_outl,
                       inits=list(dugongs_in1, dugongs_in2_A6),
                       n.chains=2,
                       n.burnin=5000,
                       n.iter=55000,
                       n.thin=1,
                       parameters.to.save=parameters.to.save,
                       DIC=TRUE)
dugongs.A6.coda <- as.mcmc.list(dugongs.A6.sim)
dugongs.A6.sim[c("pD","DIC")]
densplot(dugongs.A6.coda[,"beta"])

## ---- Question B1 ----
parameters.to.save <- c("alpha","beta","p")

beetles.B1.sim <- bugs(model="beetles-model.txt",
                       data=beetles,
                       inits=list(beetles_in1,beetles_in2),
                       n.chains=2,
                       n.burnin=10000,
                       n.iter=20000,
                       n.thin=1,
                       parameters.to.save=parameters.to.save,
                       DIC=TRUE)
beetles.B1.coda <- as.mcmc.list(beetles.B1.sim)
beetles.B1.sim[c("pD","DIC")]

## ---- Question B2 ----
beetles_model_B3a_file <- normalizePath("solutions/practical3/beetles-model-B3a.txt")
beetles_in1_B3a <- list(alpha=0, beta=0)
beetles_in2_B3a <- list(alpha=1, beta=1)
parameters.to.save <- c("alpha","beta","p")

beetles.B3a.sim <- bugs(model=beetles_model_B3a_file,
                        data=beetles,
                        inits=list(beetles_in1_B3a,beetles_in2_B3a),
                        n.chains=2,
                        n.burnin=10000,
                        n.iter=25000,
                        n.thin=1,
                        parameters.to.save=parameters.to.save,
                        DIC=TRUE)
beetles.B3a.coda <- as.mcmc.list(beetles.B3a.sim)
beetles.B3a.sim[c("pD","DIC")]

## ---- Question B3 ----
beetles_model_B3b_file <- normalizePath("solutions/practical3/beetles-model-B3b.txt")
beetles_in1_B3b <- list(alpha=0, beta=0)
beetles_in2_B3b <- list(alpha=1, beta=1)
parameters.to.save <- c("alpha","beta","p")

beetles.B3b.sim <- bugs(model=beetles_model_B3b_file,
                        data=beetles,
                        inits=list(beetles_in1_B3b,beetles_in2_B3b),
                        n.chains=2,
                        n.burnin=10000,
                        n.iter=25000,
                        n.thin=1,
                        parameters.to.save=parameters.to.save,
                        DIC=TRUE)
beetles.B3b.coda <- as.mcmc.list(beetles.B3b.sim)
beetles.B3b.sim[c("pD","DIC")]

