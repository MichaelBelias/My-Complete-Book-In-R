## Solutions for practical3

## ---- Question A1 ----
dugongs_in2_A1 <- list(alpha = 10, beta = 5, tau = 2, gamma = 0.2)
parameters.to.save <- c("alpha","beta","gamma","mu","sigma")

dugongs.A1.sim <- bugs(model="dugongs-model.txt",
                       data=dugongs,
                       inits=list(dugongs_in1, dugongs_in2_A1),
                       n.chains=2, # note 2 chains not 1 as before
                       n.burnin=10000,
                       n.iter=40000, # draw 40000, throw away first 10000
                       n.thin=1,
                       parameters.to.save=parameters.to.save,
                       DIC=FALSE)
dugongs.A1.coda <- as.mcmc.list(dugongs.A1.sim)
summary(dugongs.A1.coda[,c("alpha","beta","gamma","sigma")])

## ---- Question A1b ----
model.fit(dugongs$Y, dugongs$x, dugongs.A1.coda, "mu",
          pch=20, lwd=2, xlab="Age (years)", ylab="Length (m)",
          ylim=c(1.5,3))

## ---- Question A2 ----
dugongs_model_A2_file <- normalizePath("solutions/practical3/dugongs-model-A2.txt")
dugongs_in1_A2 <- list(alpha = 1, beta = 1, log.sigma = 1, gamma = 0.9)
dugongs_in2_A2 <- list(alpha = 1, beta = 1, log.sigma = 0.1, gamma = 0.9)
parameters.to.save <- c("alpha","beta","gamma","mu","sigma")

dugongs.A2.sim <- bugs(model=dugongs_model_A2_file,
                       data=dugongs,
                       inits=list(dugongs_in1_A2, dugongs_in2_A2),
                       n.chains=2,
                       n.burnin=10000,
                       n.iter=40000,
                       n.thin=1,
                       parameters.to.save=parameters.to.save,
                       DIC=FALSE)
dugongs.A2.coda <- as.mcmc.list(dugongs.A2.sim)
summary(dugongs.A2.coda[,c("alpha","beta","gamma","sigma")])

## ---- Question A3 ----
densplot(dugongs.A1.coda[,"gamma"])
dugongs_in2_A3 <- list(alpha = 10, beta = 5, tau = 2, gamma = 0.2)
parameters.to.save <- c("alpha","beta","gamma","mu","sigma")

dugongs.A3.sim <- bugs(model="dugongs-model.txt",
                   data=dugongs,
                   inits=list(dugongs_in1, dugongs_in2_A3),
                   n.chains=2, # note 2 chains not 1 as before
                   n.burnin=10000,
                   n.iter=40000,
                   n.thin=1,
                   parameters.to.save=parameters.to.save,
                   DIC=FALSE)
dugongs.A3.coda <- as.mcmc.list(dugongs.A3.sim)

## ---- Question A3b ----
densplot(dugongs.A3.coda[,"gamma"])

## ---- Question B1 ----
parameters.to.save <- c("alpha","beta","p")

beetles.B1.sim <- bugs(model="beetles-model.txt",
                   data=beetles,
                   inits=list(beetles_in1,beetles_in2),
                   n.chains=2,
                   n.burnin=0, # start monitoring immediate to
                               # see good convergence of centred version
                   n.iter=40000,
                   n.thin=1,
                   parameters.to.save=parameters.to.save,
                   DIC=FALSE)
beetles.B1.coda <- as.mcmc.list(beetles.B1.sim)
model.fit(x = beetles$x, y = beetles$r/beetles$n, beetles.B1.coda, "p",
          pch=20, lwd=2, xlab="Dose", ylab="Proportion dying",
          xlim = c(1.65, 1.9), ylim = c(0, 1))

traceplot(beetles.B1.coda[,c("alpha","beta")])
densplot(beetles.B1.coda[,c("alpha","beta")])

## ---- Question B2 ----
beetles_model_B2_file <- normalizePath("solutions/practical3/beetles-model-B2.txt")
parameters.to.save <- c("alpha","beta","p")

beetles.B2.sim <- bugs(model=beetles_model_B2_file,
                   data=beetles,
                   inits=list(beetles_in1,beetles_in2),
                   n.chains=2,
                   n.burnin=0, # start monitoring immediate to
                               # see poor convergence of uncentred version
                   n.iter=40000,
                   n.thin=1,
                   parameters.to.save=parameters.to.save,
                   DIC=FALSE)
beetles.B2.coda <- as.mcmc.list(beetles.B2.sim)
traceplot(beetles.B2.coda[,c("alpha","beta")])
densplot(beetles.B2.coda[,c("alpha","beta")])

## ---- Question B3a ----
beetles_model_B3a_file <- normalizePath("solutions/practical3/beetles-model-B3a.txt")
beetles_in1_B3a <- list(alpha=0, beta=0)
beetles_in2_B3a <- list(alpha=1, beta=1)
parameters.to.save <- c("alpha","beta","p")

beetles.B3a.sim <- bugs(model=beetles_model_B3a_file,
                        data=beetles,
                        inits=list(beetles_in1_B3a, beetles_in2_B3a),
                        n.chains=2,
                        n.burnin=1000,
                        n.iter=11000,
                        n.thin=1,
                        parameters.to.save=parameters.to.save,
                        DIC=FALSE)
beetles.B3a.coda <- as.mcmc.list(beetles.B3a.sim)
summary(beetles.B3a.coda[,c("alpha","beta")])

## ---- Question B3b ----
beetles_model_B3b_file <- normalizePath("solutions/practical3/beetles-model-B3b.txt")
beetles_in1_B3b <- list(alpha=0, beta=30)
beetles_in2_B3b <- list(alpha=1, beta=0)
parameters.to.save <- c("alpha","beta","p")

beetles.B3b.sim <- bugs(model=beetles_model_B3b_file,
                        data=beetles,
                        inits=list(beetles_in1_B3b, beetles_in2_B3b),
                        n.chains=2,
                        n.burnin=1000,
                        n.iter=11000,
                        n.thin=1,
                        parameters.to.save=parameters.to.save,
                        DIC=FALSE)
beetles.B3b.coda <- as.mcmc.list(beetles.B3b.sim)
summary(beetles.B3b.coda[,c("alpha","beta")])

## ---- Question B3c ----
beetles_model_B3c_file <- normalizePath("solutions/practical3/beetles-model-B3c.txt")
beetles_in1_B3b <- list(alpha=0, beta=30)
beetles_in2_B3b <- list(alpha=1, beta=0)
parameters.to.save <- c("alpha","beta","p")

beetles.B3c.sim <- bugs(model=beetles_model_B3c_file,
                        data=beetles,
                        inits=list(beetles_in1_B3b, beetles_in2_B3b),
                        n.chains=2,
                        n.burnin=1000,
                        n.iter=11000,
                        n.thin=1,
                        parameters.to.save=parameters.to.save,
                        DIC=FALSE)
beetles.B3c.coda <- as.mcmc.list(beetles.B3c.sim)
summary(beetles.B3c.coda[,c("alpha","beta")])

