## Solutions for practical3

## ---- Question A1 ----
# create a second set of initial values
dugongs_in2_A1 <- list(alpha = 10, beta = 5, tau = 2, gamma = 0.2)

dugongs.A1.jag <- jags.model(file="dugongs-model.txt",
                             data=dugongs,
                             inits=list(dugongs_in1, dugongs_in2_A1),
                             n.chains=2,
                             quiet=TRUE)

# Run 10000 burn-in samples before monitoring
update(dugongs.A1.jag, n.iter = 10000)

variable.names <- c("alpha","beta","gamma","mu","sigma")
dugongs.A1.coda <- coda.samples(dugongs.A1.jag,
                                variable.names=variable.names,
                                n.iter=30000)
summary(dugongs.A1.coda[,c("alpha","beta","gamma","sigma")])

## ---- Question A1b ----
model.fit(dugongs$Y, dugongs$x, dugongs.A1.coda, "mu",
          pch=20, lwd=2, xlab="Age (years)", ylab="Length (m)",
          ylim=c(1.5,3))

## ---- Question A2 ----
dugongs_model_A2_file <- "solutions/practical3/dugongs-model-A2.txt"
dugongs_in1_A2 <- list(alpha = 1, beta = 1, log.sigma = 1, gamma = 0.9)
dugongs_in2_A2 <- list(alpha = 1, beta = 1, log.sigma = 0.1, gamma = 0.9)

dugongs.A2.jag <- jags.model(file=dugongs_model_A2_file,
                             data=dugongs,
                             inits=list(dugongs_in1_A2, dugongs_in2_A2),
                             n.chains=2,
                             quiet=TRUE)

# Run 10000 burn-in samples before monitoring
update(dugongs.A2.jag, n.iter = 10000)

variable.names <- c("alpha","beta","gamma","mu","sigma")
dugongs.A2.coda <- coda.samples(dugongs.A2.jag,
                                variable.names=variable.names,
                                n.iter=30000)
summary(dugongs.A2.coda[,c("alpha","beta","gamma","sigma")])

## ---- Question A3 ----
densplot(dugongs.A1.coda[,"gamma"])
# create a second set of initial values
dugongs_in2_A3 <- list(alpha = 10, beta = 5, tau = 2, gamma = 0.2)

dugongs.A3.jag <- jags.model(file="dugongs-model.txt",
                             data=dugongs,
                             inits=list(dugongs_in1, dugongs_in2_A3),
                             n.chains=2,
                             quiet=TRUE)

# Run 10000 burn-in samples before monitoring
update(dugongs.A3.jag, n.iter = 10000)

variable.names <- c("alpha","beta","gamma","mu","sigma")
dugongs.A3.coda <- coda.samples(dugongs.A3.jag,
                                variable.names=variable.names,
                                n.iter=30000)

## ---- Question A3b ----
densplot(dugongs.A3.coda[,"gamma"])

## ---- Question B1 ----
# Use beetles_in1_jags and beetles_in2_jags, which are
# slightly less extreme than the Win/OpenBUGS versions
# otherwise JAGS will not run
beetles.B1.jag <- jags.model(file="beetles-model.txt",
                             data=beetles,
                             inits=list(beetles_in1_jags,
                                        beetles_in2_jags),
                             n.chains=2,
                             quiet=TRUE)
# Start monitoring immediately to see
# good convergence of centred version
# update(beetles.B1.jag, n.iter = 10000)
variable.names <- c("alpha","beta","p")
beetles.B1.coda <- coda.samples(beetles.B1.jag,
                                variable.names=variable.names,
                                n.iter=40000)
model.fit(x = beetles$x, y = beetles$r/beetles$n, beetles.B1.coda, "p",
          pch=20, lwd=2, xlab="Dose", ylab="Proportion dying",
          xlim = c(1.65, 1.9), ylim = c(0, 1))

traceplot(beetles.B1.coda[,c("alpha","beta")])
densplot(beetles.B1.coda[,c("alpha","beta")])

## ---- Question B2 ----
beetles_model_B2_file <- "solutions/practical3/beetles-model-B2.txt"

# Use beetles_in1_jags and beetles_in2_jags, which are
# slightly less extreme than the Win/OpenBUGS versions
# otherwise JAGS will not run
beetles.B2.jag <- jags.model(file=beetles_model_B2_file,
                             data=beetles,
                             inits=list(beetles_in1_jags,
                                        beetles_in2_jags),
                             n.chains=2,
                             quiet=TRUE)
# Start monitoring immediately to see poor convergence of uncentred version
# update(beetles.B2.jag, n.iter = 10000)
variable.names <- c("alpha","beta","p")
beetles.B2.coda <- coda.samples(beetles.B2.jag,
                                variable.names=variable.names,
                                n.iter=40000)
traceplot(beetles.B2.coda[,c("alpha","beta")])
densplot(beetles.B2.coda[,c("alpha","beta")])

## ---- Question B3a ----
beetles_model_B3a_file <- "solutions/practical3/beetles-model-B3a.txt"
beetles_in1_B3a <- list(alpha=0, beta=0)
beetles_in2_B3a <- list(alpha=1, beta=1)

beetles.B3a.jag <- jags.model(file=beetles_model_B3a_file,
                              data=beetles,
                              inits=list(beetles_in1_B3a,beetles_in2_B3a),
                              n.chains=2,
                              quiet=TRUE)
update(beetles.B3a.jag, n.iter = 1000)
variable.names <- c("alpha","beta","p")
beetles.B3a.coda <- coda.samples(beetles.B3a.jag,
                                 variable.names=variable.names,
                                 n.iter=10000)
summary(beetles.B3a.coda[,c("alpha","beta")])

## ---- Question B3b ----
beetles_model_B3b_file <- "solutions/practical3/beetles-model-B3b.txt"
beetles_in1_B3b <- list(alpha=0, beta=30)
beetles_in2_B3b <- list(alpha=1, beta=0)

beetles.B3b.jag <- jags.model(file=beetles_model_B3b_file,
                              data=beetles,
                              inits=list(beetles_in1_B3b,beetles_in2_B3b),
                              n.chains=2,
                              quiet=TRUE)
update(beetles.B3b.jag, n.iter = 1000)
variable.names <- c("alpha","beta","p")
beetles.B3b.coda <- coda.samples(beetles.B3b.jag,
                                 variable.names=variable.names,
                                 n.iter=10000)
summary(beetles.B3b.coda[,c("alpha","beta")])

## ---- Question B3c ----
beetles_model_B3c_file <- "solutions/practical3/beetles-model-B3c.txt"
beetles_in1_B3b <- list(alpha=0, beta=30)
beetles_in2_B3b <- list(alpha=1, beta=0)

beetles.B3c.jag <- jags.model(file=beetles_model_B3c_file,
                              data=beetles,
                              inits=list(beetles_in1_B3b,beetles_in2_B3b),
                              n.chains=2,
                              quiet=TRUE)
update(beetles.B3c.jag, n.iter = 1000)
variable.names <- c("alpha","beta","p")
beetles.B3c.coda <- coda.samples(beetles.B3c.jag,
                                 variable.names=variable.names,
                                 n.iter=10000)
summary(beetles.B3c.coda[,c("alpha","beta")])

