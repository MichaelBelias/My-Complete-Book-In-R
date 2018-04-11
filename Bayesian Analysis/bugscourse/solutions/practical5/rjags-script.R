## Solutions for practical5

## ---- Question A1 ----
thm.A1.jag <- jags.model(file="thm-model.txt",
                         data=thm,
                         inits=list(thm_in1,thm_in2),
                         n.chains=2,
                         quiet=TRUE)
update(thm.A1.jag, n.iter = 5000)
variable.names <- c("mu","sigma2","psi2","vpc","theta")
thm.A1.coda <- coda.samples(thm.A1.jag,
                            variable.names=variable.names,
                            n.iter=10000)
thm.A1.pD <- dic.samples(thm.A1.jag, 10000, type="pD")
thm.A1.pD
summary(thm.A1.coda[,c("mu","sigma2","psi2","vpc")])

## ---- Question A1b ----
bugs.boxplot(thm.A1.coda, "theta",ordered=TRUE)

## ---- Question A2 ----
thm_model_A2_file <- "solutions/practical5/thm-model-A2.txt"
thm_in1_A2 <- list(mu = 50, psi = 1, nu = 2, phi = 0.5)
thm_in2_A2 <- list(mu = 200, psi = 10, nu = 1, phi = 0.25)
thm.A2.jag <- jags.model(file=thm_model_A2_file,
                         data=thm,
                         inits=list(thm_in1_A2,thm_in2_A2),
                         n.chains=2,
                         quiet=TRUE)
update(thm.A2.jag, n.iter = 5000)
variable.names <- c("phi2","mean.var","mu","nu","psi2",
                    "sigma2","vpc","theta")
thm.A2.coda <- coda.samples(thm.A2.jag,
                            variable.names=variable.names,
                            n.iter=10000,
                            thin=100 # retain only every 100th iteration
                            )
thm.A2.pD <- dic.samples(thm.A2.jag, 10000, thin = 100, type="pD")
thm.A2.pD
summary(thm.A2.coda[,c("phi2","mean.var","mu","nu",
                       "psi2","vpc")])

## ---- Question A2b ----
bugs.boxplot(thm.A2.coda, "theta",ordered=TRUE)

## ---- Question A2c ----
bugs.boxplot(thm.A2.coda, "sigma2",ordered=TRUE)

## ---- Question A4 ----
thm_model_A4_file <- "solutions/practical5/thm-x-model-A4.txt"
thm_in1_A4 <- list(mu = 50, psi = 1, alpha = 2, beta = 0.5)
thm_in2_A4 <- list(mu = 200, psi = 10, alpha = 1, beta = 0.1)
thm_x_all <- c(thm_x, thm) # append thm_x data to the original thm data

thm.A4.jag <- jags.model(file=thm_model_A4_file,
                         data=thm_x_all,
                         inits=list(thm_in1_A4,thm_in2_A4),
                         n.chains=2,
                         quiet=TRUE)
update(thm.A4.jag, n.iter = 10000)
variable.names <- c("beta","mean.var","mu","psi2","vpc","theta")
thm.A4.coda <- coda.samples(thm.A4.jag,
                            variable.names=variable.names,
                            n.iter=30000)
thm.A4.pD <- dic.samples(thm.A4.jag, 30000, type="pD")
thm.A4.pD
summary(thm.A4.coda[,c("beta","mean.var","mu","psi2","vpc")])

## ---- Question A4b ----
bugs.boxplot(thm.A4.coda, "theta",ordered=TRUE)

