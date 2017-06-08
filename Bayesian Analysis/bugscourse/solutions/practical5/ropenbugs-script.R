## Solutions for practical5

## ---- Question A1 ----
parameters.to.save <- c("mu","sigma2","psi2","vpc","theta")

thm.A1.sim <- bugs(model="thm-model.txt",
                   data=thm,
                   inits=list(thm_in1,thm_in2),
                   n.chains=2,
                   n.burnin=5000,
                   n.iter=15000,
                   n.thin=1,
                   parameters.to.save=parameters.to.save,
                   DIC=TRUE)
thm.A1.coda <- as.mcmc.list(thm.A1.sim)
thm.A1.sim[c("pD","DIC")]
summary(thm.A1.coda[,c("mu","sigma2","psi2","vpc")])

## ---- Question A1b ----
bugs.boxplot(thm.A1.coda, "theta",ordered=TRUE)

## ---- Question A2 ----
thm_model_A2_file <- normalizePath("solutions/practical5/thm-model-A2.txt")
thm_in1_A2 <- list(mu = 50, psi = 1, nu = 2, phi = 0.5)
thm_in2_A2 <- list(mu = 200, psi = 10, nu = 1, phi = 0.25)
parameters.to.save <- c("phi2","mean.var","mu","nu",
                        "psi2","sigma2","vpc","theta")

thm.A2.sim <- bugs(model=thm_model_A2_file,
                   data=thm,
                   inits=list(thm_in1_A2,thm_in2_A2),
                   n.chains=2,
                   n.burnin=5000,
                   n.iter=15000,
                   n.thin=100, # retain only every 100th iteration
                   parameters.to.save=parameters.to.save,
                   DIC=TRUE)
thm.A2.coda <- as.mcmc.list(thm.A2.sim)
thm.A2.sim[c("pD","DIC")]
summary(thm.A2.coda[,c("phi2","mean.var","mu","nu",
                       "psi2","vpc")])

## ---- Question A2b ----
bugs.boxplot(thm.A2.coda, "theta",ordered=TRUE)

## ---- Question A2c ----
bugs.boxplot(thm.A2.coda, "sigma2",ordered=TRUE)

## ---- Question A4 ----
thm_model_A4_file <- normalizePath("solutions/practical5/thm-x-model-A4.txt")
thm_in1_A4 <- list(mu = 50, psi = 1, alpha = 2, beta = 0.5)
thm_in2_A4 <- list(mu = 200, omega = 10, alpha = 1, beta = 0.1)
thm_x_all <- c(thm_x, thm) # append thm_x data to the original thm data
parameters.to.save <- c("beta","mean.var","mu","psi2","vpc","theta")

thm.A4.sim <- bugs(model=thm_model_A4_file,
                   data=thm_x_all,
                   inits=list(thm_in1_A4,thm_in2_A4),
                   n.chains=2,
                   n.burnin=10000,
                   n.iter=40000,
                   n.thin=1,
                   parameters.to.save=parameters.to.save,
                   DIC=TRUE)
thm.A4.coda <- as.mcmc.list(thm.A4.sim)
thm.A4.sim[c("pD","DIC")]
summary(thm.A4.coda[,c("beta","mean.var","mu","psi2","vpc")])

## ---- Question A4b ----
bugs.boxplot(thm.A4.coda, "theta",ordered=TRUE)

