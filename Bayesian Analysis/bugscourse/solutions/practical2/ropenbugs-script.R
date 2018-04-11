## Solutions for practical2

## ---- Question A1 ----
parameters.to.save <- c("theta","r.pred","P.crit")
drug.A1.sim <- bugs(model="drug-model.txt",
                    data=drug,
                    inits=drug_in1,
                    n.chains=1,
                    n.burnin=0,
                    n.iter=30000,
                    n.thin=1,
                    parameters.to.save=parameters.to.save,
                    DIC=FALSE)
drug.A1.coda <- as.mcmc.list(drug.A1.sim)
traceplot(drug.A1.coda)

## ---- Question A1b ----
summary(drug.A1.coda)

## ---- Question A1c ----
par(mfrow=c(1,2))
densplot(drug.A1.coda[, "theta"])
densplot.discrete(drug.A1.coda[, "r.pred"])

## ---- Question A3 ----
# change a and b to 1 in the data
drug_data_A3 <- drug
drug_data_A3$a <- 1
drug_data_A3$b <- 1

drug.A3.sim <- bugs(model="drug-model.txt",
                    data=drug_data_A3,
                    inits=drug_in1,
                    n.chains=1,
                    n.burnin=0,
                    n.iter=30000,
                    n.thin=1,
                    parameters.to.save=c("theta","r.pred","P.crit"),
                    DIC=FALSE)
drug.A3.coda <- as.mcmc.list(drug.A3.sim)
summary(drug.A3.coda)

## ---- Question B1 ----
thm_model_B1_file <- normalizePath("solutions/practical2/thm-model-B1.txt")
thm_in1_B1 <- list(theta = 100)

# No data is needed here because it is included in the model
thm.B1.sim <- bugs(model=thm_model_B1_file,
                   data="nodata.txt",
                   inits=thm_in1_B1,
                   n.chains=1,
                   n.burnin=0,
                   n.iter=10000,
                   n.thin=1,
                   parameters.to.save=c("ypred","pmean.130","py.130","theta"),
                   DIC=FALSE)
thm.B1.coda <- as.mcmc.list(thm.B1.sim)
summary(thm.B1.coda)

## ---- Question B2a ----
thm_model_B2a_file <- normalizePath("solutions/practical2/thm-model-B2a.txt")
thm_in1_B1 <- list(theta = 100)

# No data is needed here because it is included in the model
thm.B2a.sim <- bugs(model=thm_model_B2a_file,
                   data="nodata.txt",
                   inits=thm_in1_B1,
                   n.chains=1,
                   n.burnin=0,
                   n.iter=10000,
                   n.thin=1,
                   parameters.to.save=c("ypred","pmean.130","py.130","theta"),
                   DIC=FALSE)
thm.B2a.coda <- as.mcmc.list(thm.B2a.sim)
summary(thm.B2a.coda)

## ---- Question B2b ----
thm_model_B2b_file <- normalizePath("solutions/practical2/thm-model-B2b.txt")
thm_in1_B1 <- list(theta = 100)

# No data is needed here because it is included in the model
thm.B2b.sim <- bugs(model=thm_model_B2b_file,
                   data="nodata.txt",
                   inits=thm_in1_B1,
                   n.chains=1,
                   n.burnin=0,
                   n.iter=10000,
                   n.thin=1,
                   parameters.to.save=c("ypred","pmean.130","py.130","theta"),
                   DIC=FALSE)
thm.B2b.coda <- as.mcmc.list(thm.B2b.sim)
summary(thm.B2b.coda)

## ---- Question B3 ----
par(mfrow=c(1,2))
densplot(thm.B1.coda[,c("theta","ypred")])
densplot(thm.B2a.coda[,c("theta","ypred")])
densplot(thm.B2b.coda[,c("theta","ypred")])

