## Solutions for practical2

## ---- Question A1 ----
drug.A1.jag <- jags.model(file="drug-model.txt",
                          data=drug,
                          inits=drug_in1,
                          quiet=TRUE)
variable.names <- c("theta","r.pred","P.crit")
drug.A1.coda <- coda.samples(drug.A1.jag,
                             variable.names=variable.names,
                             n.iter=30000)
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

drug.A3.jag <- jags.model(file="drug-model.txt",
                          data=drug_data_A3,
                          inits=drug_in1,
                          quiet=TRUE)
drug.A3.coda <- coda.samples(drug.A3.jag,
                             variable.names=c("theta","r.pred","P.crit"),
                             n.iter=30000)
summary(drug.A3.coda)

## ---- Question B1 ----
thm_model_B1_file <- "solutions/practical2/thm-model-B1-jags.txt"
thm_data_B1 <- list(y = rep(130, 2))
thm_in1_B1 <- list(theta = 100)

thm.B1.jag <- jags.model(file=thm_model_B1_file,
                         data=thm_data_B1,
                         inits=thm_in1_B1,
                         quiet=TRUE)
thm.B1.coda <- coda.samples(thm.B1.jag,
                            variable.names=c("ypred","pmean.130",
                                             "py.130","theta"),
                            n.iter=10000)
summary(thm.B1.coda)

## ---- Question B2a ----
thm_model_B2a_file <- "solutions/practical2/thm-model-B2a-jags.txt"
thm_data_B1 <- list(y = rep(130, 2))
thm_in1_B1 <- list(theta = 100)

thm.B2a.jag <- jags.model(file=thm_model_B2a_file,
                         data=thm_data_B1,
                         inits=thm_in1_B1,
                         quiet=TRUE)
thm.B2a.coda <- coda.samples(thm.B2a.jag,
                            variable.names=c("ypred","pmean.130",
                                             "py.130","theta"),
                            n.iter=10000)
summary(thm.B2a.coda)

## ---- Question B2b ----
thm_model_B2b_file <- "solutions/practical2/thm-model-B2b-jags.txt"
thm_data_B1 <- list(y = rep(130, 2))
thm_in1_B1 <- list(theta = 100)

thm.B2b.jag <- jags.model(file=thm_model_B2b_file,
                         data=thm_data_B1,
                         inits=thm_in1_B1,
                         quiet=TRUE)
thm.B2b.coda <- coda.samples(thm.B2b.jag,
                            variable.names=c("ypred","pmean.130",
                                             "py.130","theta"),
                            n.iter=10000)
summary(thm.B2b.coda)

## ---- Question B3 ----
par(mfrow=c(1,2))
densplot(thm.B1.coda[,c("theta","ypred")])
densplot(thm.B2a.coda[,c("theta","ypred")])
densplot(thm.B2b.coda[,c("theta","ypred")])

