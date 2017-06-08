library(rjags)

##################################################
### USEFUL PLOT UTILITY FUNCTIONS

## Get CODA samples of a vector-valued variable "varname".
## Rows are MCMC iterations, with multiple chains stacked
## on top of each other
## Columns are indices of the vector
## The functions below depend on this
get.coda.matrix <- function(coda, varname){
    var.cols <- grep(paste("^",varname,"\\[",sep=""), varnames(coda))
    m <- do.call("rbind", coda[,var.cols])
    indnames <- gsub(paste0(varname, "\\[([0-9]+)\\]"), "\\1", colnames(m))
    m[, order(as.numeric(indnames))]
}


## R implementation of the Win/OpenBUGS Inference->Compare->Model Fit
## function for plotting posterior fitted values from a regression model
model.fit <- function(y, # Data for outcome in regression (numeric)
                      x, # Data for predictor in regression (numeric)
                      coda, # CODA R object from BUGS/JAGS model fit
                      mu, # Variable name (character) for the
                          # fitted values (E(Y)) in the
                          # regression model
                      indices, # Indices of "mu" to use,
                      quantiles=c(0.025, 0.5, 0.975), # Quantiles of
                                                      # the fitted
                                                      # values to plot
                      ... # Other arguments to control the
                          # appearance of the plot, passed
                          # to the R functions plot()
                          # (plot.default()) and lines().
                          # See example below
                      ){
    res <- get.coda.matrix(coda, mu)
    res <- res[, indices, drop = FALSE]
    qs <- apply(res, 2, quantile, quantiles)
    lines(x[indices], qs[1,], col="lightblue", ...)
    lines(x[indices], qs[2,], col="red", ...)
    lines(x[indices], qs[3,], col="lightblue", ...)
}


##################################################
### QUESTION 1:  LONGITUDINAL DATA

teachers.data <- list(
  id = c(1, 1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5,
         6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11,
         11, 11, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15,
         16, 16, 17, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 20, 21,
         21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 25, 25, 25,
         25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 28, 29, 29, 29, 30, 30,
         31, 31, 31, 32, 32, 33, 33, 33, 34, 34, 34, 35, 35, 36, 36, 37,
         37, 37, 38, 38, 38, 39, 39, 39, 39, 40, 40, 41, 41, 41, 42, 42,
         42, 43, 43, 43, 44, 44, 44, 45, 46, 46, 46, 47, 48, 48, 48, 48,
         49, 49, 49, 50, 50, 50, 51, 51, 51),
  experience = c(1,2,4,1,1,3,2,3,4,1,2,3,4,1,2,3,4,3,4,1,2,3,4,1,2,3,4,1,2,3,4,
                 1,2,3,1,2,3,1,2,3,4,1,2,3,4,1,2,3,1,3,1,2,4,1,2,3,1,2,3,1,2,3,
                 4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,1,2,3,4,2,3,4,1,2,3,1,2,3,4,1,2,
                 3,1,4,2,3,4,1,2,1,2,4,1,2,4,1,3,1,2,1,2,4,1,3,4,1,2,3,4,1,4,1,
                 3,4,1,2,4,1,2,3,1,3,4,4,1,3,4,1,1,2,3,4,1,2,3,1,2,3,1,2,3),
  score = c(0.409999996423721, 1.04999995231628, 0.910000026226044,
            0.639999985694885, 1.12999999523163, 1.35000002384186,
            -0.400000005960464, 0.270000010728836, 0.259999990463257,
            1.01999998092651, 0.870000004768372, 0.980000019073486,
            0.970000028610229, 0.219999998807907, 0.28999999165535,
            0.829999983310699, 0.540000021457672, 0.330000013113022,
            0.680000007152557, 0.810000002384186, 0.769999980926514,
            0.879999995231628, 0.879999995231628, 0.319999992847443,
            0.829999983310699, 0.00999999977648258, 0.200000002980232,
            0.119999997317791, -0.129999995231628, 0.200000002980232,
            0.310000002384186, 1.11000001430511, 0.779999971389771,
            0.829999983310699, -0.620000004768372, -0.360000014305115,
            0.239999994635582, -0.46000000834465, 0.330000013113022,
            0.219999998807907, 0.430000007152557, 0.620000004768372,
            0.400000005960464, 0.180000007152557, 0.159999996423721,
            -0.519999980926514, 0.759999990463257, 0.819999992847443,
            -0.910000026226044, -0.100000001490116, 0.409999996423721,
            0.150000005960464, 0.219999998807907, -0.259999990463257, 0.5,
            0.430000007152557, 1.21000003814697, 0.479999989271164,
            0.620000004768372, 0.939999997615814, 1.00999999046326,
            0.730000019073486, 1, 0.300000011920929, 0.850000023841858,
            0.589999973773956, 0.740000009536743, 0.379999995231628,
            0.800000011920929, 0.800000011920929, 0.930000007152557,
            0.589999973773956, 0.800000011920929, 0.639999985694885,
            0.509999990463257, 0.699999988079071, 1.22000002861023,
            0.170000001788139, 0.819999992847443, 0.579999983310699,
            0.699999988079071, 1.10000002384186, 0.970000028610229,
            1.05999994277954, 1.4099999666214, 1.39999997615814,
            0.479999989271164, 0.899999976158142, 0.980000019073486,
            0.839999973773956, 0.790000021457672, 0.670000016689301,
            1, 0.980000019073486, 0.689999997615814, 0.300000011920929,
            1.01999998092651, 0.819999992847443, 1.07000005245209,
            1.21000003814697, 0.889999985694885, 1.16999995708466,
            1.14999997615814, 0.280000001192093, 1.03999996185303,
            0.850000023841858, 0.560000002384186, 1.22000002861023,
            1.13999998569489, 0.629999995231628, 0.990000009536743,
            0.990000009536743, 0.709999978542328, 0.990000009536743,
            0.509999990463257, 0.589999973773956, 0.330000013113022,
            -0.0199999995529652, 0.300000011920929, 0.629999995231628,
            0.109999999403954, 0.360000014305115, 0.280000001192093,
            0.370000004768372, 0.209999993443489, 0.0299999993294477,
            0.860000014305115, 0.769999980926514, 1.00999999046326,
            0.230000004172325, 0.439999997615814, 0.0599999986588955,
            1.01999998092651, 0.779999971389771, 1, 0.680000007152557,
            0.920000016689301, 0.970000028610229, 0.910000026226044,
            0.730000019073486, 0.629999995231628, 0.620000004768372,
            0.620000004768372, 0.550000011920929, 0.660000026226044,
            0.959999978542328, 0.939999997615814, 0.920000016689301,
            0.910000026226044, 0.920000016689301, 1.52999997138977,
            1.60000002384186, 1.44000005722046),
  # gender = c(1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
  #            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  #            0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  #            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  #            0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  #            0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
  #            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1,
  #            1, 1, 0, 0, 0),
  Nteachers = 51,
  Nmeas = 153)

# Question 1: standard linear regression
teachers.model <- "model {
  for (i in 1:Nmeas){
    score[i] ~ dnorm(mu[i], score.prec)

    # it might be worth centering the covariates, but it doesn't make much
    # difference here
    mu[i] <- beta1 + beta2 * experience[i]
    score.pred[i] ~ dnorm(mu[i], score.prec)
  }
  beta1 ~ dnorm(0, 0.0001)
  beta2 ~ dnorm(0, 0.0001)

  score.prec <- 1/pow(score.sd, 2)
  score.sd ~ dunif(0, 100)
}"

teachers.inits <-
  list(list(score.sd = 1,
            beta1 = 0,
            beta2 = 0,
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("beta1", "beta2", "score.sd", "mu", "score.pred")

teachers.jag <- jags.model(textConnection(teachers.model),
                                    data = teachers.data,
                                    inits = teachers.inits,
                                    n.chains = 1)
teachers.coda <- coda.samples(model = teachers.jag,
                                       variable.names = variable.names,
                                       n.iter = 50000)

summary(teachers.coda)

# Model fit for individual 13
id13 <- which(teachers.data$id == 13)
plot(x = teachers.data$experience[id13],
     y = teachers.data$score[id13],
     ylim = range(teachers.data$score))
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers.coda,
          "mu",
          indices = id13)

# Predictions for individual 13
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers.coda,
          "score.pred",
          indices = id13)

# Model fit for individual 16
id16 <- which(teachers.data$id == 16)
plot(x = teachers.data$experience[id16],
     y = teachers.data$score[id16],
     ylim = range(teachers.data$score))
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers.coda,
          "mu",
          indices = id16)

# Predictions for individual 16
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers.coda,
          "score.pred",
          indices = id16)


# Code to look at all data at the same time
plot(teachers.data$experience, teachers.data$score)
for (id in unique(teachers.data$id)){
  indices <- teachers.data$id == id
  lines(teachers.data$experience[indices], teachers.data$score[indices])
}
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers.coda,
          "mu",
          indices = teachers.data$id == 5)
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers.coda,
          "score.pred",
          indices = teachers.data$id == 5)

summary(teachers.coda[, c("beta1", "beta2")])

# Question 2: random intercept model

teachers.model2 <- "model{
for (i in 1:Nmeas){ 
    score[i] ~ dnorm(mu[i], score.prec)
    mu[i] <- beta1[id[i]] + beta2 * experience[i]
    score.pred[i] ~ dnorm(mu[i], score.prec)
  }
  for (j in 1:Nteachers){
  beta1[j] ~ dnorm(beta1.m,beta1.prec)
  }

beta1.m ~ dnorm(0, 0.0001)
beta2 ~ dnorm(0, 0.0001)

beta1.prec <- 1/pow(sd.beta1,2)
sd.beta1 ~ dunif(0, 100)



score.prec <- 1/pow(score.sd, 2)
score.sd ~ dunif(0, 100)


}"

teachers.inits2 <-
  list(list(score.sd = 1,
            sd.beta1=1,
            beta1.m=0,
            beta2 = 0,
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("beta1", "beta2", "score.sd","sd.beta1", "mu", "score.pred")

teachers.jag2 <- jags.model(textConnection(teachers.model2),
                                    data = teachers.data,
                                    inits = teachers.inits2,
                                    n.chains = 1)
teachers.coda2 <- coda.samples(model = teachers.jag2,
                                       variable.names = variable.names,
                                       n.iter = 50000)

summary(teachers.coda2)



# Question 3: random slope and intercept model

teachers.model3 <- "model{
for (i in 1:Nmeas){ 
score[i] ~ dnorm(mu[i], score.prec)
mu[i] <- beta1[id[i]] + beta2[id[i]] * experience[i]
score.pred[i] ~ dnorm(mu[i], score.prec)
}
for (j in 1:Nteachers){
beta1[j] ~ dnorm(beta1.m,beta1.prec)
beta2[j] ~ dnorm(beta2.m,beta2.prec)
}

beta1.m ~ dnorm(0, 0.0001)
beta2.m ~ dnorm(0, 0.0001)

beta1.prec <- 1/pow(sd.beta1,2)
sd.beta1 ~ dunif(0, 100)

beta2.prec <- 1/pow(sd.beta2,2)
sd.beta2 ~ dunif(0, 100)


score.prec <- 1/pow(score.sd, 2)
score.sd ~ dunif(0, 100)


}"

teachers.inits3 <-
  list(list(score.sd = 1,
            sd.beta1=1,
            sd.beta2=1,
            beta1.m=0,
            beta2.m=0,
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("beta1", "beta2", "score.sd","sd.beta1", "mu", "score.pred")

teachers.jag3<- jags.model(textConnection(teachers.model3),
                            data = teachers.data,
                            inits = teachers.inits3,
                            n.chains = 1)
teachers.coda3 <- coda.samples(model = teachers.jag3,
                               variable.names = variable.names,
                               n.iter = 50000)

summary(teachers.coda2)




# Question 4: Predictions for a new teacher

##################################################
### QUESTION 2:  MODEL CHECKING

election.data <- list(
  # Regions are (in order):
  # "Northern Ireland", "East", "East Midlands", "London", "North East",
  # "North West", "Scotland", "South East", "South West", "Wales",
  # "West Midlands", "Yorkshire and the Humber".
  y = round(c(0, 1445946, 969379, 1233386, 300883, 1050124, 434097, 2234360,
        1319994, 407813, 1097750, 796822)/1000),
  n = round(c(718103, 2948623, 2230402, 3536251, 1188153, 3364055, 2910465,
        4340668, 2836294, 1498063, 2628579, 2444177)/1000),
  N = 12
)

# Question 1: leaving out Northern Ireland (y[i]=1)
election.model <- "model {
  for (i in 2:N){
    y[i] ~ dbin(p[i], n[i])
  }
  for (i in 1:N){
    logit(p[i]) <- phi[i]
    phi[i] ~ dnorm(mu, tau)
  }
  mu ~ dunif(-10, 10)
  tau <- 1/pow(sd, 2)
  sd ~ dunif(0, 10)

  y.pred[1] ~ dbin(p[1], n[1])

  # mid p-value for unit 1 (Northern Ireland)
  P.value <- step(y.pred[1] - y[1] - 0.00001) +
                  0.5 * equals(y.pred[1], y[1])
}"

election.inits <-
  list(list(mu = 0,
            sd = 1,
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("mu", "tau", "p", "P.value")

election.jag <- jags.model(textConnection(election.model),
                           data = election.data,
                           inits = election.inits,
                           n.chains = 1)
election.coda <- coda.samples(model = election.jag,
                              variable.names = variable.names,
                              n.iter = 50000)
summary(election.coda[,"P.value"])

#            Mean        SD  Naive SE Time-series SE
# P.value  1.0000 0.0000000 0.000e+00      0.000e+00
#
# p-value is 1, suggesting very strong evidence that the model
# is not adequate for Northern Ireland (as expected!)

# Question 2: mixed approximate p-values

# Question 3(a): prior predictive

# Question 3(b): posterior predictive

# Question 4: using DUP results

##################################################
### QUESTION 3: JOINT MODELLING

install.packages("JM")  ### For the AIDS data
library(JM)

# (Start formatting data into a suitable data list)
install.packages("reshape2")
library(reshape2)

# Just use the first 20 patients
aids <- droplevels(subset(aids, patient %in% as.character(1:20)))

# clearer/shorter names
colnames(aids)[2] <- "surv.time"
colnames(aids)[5] <-  "t.obs"

aids.CD4.df <- dcast(patient ~ t.obs, value.var = "CD4", data = aids)[, -1]
aids.drug.df <- dcast(patient ~ t.obs, value.var = "drug", data = aids)[, 1:2]
aids.drug.ddI <- ifelse(aids.drug.df[, 2] == "ddI", 1, 0)

aids.t.obs <- sort(unique(aids$t.obs))
aids.t.pred <- seq(from = 0, to = max(aids.t.obs), by = 0.5)
aids.n.t.preds <- length(aids.t.pred)
aids.npredpatients <- 10

aids.npatients <- nrow(aids.CD4.df)
aids.n.t.obs <- length(aids.t.obs)

reps <- 1
aids.pred.cov <- expand.grid(drug = rep(c("ddC", "ddI"), each = reps))
aids.pred.drug.ddI <- ifelse(aids.pred.cov$drug == "ddI", 1, 0)
aids.npredpatients <- nrow(aids.pred.cov)

aids.death.df <- dcast(patient ~ t.obs, value.var = "death", data = aids)[, 1:2]
aids.is.cens <- ifelse(aids.death.df[, 2] == 0, 1, 0)

aids.time.df <- dcast(patient ~ t.obs, value.var = "surv.time", data = aids)[, 1:2]
aids.surv.time <- ifelse(aids.is.cens, NA, aids.time.df[, 2])
aids.break.time <- aids.time.df[, 2]

aids.data <- list(# Observations
                  I = aids.npatients,
                  CD4 = aids.CD4.df,
                  drug.ddI = aids.drug.ddI,
                  t.obs = aids.t.obs,
                  J = aids.n.t.obs,
                  is.cen = aids.is.cens,
                  surv.time = aids.surv.time,
                  break.time = aids.break.time,
                  # Predictions
                  K = aids.npredpatients,
                  drug.ddI.pred = aids.pred.drug.ddI,
                  t.pred = aids.t.pred,
                  L = aids.n.t.preds)
# (End formatting data into a suitable data list)

aids.model <- "
model {
  for (i in 1:I){
    for (j in 1:J){
      CD4[i, j] ~ dnorm(CD4.mu[i, j], CD4.sigma.prec)
      CD4.mu[i, j] <- beta[i, 1] +
                      beta[i, 2] * t.obs[j] +
                      beta3 * t.obs[j] * drug.ddI[i]
    }

    is.cen[i] ~ dinterval(surv.time[i], break.time[i])
    surv.time[i] ~ dweib(1, surv.mu[i])

    log(surv.mu[i]) <- beta4 +
                       beta5 * drug.ddI[i] +
                       r1 * beta[i, 1] +
                       r2 * beta[i, 2]

    beta[i, 1:2] ~ dmnorm(beta.mu[], beta.sigma.prec[,])
  }

  for (k in 1:K){
    for (l in 1:L){
      CD4.pred[k, l] ~ dnorm(CD4.pred.mu[k, l], CD4.sigma.prec)

      CD4.pred.mu[k, l] <- beta.pred[k, 1] +
                           beta.pred[k, 2] * t.pred[l] +
                           beta3 * t.pred[l] * drug.ddI.pred[k]
    }
    beta.pred[k, 1:2] ~ dmnorm(beta.mu[], beta.sigma.prec[,])
  }

  # priors
  beta.mu[1] ~ dnorm(0, 0.0001)
  beta.mu[2] ~ dnorm(0, 0.0001)
  beta.sigma.prec[1:2, 1:2] ~ dwish(R[,], d)
  beta.sigma.cov[1:2, 1:2] <- inverse(beta.sigma.prec[,])
  d <- 2
  R[1, 1] <- 1
  R[2, 1] <- 0
  R[1, 2] <- 0
  R[2, 2] <- 1

  beta3 ~ dnorm(0, 0.0001)
  beta4 ~ dnorm(0, 0.0001)
  beta5 ~ dnorm(0, 0.0001)
  CD4.sigma.prec <- 1/pow(CD4.sigma.sd, 2)
  CD4.sigma.sd ~ dunif(0, 100)

  r1 ~ dnorm(0, 0.0001)
  r2 ~ dnorm(0, 0.0001)
}
"
variable.names <- c("beta.mu", "beta.sigma.cov", "CD4.pred", "CD4.pred.mu",
                    "r1", "r2", "beta4", "beta5")

aids.inits <- list(list(
  surv.time = ifelse(is.na(aids.data$surv.time), aids.data$break.time + 1, NA),
  r1 = -0.2,
  r2 = -2,
  beta.mu = c(5, -0.2)
  ))

aids.jag <- jags.model(textConnection(aids.model),
                       data = aids.data,
                       inits = aids.inits,
                       n.chains = 1)
update(aids.jag, n.iter = 50000)
aids.coda <- coda.samples(model = aids.jag,
                          variable.names = variable.names,
                          n.iter = 100000,
                          thin = 30)

# Question 2: longitudinal only
