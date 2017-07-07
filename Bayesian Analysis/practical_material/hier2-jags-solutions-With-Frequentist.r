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


##### Create a data.frame of the data
IPD  = as.data.frame(teachers.data[which(lapply(teachers.data, length)>1)])

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

### 

summary(teachers.coda[, c("beta1", "beta2")])

fit = lm(score ~ experience , data= IPD)


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
teachers.model2 <- "model {
for (i in 1:Nmeas){
score[i] ~ dnorm(mu[i], score.prec)
mu[i] <- beta1[id[i]] + beta2 * experience[i]
score.pred[i] ~ dnorm(mu[i], score.prec)
}
for (j in 1:Nteachers){
beta1[j] ~ dnorm(beta1.mean, beta1.prec)
}
beta1.mean ~ dnorm(0, 0.0001)
beta2 ~ dnorm(0, 0.0001)

beta1.prec <- 1/pow(beta1.sd, 2)
beta1.sd ~ dunif(0, 100)

score.prec <- 1/pow(score.sd, 2)
score.sd ~ dunif(0, 100)
}"

teachers2.inits <-
  list(list(score.sd = 1,
            beta1.sd = 1,
            beta2 = 0,
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("beta1", "beta1.mean", "beta1.sd", "beta2", "score.sd",
                    "mu", "score.pred")

teachers2.jag <- jags.model(textConnection(teachers.model2),
                            data = teachers.data,
                            inits = teachers2.inits,
                            n.chains = 1)
teachers2.coda <- coda.samples(model = teachers2.jag,
                               variable.names = variable.names,
                               n.iter = 50000)


id13 <- which(teachers.data$id == 13)
plot(x = teachers.data$experience[id13],
     y = teachers.data$score[id13],
     ylim = range(teachers.data$score))
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers2.coda,
          "mu",
          indices = id13)
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers2.coda,
          "score.pred",
          indices = id13)

id16 <- which(teachers.data$id == 16)
plot(x = teachers.data$experience[id16],
     y = teachers.data$score[id16],
     ylim = range(teachers.data$score))
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers2.coda,
          "mu",
          indices = id16)
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers2.coda,
          "score.pred",
          indices = id16)

summary(teachers2.coda[,c("beta1.mean", "beta1.sd", "beta2")])

# beta1.mean       0.4639296  0.5597252  0.60948  0.6583536 0.75229
# beta1.sd         0.2781638  0.3258148  0.35466  0.3864180 0.45401
# beta2           -0.0246069  0.0027522  0.01681  0.0306800 0.05795

# The model fits show much closer agreement between the observations and
# the model. This is not really suprising, since there is clear heterogeneity
# in the teachers' ratings at experience=0.
#
# The score.pred predictive check is a "posterior predictive" method, and so
# will generally be rather conservative: we are using the posterior distribution
# of mu[i] to make these predictions and this posterior distribution
# "has access to" data point score[i]. So the prediction is based upon a
# quantity (mu[i]) that is very closely based upon the observed data (score[i]).
# The prediction score.pred[i] could be viewed as "cheating" in a sense, and
# so will tend to suggest the model is better than a "fairer" check might
# suggest.

# Question 3: random slope and intercept model
teachers3.model <- "model {
for (i in 1:Nmeas){
score[i] ~ dnorm(mu[i], score.prec)
mu[i] <- beta[id[i], 1] + beta[id[i], 2] * experience[i]
score.pred[i] ~ dnorm(mu[i], score.prec)
}
for (j in 1:Nteachers){
beta[j, 1:2] ~ dmnorm(beta.mean[], beta.prec[,])
}

beta.mean[1] ~ dnorm(0, 0.0001)
beta.mean[2] ~ dnorm(0, 0.0001)
beta.prec[1:2, 1:2] ~ dwish(R[,], d)
beta.cov[1:2, 1:2] <- inverse(beta.prec[1:2, 1:2])
d <- 2
R[1, 1] <- 1
R[2, 1] <- 0
R[1, 2] <- 0
R[2, 2] <- 1

score.prec <- 1/pow(score.sd, 2)
score.sd ~ dunif(0, 100)
}"

teachers3.inits <-
  list(list(score.sd = 1,
            beta.mean = c(0, 0),
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("beta", "beta.mean", "beta.cov", "score.sd", "mu",
                    "score.pred")

teachers3.jag <- jags.model(textConnection(teachers3.model),
                            data = teachers.data,
                            inits = teachers3.inits,
                            n.chains = 1)
teachers3.coda <- coda.samples(model = teachers3.jag,
                               variable.names = variable.names,
                               n.iter = 50000)

summary(teachers3.coda[,"score.pred[13]"])

# Question 4: Predictions for a new teacher
teachers4.model <- "model {
for (i in 1:Nmeas){
score[i] ~ dnorm(mu[i], score.prec)
mu[i] <- beta[id[i], 1] + beta[id[i], 2] * experience[i]
score.pred[i] ~ dnorm(mu[i], score.prec)
}
for (j in 1:Nteachers){
beta[j, 1:2] ~ dmnorm(beta.mean[], beta.prec[,])
}

beta.mean[1] ~ dnorm(0, 0.0001)
beta.mean[2] ~ dnorm(0, 0.0001)
beta.prec[1:2, 1:2] ~ dwish(R[,], d)
beta.cov[1:2, 1:2] <- inverse(beta.prec[1:2, 1:2])
d <- 2
R[1, 1] <- 1
R[2, 1] <- 0
R[1, 2] <- 0
R[2, 2] <- 1

score.prec <- 1/pow(score.sd, 2)
score.sd ~ dunif(0, 100)

# Predictions
beta.new[1:2] ~ dmnorm(beta.mean[], beta.prec[,])
for (i in 1:4){
score.new[i] ~ dnorm(mu.new[i], score.prec)
# mu.new[i] = mean after i-1 years of experience
mu.new[i] <- beta.new[1] + beta.new[2] * (i - 1)
}
}"

teachers4.inits <-
  list(list(score.sd = 1,
            beta.mean = c(0, 0),
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("beta", "beta.mean", "beta.cov", "score.sd", "mu",
                    "beta.new", "mu.new", "score.pred")

teachers4.jag <- jags.model(textConnection(teachers4.model),
                            data = teachers.data,
                            inits = teachers4.inits,
                            n.chains = 1)
teachers4.coda <- coda.samples(model = teachers4.jag,
                               variable.names = variable.names,
                               n.iter = 50000)
summary(teachers4.coda[,"score.pred[13]"])

id13 <- which(teachers.data$id == 13)
plot(x = teachers.data$experience[id13],
     y = teachers.data$score[id13],
     ylim = range(teachers.data$score))
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers4.coda,
          "mu",
          indices = id13)

model.fit(teachers.data$score,
          teachers.data$experience,
          teachers4.coda,
          "score.pred",
          indices = id13)

id16 <- which(teachers.data$id == 16)
plot(x = teachers.data$experience[id16],
     y = teachers.data$score[id16],
     ylim = range(teachers.data$score))
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers4.coda,
          "mu",
          indices = id16)
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers4.coda,
          "score.pred",
          indices = id16)

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


##### Create a data.frame of the data
IPD  = as.data.frame(election.data[which(lapply(election.data, length)>1)])

# Question 1: leaving out Northern Ireland
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
election.model <- "model {
for (i in 1:N){
y[i] ~ dbin(p[i], n[i])
logit(p[i]) <- phi[i]
phi[i] ~ dnorm(mu, tau)
}
mu ~ dunif(-10, 10)
tau <- 1/pow(sd, 2)
sd ~ dunif(0, 10)

for (i in 1:N){
phi.pred[i] ~ dnorm(mu, tau)
logit(p.pred[i]) <- phi.pred[i]
y.pred[i] ~ dbin(p.pred[i], n[i])
P.mixed[i] <- step(y.pred[i] - y[i] - 0.00001) +
0.5 * equals(y.pred[i], y[i])
}
}"

election.inits <-
  list(list(mu = 0,
            sd = 1,
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("mu", "tau", "p", "P.mixed", "p.pred")

election.jag <- jags.model(textConnection(election.model),
                           data = election.data,
                           inits = election.inits,
                           n.chains = 1)
election.coda <- coda.samples(model = election.jag,
                              variable.names = variable.names,
                              n.iter = 50000)
summary(election.coda)

#                  Mean       SD  Naive SE Time-series SE
# P.mixed[1]   0.994350 0.052850 2.364e-04      2.664e-04
# P.mixed[2]   0.300420 0.458402 2.050e-03      2.095e-03
# P.mixed[3]   0.337680 0.472850 2.115e-03      2.143e-03
# P.mixed[4]   0.405950 0.491024 2.196e-03      2.196e-03
# P.mixed[5]   0.501330 0.499808 2.235e-03      2.292e-03
# P.mixed[6]   0.442050 0.496560 2.221e-03      2.247e-03
# P.mixed[7]   0.631930 0.482114 2.156e-03      2.238e-03
# P.mixed[8]   0.284890 0.451350 2.018e-03      2.054e-03
# P.mixed[9]   0.314020 0.464064 2.075e-03      2.116e-03
# P.mixed[10]  0.479540 0.499386 2.233e-03      2.288e-03
# P.mixed[11]  0.354630 0.478327 2.139e-03      2.174e-03
# P.mixed[12]  0.428620 0.494803 2.213e-03      2.245e-03

# The p-value for Northern Ireland is about 0.99, so there is slight
# attenuation of the p-value. The Northern Ireland data are a relatively
# large proportion of the data in this example, and it is pulling the
# entire posterior distribution towards it, as any data point would do.

# Question 3(a): prior predictive
election.prior.model <- "model {
# Remove the data from the model so that we just sample from the prior
# for (i in 1:N){
#   y[i] ~ dbin(p[i], n[i])
#   logit(p[i]) <- phi[i]
#   phi[i] ~ dnorm(mu, tau)
# }

mu ~ dunif(-10, 10)
tau <- 1/pow(sd, 2)
sd ~ dunif(0, 10)

for (i in 1:N){
phi.pred[i] ~ dnorm(mu, tau)
logit(p.pred[i]) <- phi.pred[i]
y.pred[i] ~ dbin(p.pred[i], n[i])
P.prior[i] <- step(y.pred[i] - y[i] - 0.00001) +
0.5 * equals(y.pred[i], y[i])
}
}"

election.prior.inits <-
  list(list(mu = 0,
            sd = 1,
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("mu", "tau", "P.prior", "p.pred")

election.prior.jag <- jags.model(textConnection(election.prior.model),
                                 data = election.data,
                                 inits = election.prior.inits,
                                 n.chains = 1)
election.prior.coda <- coda.samples(model = election.prior.jag,
                                    variable.names = variable.names,
                                    n.iter = 50000)
summary(election.prior.coda)

#                   Mean        SD  Naive SE Time-series SE
# P.prior[1]   9.017e-01 1.987e-01 8.887e-04      8.887e-04
# P.prior[2]   4.985e-01 5.000e-01 2.236e-03      2.253e-03
# P.prior[3]   5.127e-01 4.998e-01 2.235e-03      2.235e-03
# P.prior[4]   5.283e-01 4.992e-01 2.232e-03      2.232e-03
# P.prior[5]   5.504e-01 4.974e-01 2.224e-03      2.224e-03
# P.prior[6]   5.365e-01 4.987e-01 2.230e-03      2.230e-03
# P.prior[7]   5.775e-01 4.939e-01 2.209e-03      2.209e-03
# P.prior[8]   4.951e-01 5.000e-01 2.236e-03      2.263e-03
# P.prior[9]   5.031e-01 5.000e-01 2.236e-03      2.236e-03
# P.prior[10]  5.444e-01 4.980e-01 2.227e-03      2.227e-03
# P.prior[11]  5.148e-01 4.998e-01 2.235e-03      2.212e-03
# P.prior[12]  5.333e-01 4.989e-01 2.231e-03      2.232e-03

# The prior predictive p-value does clearly identify Northern Ireland as an
# outlier in this case, although not quite as strongly as other approaches.

# Question 3(b): posterior predictive
election.post.model <- "model {
for (i in 1:N){
y[i] ~ dbin(p[i], n[i])
logit(p[i]) <- phi[i]
phi[i] ~ dnorm(mu, tau)

# Draw a prediction for y from the posterior distribution
y.post[i] ~ dbin(p[i], n[i])
P.post[i] <- step(y.post[i] - y[i] - 0.00001) +
0.5 * equals(y.post[i], y[i])
}
mu ~ dunif(-10, 10)
tau <- 1/pow(sd, 2)
sd ~ dunif(0, 10)
}"

election.post.inits <-
  list(list(mu = 0,
            sd = 1,
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("mu", "tau", "p", "P.post")

election.post.jag <- jags.model(textConnection(election.post.model),
                                data = election.data,
                                inits = election.post.inits,
                                n.chains = 1)
election.post.coda <- coda.samples(model = election.post.jag,
                                   variable.names = variable.names,
                                   n.iter = 50000)
summary(election.post.coda)

#                 Mean       SD  Naive SE Time-series SE
# P.post[1]   0.817870 0.240613 1.076e-03      1.571e-03
# P.post[2]   0.498800 0.497357 2.224e-03      2.433e-03
# P.post[3]   0.493970 0.496793 2.222e-03      2.423e-03
# P.post[4]   0.495690 0.497244 2.224e-03      2.417e-03
# P.post[5]   0.498380 0.495149 2.214e-03      2.424e-03
# P.post[6]   0.495780 0.497420 2.225e-03      2.481e-03
# P.post[7]   0.501160 0.496411 2.220e-03      2.411e-03
# P.post[8]   0.501840 0.497927 2.227e-03      2.407e-03
# P.post[9]   0.502560 0.496999 2.223e-03      2.427e-03
# P.post[10]  0.491500 0.495452 2.216e-03      2.431e-03
# P.post[11]  0.495940 0.497311 2.224e-03      2.401e-03
# P.post[12]  0.497650 0.496955 2.222e-03      2.464e-03

# Using the posterior predictive p-value, Northern Ireland has the most extreme
# p-value, but nowhere near as extreme a p-value as under the full leave-one-
# out approach or the mixed p-value. This is an example of the conservativeness
# of the posterior predictive p-value

# Question 4: using DUP results
election.data$y[1] <- round(184260/1000)

# Leaving out Northern Ireland
election.model <- "model {
for (i in 2:N){
y[i] ~ dbin(p[i], n[i])
logit(p[i]) <- phi[i]
phi[i] ~ dnorm(mu, tau)
}
mu ~ dunif(-10, 10)
tau <- 1/pow(sd, 2)
sd ~ dunif(0, 10)

phi.pred[1] ~ dnorm(mu, tau)
logit(p.pred[1]) <- phi.pred[1]
y.pred[1] ~ dbin(p.pred[1], n[1])
P.value[1] <- step(y.pred[1] - y[1] - 0.00001) +
0.5 * equals(y.pred[1], y[1])
}"

election.inits <-
  list(list(mu = 0,
            sd = 1,
            .RNG.seed = 1,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("mu", "tau", "p", "P.value", "p.pred")

election.jag <- jags.model(textConnection(election.model),
                           data = election.data,
                           inits = election.inits,
                           n.chains = 1)
election.coda <- coda.samples(model = election.jag,
                              variable.names = variable.names,
                              n.iter = 50000)

#            Mean       SD  Naive SE Time-series SE
# P.value  0.7702 0.419684 1.877e-03      1.944e-03

# approximate mixed p-value
election.model <- "model {
for (i in 1:N){
y[i] ~ dbin(p[i], n[i])
logit(p[i]) <- phi[i]
phi[i] ~ dnorm(mu, tau)
}
mu ~ dunif(-10, 10)
tau <- 1/pow(sd, 2)
sd ~ dunif(0, 10)

for (i in 1:N){
phi.pred[i] ~ dnorm(mu, tau)
logit(p.pred[i]) <- phi.pred[i]
y.pred[i] ~ dbin(p.pred[i], n[i])
P.mixed[i] <- step(y.pred[i] - y[i] - 0.00001) +
0.5 * equals(y.pred[i], y[i])
}
}"

election.inits <-
  list(list(mu = 0,
            sd = 1,
            .RNG.seed = 87,
            .RNG.name = "base::Mersenne-Twister"))

variable.names <- c("mu", "tau", "p", "P.mixed", "p.pred")

election.jag <- jags.model(textConnection(election.model),
                           data = election.data,
                           inits = election.inits,
                           n.chains = 1)
election.coda <- coda.samples(model = election.jag,
                              variable.names = variable.names,
                              n.iter = 50000)

#                Mean       SD  Naive SE Time-series SE
# P.mixed[1]   0.7578 0.427426 1.912e-03      1.960e-03
# P.mixed[2]   0.1618 0.368110 1.646e-03      1.711e-03
# P.mixed[3]   0.2662 0.441690 1.975e-03      2.010e-03
# P.mixed[4]   0.4914 0.499710 2.235e-03      2.316e-03
# P.mixed[5]   0.7684 0.421121 1.883e-03      1.943e-03
# P.mixed[6]   0.5989 0.489915 2.191e-03      2.226e-03
# P.mixed[7]   0.9575 0.201336 9.004e-04      9.287e-04
# P.mixed[8]   0.1246 0.330139 1.476e-03      1.545e-03
# P.mixed[9]   0.1990 0.399048 1.785e-03      1.828e-03
# P.mixed[10]  0.7121 0.452218 2.022e-03      2.049e-03
# P.mixed[11]  0.2995 0.457847 2.048e-03      2.068e-03
# P.mixed[12]  0.5545 0.496639 2.221e-03      2.239e-03

# p-value is now not extreme for Northern Ireland.
# The approximate check suggests that region 7 (which is Scotland) is now
# much more extreme.

##################################################
### QUESTION 3: JOINT MODELLING

#install.packages("JM")  ### For the AIDS data
library(JM)

# (Start formatting data into a suitable data list)
#install.packages("reshape2")
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

#                          2.5%      25%      50%      75%   97.5%
# beta.mu[1]            4.54706  5.95372  6.67186  7.38000  8.8251
# beta.mu[2]           -0.36130 -0.17549 -0.07358  0.02763  0.2431
# beta.sigma.cov[1,1]  11.28538 16.40801 20.57330 26.31537 42.9795
# beta.sigma.cov[2,1]  -1.87571 -0.81323 -0.44794 -0.13000  0.5094
# beta.sigma.cov[1,2]  -1.87571 -0.81323 -0.44794 -0.13000  0.5094
# beta.sigma.cov[2,2]   0.06178  0.10274  0.13738  0.18927  0.3638
# beta4                -3.98883 -1.70487 -0.79501  0.20210  2.4637
# beta5                -2.23006 -0.53051  0.15274  0.82236  2.4099
# r1                   -1.46464 -0.87887 -0.62395 -0.42551 -0.1364
# r2                   -5.54341 -1.68233 -0.02328  2.04557  8.9820

# The 95% CI for r1 is negative, suggesting that a higher initial value for
# CD4 leads to a lower "surv.mu", indicating longer expected survival.
# The 95% CI for r2 is very wide, suggesting little relationship with the
# individual's random slope.


# Question 2: longitudinal only
aids.model <- "
model {
for (i in 1:I){
for (j in 1:J){
CD4[i, j] ~ dnorm(CD4.mu[i, j], CD4.sigma.prec)
CD4.mu[i, j] <- beta[i, 1] +
beta[i, 2] * t.obs[j] +
beta3 * t.obs[j] * drug.ddI[i]
}

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
CD4.sigma.prec <- 1/pow(CD4.sigma.sd, 2)
CD4.sigma.sd ~ dunif(0, 100)
}
"
variable.names <- c("beta.mu", "beta.sigma.cov", "CD4.pred", "CD4.pred.mu")

aids.inits <- list(list(
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

#                          2.5%      25%      50%      75%   97.5%
# beta.mu[1]            4.58326  6.02217  6.72908  7.45808  8.8259
# beta.mu[2]           -0.34781 -0.14470 -0.04699  0.05613  0.2788
# beta.sigma.cov[1,1]  11.17727 16.72186 20.75894 26.10114 41.8508
# beta.sigma.cov[2,1]  -2.03920 -0.90788 -0.52990 -0.19612  0.4044
# beta.sigma.cov[1,2]  -2.03920 -0.90788 -0.52990 -0.19612  0.4044
# beta.sigma.cov[2,2]   0.06556  0.11263  0.15264  0.20865  0.4212
