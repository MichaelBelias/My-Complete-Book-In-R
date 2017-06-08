library(R2OpenBUGS)
library(coda)

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
teachers.model <- function(){
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
}

ini <- list(list(score.sd = 1, beta1 = 0, beta2 = 0))
parameters.to.save <- c("beta1", "beta2", "score.sd", "mu", "score.pred")

teachers.sim <- bugs(data = teachers.data,
                     inits = ini,
                     model = teachers.model,
                     n.chains = 1,
                     n.burnin = 1000,
                     n.iter = 50000,
                     n.thin = 1,
                     parameters.to.save = parameters.to.save,
                     DIC = FALSE)
teachers.coda <- as.mcmc.list(teachers.sim)

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

# 2.5%      25%     50%     75%  97.5%
# beta1  0.46950  0.57660 0.63290 0.68750 0.7939
# beta2 -0.05902 -0.01864 0.00258 0.02408 0.0647

# Question 2: random intercept model
teachers2.model <- function(){
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
}

ini <- list(list(score.sd = 1, beta1.sd = 1, beta2 = 0))
parameters.to.save <- c("beta1", "beta1.mean", "beta1.sd", "beta2", "score.sd",
                        "mu", "score.pred")

teachers2.sim <- bugs(data = teachers.data,
                     inits = ini,
                     model = teachers2.model,
                     n.chains = 1,
                     n.burnin = 1000,
                     n.iter = 50000,
                     n.thin = 1,
                     parameters.to.save = parameters.to.save,
                     DIC = FALSE)
teachers2.coda <- as.mcmc.list(teachers2.sim)

# Model fit for individual 13
id13 <- which(teachers.data$id == 13)
plot(x = teachers.data$experience[id13],
     y = teachers.data$score[id13],
     ylim = range(teachers.data$score))
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers2.coda,
          "mu",
          indices = id13)

# Predictions for individual 13
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers2.coda,
          "score.pred",
          indices = id13)

# Model fit for individual 16
id16 <- which(teachers.data$id == 16)
plot(x = teachers.data$experience[id16],
     y = teachers.data$score[id16],
     ylim = range(teachers.data$score))
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers2.coda,
          "mu",
          indices = id16)

# Predictions for individual 16
model.fit(teachers.data$score,
          teachers.data$experience,
          teachers2.coda,
          "score.pred",
          indices = id16)

#                       2.5%       25%      50%       75%   97.5%
# beta1.mean       0.4637975  0.559000  0.60870  0.659400 0.75590
# beta1.sd         0.2775000  0.326000  0.35420  0.386100 0.45740
# beta2           -0.0254012  0.002724  0.01686  0.031080 0.05758

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
teachers3.model <- function(){
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
}

ini <- list(list(score.sd = 1, beta.mean = c(0, 0)))
parameters.to.save <- c("beta", "beta.mean", "beta.cov", "score.sd", "mu",
                        "score.pred")

teachers3.sim <- bugs(data = teachers.data,
                     inits = ini,
                     model = teachers3.model,
                     n.chains = 1,
                     n.burnin = 1000,
                     n.iter = 50000,
                     n.thin = 1,
                     parameters.to.save = parameters.to.save,
                     DIC = FALSE)
teachers3.coda <- as.mcmc.list(teachers3.sim)

summary(teachers3.coda)

#                     2.5%       25%       50%       75%    97.5%
# beta.cov[1,1]    0.217500  0.309900  0.372600  0.445300  0.62670
# beta.cov[1,2]   -0.157402 -0.102300 -0.079530 -0.060737 -0.03175
# beta.cov[2,1]   -0.157402 -0.102300 -0.079530 -0.060737 -0.03175
# beta.cov[2,2]    0.031930  0.043150  0.050830  0.060230  0.08365
# beta.mean[1]     0.380098  0.509000  0.575350  0.641600  0.77050
# beta.mean[2]    -0.040630  0.008755  0.034045  0.059150  0.10770

# The covariance between beta[1] and beta[2] has median -0.08
# (95% CI -0.16, -0.03, suggesting a higher rate of improvement for teachers
# whose initial rating is low

# Question 4: Predictions for a new teacher
teachers4.model <- function(){
  for (i in 1:Nmeas){
    score[i] ~ dnorm(mu[i], score.prec)
    mu[i] <- beta[id[i], 1] + beta[id[i], 2] * experience[i]
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
}

ini <- list(list(score.sd = 1, beta.mean = c(0, 0)))
parameters.to.save <- c("beta", "beta.mean", "beta.cov", "score.sd", "mu",
                        "beta.new", "mu.new")

teachers4.sim <- bugs(data = teachers.data,
                     inits = ini,
                     model = teachers4.model,
                     n.chains = 1,
                     n.burnin = 1000,
                     n.iter = 50000,
                     n.thin = 1,
                     parameters.to.save = parameters.to.save,
                     DIC = FALSE)
teachers4.coda <- as.mcmc.list(teachers4.sim)

summary(teachers4.coda)

# beta.new[1]   -0.679102  0.156200  0.5785000  0.987300  1.81400
# beta.new[2]   -0.421602 -0.119300  0.0330300  0.187325  0.49580
# mu.new[1]     -0.679102  0.156200  0.5785000  0.987300  1.81400
# mu.new[2]     -0.437702  0.259500  0.6107000  0.956700  1.64600
# mu.new[3]     -0.386620  0.300375  0.6414000  0.988600  1.66300
# mu.new[4]     -0.517102  0.275675  0.6773000  1.081000  1.87000

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

# Question 1: leaving out Northern Ireland
election.model <- function(){
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

  # mid p-value for unit 1 (Northern Ireland)
  P.value <- step(y.pred[1] - y[1] - 0.00001) +
                  0.5 * equals(y.pred[1], y[1])
}

election.inits <- list(list(mu = 0, sd = 1))
parameters.to.save <- c("mu", "tau", "P.value")

election.sim <- bugs(data = election.data,
                     inits = election.inits,
                     model = election.model,
                     n.chains = 1,
                     n.burnin = 1000,
                     n.iter = 10000,
                     n.thin = 1,
                     parameters.to.save = parameters.to.save,
                     DIC = FALSE)
election.coda <- as.mcmc.list(election.sim)

summary(election.coda[,"P.value"])

#            Mean        SD  Naive SE Time-series SE
# P.value  1.0000 0.0000000 0.000e+00      0.000e+00
#
# p-value is 1, suggesting very strong evidence that the model
# is not adequate for Northern Ireland (as expected!)

# Question 2: mixed approximate p-values
election.model <- function(){
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
}

election.inits <- list(list(mu = 0, sd = 1))
parameters.to.save <- c("mu", "tau", "P.mixed")

election.sim <- bugs(data = election.data,
                     inits = election.inits,
                     model = election.model,
                     n.chains = 1,
                     n.burnin = 1000,
                     n.iter = 10000,
                     n.thin = 1,
                     parameters.to.save = parameters.to.save,
                     DIC = FALSE)
election.coda <- as.mcmc.list(election.sim)
summary(election.coda)

#                Mean      SD  Naive SE Time-series SE
# P.mixed[1]   0.9940 0.05445 0.0005739      0.0005832
# P.mixed[10]  0.4816 0.49955 0.0052657      0.0052657
# P.mixed[11]  0.3523 0.47772 0.0050356      0.0051827
# P.mixed[12]  0.4241 0.49418 0.0052091      0.0052091
# P.mixed[2]   0.2921 0.45476 0.0047936      0.0047178
# P.mixed[3]   0.3394 0.47343 0.0049904      0.0050848
# P.mixed[4]   0.3984 0.48961 0.0051609      0.0051609
# P.mixed[5]   0.4974 0.49966 0.0052669      0.0052669
# P.mixed[6]   0.4362 0.49585 0.0052267      0.0052628
# P.mixed[7]   0.6243 0.48415 0.0051034      0.0053207
# P.mixed[8]   0.2747 0.44637 0.0047052      0.0047052
# P.mixed[9]   0.3060 0.46080 0.0048572      0.0048833

# The p-value for Northern Ireland is about 0.99, so there is slight
# attenuation of the p-value. The Northern Ireland data are a relatively
# large proportion of the data in this example, and it is pulling the
# entire posterior distribution towards it, as any data point would do.

# Question 3(a): prior predictive
election.prior.model <- function(){
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
}

election.prior.inits <- list(list(mu = 0, sd = 1))

parameters.to.save <- c("mu", "tau", "P.prior")
election.prior.sim <- bugs(data = election.data,
                           inits = election.prior.inits,
                           model = election.prior.model,
                           n.chains = 1,
                           n.burnin = 1000,
                           n.iter = 10000,
                           n.thin = 1,
                           parameters.to.save = parameters.to.save,
                           DIC = FALSE)
election.prior.coda <- as.mcmc.list(election.prior.sim)
summary(election.prior.coda)

#                  Mean        SD  Naive SE Time-series SE
# P.prior[1]    0.90406    0.1969  0.002076       0.002076
# P.prior[10]   0.54022    0.4984  0.005253       0.005253
# P.prior[11]   0.50856    0.5000  0.005270       0.005270
# P.prior[12]   0.53433    0.4988  0.005258       0.005258
# P.prior[2]    0.50461    0.5000  0.005270       0.005195
# P.prior[3]    0.51194    0.4999  0.005269       0.005269
# P.prior[4]    0.52789    0.4992  0.005263       0.005363
# P.prior[5]    0.54506    0.4980  0.005249       0.005091
# P.prior[6]    0.53256    0.4990  0.005260       0.005260
# P.prior[7]    0.57444    0.4944  0.005211       0.005211
# P.prior[8]    0.49700    0.5000  0.005271       0.005271
# P.prior[9]    0.50644    0.5000  0.005270       0.005270

# The prior predictive p-value does clearly identify Northern Ireland as an
# outlier in this case, although not quite as strongly as other approaches.

# Question 3(b): posterior predictive
election.post.model <- function(){
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
}

election.post.inits <- list(list(mu = 0, sd = 1))

parameters.to.save <- c("mu", "tau", "P.post")

election.post.sim <- bugs(data = election.data,
                          inits = election.post.inits,
                          model = election.post.model,
                          n.chains = 1,
                          n.burnin = 1000,
                          n.iter = 10000,
                          n.thin = 1,
                          parameters.to.save = parameters.to.save,
                          DIC = FALSE)
election.post.coda <- as.mcmc.list(election.post.sim)
summary(election.post.coda)

#               Mean     SD  Naive SE Time-series SE
# P.post[1]   0.8212 0.2397  0.002526       0.002955
# P.post[10]  0.4951 0.4955  0.005223       0.005444
# P.post[11]  0.4983 0.4969  0.005238       0.005238
# P.post[12]  0.4899 0.4974  0.005243       0.005243
# P.post[2]   0.4883 0.4972  0.005241       0.005154
# P.post[3]   0.4969 0.4969  0.005238       0.005335
# P.post[4]   0.4883 0.4975  0.005244       0.005244
# P.post[5]   0.4889 0.4952  0.005219       0.005219
# P.post[6]   0.4963 0.4969  0.005237       0.005237
# P.post[7]   0.4972 0.4961  0.005230       0.005230
# P.post[8]   0.4860 0.4979  0.005249       0.005249
# P.post[9]   0.5036 0.4976  0.005245       0.005392

# Using the posterior predictive p-value, Northern Ireland has the most extreme
# p-value, but nowhere near as extreme a p-value as under the full leave-one-
# out approach or the mixed p-value. This is an example of the conservativeness
# of the posterior predictive p-value

# Question 4: using DUP results
election.data$y[1] <- round(184260/1000)

# Leaving out Northern Ireland
election.model <- function(){
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
}

election.inits <- list(list(mu = 0, sd = 1))
parameters.to.save <- c("mu", "tau", "P.value")

election.sim <- bugs(data = election.data,
                     inits = election.inits,
                     model = election.model,
                     n.chains = 1,
                     n.burnin = 1000,
                     n.iter = 10000,
                     n.thin = 1,
                     parameters.to.save = parameters.to.save,
                     DIC = FALSE)
election.coda <- as.mcmc.list(election.sim)

# Mean     SD Naive SE Time-series SE
# mu         -0.6016 0.1947 0.002053        0.00209
# P.value[1]  0.7725 0.4183 0.004410        0.00441

# approximate mixed p-value
election.model <- function(){
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
}

election.inits <- list(list(mu = 0, sd = 1))
parameters.to.save <- c("mu", "tau", "P.mixed")

election.sim <- bugs(data = election.data,
                     inits = election.inits,
                     model = election.model,
                     n.chains = 1,
                     n.burnin = 1000,
                     n.iter = 10000,
                     n.thin = 1,
                     parameters.to.save = parameters.to.save,
                     DIC = FALSE)
election.coda <- as.mcmc.list(election.sim)

#                Mean     SD Naive SE Time-series SE
# P.mixed[1]   0.7590 0.4267 0.004498       0.004364
# P.mixed[10]  0.7233 0.4472 0.004714       0.004714
# P.mixed[11]  0.3061 0.4606 0.004855       0.004855
# P.mixed[12]  0.5594 0.4961 0.005230       0.005230
# P.mixed[2]   0.1599 0.3665 0.003863       0.003863
# P.mixed[3]   0.2687 0.4430 0.004670       0.004572
# P.mixed[4]   0.4871 0.4996 0.005267       0.005267
# P.mixed[5]   0.7717 0.4191 0.004417       0.004417
# P.mixed[6]   0.6063 0.4884 0.005148       0.005148
# P.mixed[7]   0.9565 0.2037 0.002147       0.002184
# P.mixed[8]   0.1283 0.3344 0.003525       0.003586
# P.mixed[9]   0.2043 0.4029 0.004247       0.004314

# p-value is now not extreme for Northern Ireland.
# The approximate check suggests that region 7 (which is Scotland) is now
# much more extreme.


##################################################
### QUESTION 3: JOINT MODELLING

install.packages("JM")  ### For the AIDS data
library(JM)

# (Start formatting data into a suitable data list)
install.packages("reshape2")
library(reshape2)

# Just use the first 20 patients
aids <- droplevels(subset(aids, patient %in% as.character(1:20)))
aids.id <- droplevels(subset(aids.id, patient %in% as.character(1:20)))

# clearer/shorter names
colnames(aids)[2] <- "surv.time"
colnames(aids)[5] <-  "t.obs"

aids.CD4.df <- dcast(patient ~ t.obs, value.var = "CD4", data = aids)[, -1]
aids.drug.df <- dcast(patient ~ t.obs, value.var = "drug", data = aids)[, 1:2]
aids.drug.ddI <- ifelse(aids.drug.df[, 2] == "ddI", 1, 0)

aids.npatients <- nrow(aids.CD4.df)

aids.t.obs <- sort(unique(aids$t.obs))
aids.n.t.obs <- length(aids.t.obs)

aids.t.pred <- seq(from = 0, to = max(aids.t.obs), by = 0.5)
aids.n.t.preds <- length(aids.t.pred)
aids.npredpatients <- 10

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
                  CD4 = as.matrix(aids.CD4.df),
                  drug.ddI = aids.drug.ddI,
                  t.obs = aids.t.obs,
                  J = aids.n.t.obs,
                  surv.c = ifelse(aids.id$death, 0, aids.id$Time),
                  surv.time = ifelse(aids.id$death, aids.id$Time, NA),
                  break.time = aids.break.time,
                  # Predictions
                  K = aids.npredpatients,
                  drug.ddI.pred = aids.pred.drug.ddI,
                  t.pred = aids.t.pred,
                  L = aids.n.t.preds)
# (End formatting data into a suitable data list)

aids.inits <- list(list(
    CD4 = matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, 6.809, 8.011, NA, NA, NA, 5.675, NA, NA, NA,
                   NA, 9.448, NA, NA, NA, NA, NA, NA, NA, 3.495, NA, NA, NA, NA,
                   NA, NA, 10.19, NA, NA, 6.381, NA, 9.398, NA, NA, 2.547, NA, NA,
                   NA, NA, 6.655, NA, NA, NA, 8.812, NA, NA, 14.99, NA, -2.036,
                   12.06, NA, 13.06, NA, 2.432, 1.308, NA, 5.809, NA, NA, 7.694,
                   NA, 12.73, NA, 6.229, 3.44, 5.212, 15.13, 4.421, -6.638, 18.13,
                   8.331, 13.58, NA, 3.661, 2.666, 6.087, 9.511, 5.65, 4.515,
                   15.04, 15.81), nrow = 20, ncol = 5),
    CD4.pred = matrix(c(7.791, 13.16, 10.91, 13.56, 9.057, 13.83, 7.352, 14.11,
                        12.46, 11, 7.054, 14.72, 8.335, 11.33, 7.58, 12.82, 11.15,
                        14.48, 9.761, 11.07, 7.793, 13.77, 9.651, 13.3, 7.171,
                        13.13, 8.161, 14.4, 8.072, 13.27, 8.685, 12.5, 9.816,
                        14.54, 10.4, 12.24, 10.75, 15.17, 11.83, 13.18, 9.485,
                        11.44, 8.348, 13.07, 9.595, 12.85, 10.35, 9.926, 9.565,
                        13.17, 11.48, 11.74, 10.28, 11.32, 7.707, 10.35, 12.28,
                        10.96, 10.9, 11.49, 10.62, 10.72, 10.47, 15.29, 9.581,
                        14.77, 10.9, 11.9, 9.788, 10.87, 9.039, 15.07, 10.63,
                        12.46), nrow = 2, ncol = 37),
    CD4.sigma.sd = 1.469,
    beta = matrix(c(7.59, 6.547, 5.031, 4.346, 9.251, 5.613, 6.525, 2.635, 4.053,
                    10.55, 8.262, 13.92, 3.59, 2.925, 10.02, 3.226, 4.638, 4.588,
                    4.503, 19.19, 0.2846, 0.09968, 0.2863, -0.04011, -0.01887,
                    0.6014, -0.1435, -0.3375, 0.6405, 0.06098, 0.3629, -0.06609,
                    0.09016, 0.07878, 0.04751, 0.2625, 0.05186, 0.02015, 0.6716,
                    -0.1723), nrow = 20, ncol = 2),
    beta.mu = c(6.55, 0.2097),
    beta.pred = matrix(c(8.684, 13.55, 0.1315, 0.06223), nrow = 2, ncol = 2),
    beta.sigma.prec = matrix(c(0.08727, -0.008303, -0.008303, 10.23),
                             nrow = 2, ncol = 2),
    beta3 = -0.1595,
    beta4 = -1.163,
    beta5 = 0.9403,
    r1 = -0.7182,
    r2 = 3.727,
    surv.time = c(614.9, 250.2, NA, 92.54, 654.2, NA, 853.7, NA, NA, 634.9, NA,
                  100, NA, NA, 98.77, NA, 20.12, 69.08, NA, 100)
))
parameters.to.save <- c("beta.mu", "beta.sigma.cov", "CD4.pred", "CD4.pred.mu",
                        "r1", "r2", "beta4", "beta5")

# takes about 4 mins
aids.sim <- bugs(data = aids.data,
                 inits = aids.inits,
                 model = "hier2-joint.aids.txt",
                 n.chains = 1,
                 n.burnin = 10000,
                 n.iter = 25000,
                 n.thin = 20,
                 parameters.to.save = parameters.to.save,
                 DIC = FALSE)

aids.coda <- as.mcmc.list(aids.sim)

#                          2.5%      25%       50%      75%   97.5%
# beta.mu[1]            4.48615  5.96875  6.682000  7.40350  8.7189
# beta.mu[2]           -0.35974 -0.15520 -0.050630  0.03639  0.2068
# beta.sigma.cov[1,1]  11.48450 17.08250 21.130000 26.54500 42.2110
# beta.sigma.cov[1,2]  -1.82675 -0.85760 -0.459750 -0.12433  0.5074
# beta.sigma.cov[2,1]  -1.82675 -0.85760 -0.459750 -0.12433  0.5074
# beta.sigma.cov[2,2]   0.06368  0.10415  0.141250  0.19368  0.3930
# beta4                -3.54343 -1.52625 -0.549450  0.33340  2.0448
# beta5                -1.97447 -0.59390  0.118150  0.83110  2.5610
# r1                   -1.28020 -0.90638 -0.703200 -0.50495 -0.1200
# r2                   -5.73175 -1.45175  0.002972  2.39900  9.1626

# The 95% CI for r1 is negative, suggesting that a higher initial value for
# CD4 leads to a lower "surv.mu", indicating longer expected survival.
# The 95% CI for r2 is very wide, suggesting little relationship with the
# individual's random slope.

# Question 2
aids.model <- function(){
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

aids.inits <- list(list(
  CD4 = matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                 NA, NA, NA, NA, 6.809, 8.011, NA, NA, NA, 5.675, NA, NA, NA,
                 NA, 9.448, NA, NA, NA, NA, NA, NA, NA, 3.495, NA, NA, NA, NA,
                 NA, NA, 10.19, NA, NA, 6.381, NA, 9.398, NA, NA, 2.547, NA, NA,
                 NA, NA, 6.655, NA, NA, NA, 8.812, NA, NA, 14.99, NA, -2.036,
                 12.06, NA, 13.06, NA, 2.432, 1.308, NA, 5.809, NA, NA, 7.694,
                 NA, 12.73, NA, 6.229, 3.44, 5.212, 15.13, 4.421, -6.638, 18.13,
                 8.331, 13.58, NA, 3.661, 2.666, 6.087, 9.511, 5.65, 4.515,
                 15.04, 15.81), nrow = 20, ncol = 5),
  CD4.pred = matrix(c(7.791, 13.16, 10.91, 13.56, 9.057, 13.83, 7.352, 14.11,
                      12.46, 11, 7.054, 14.72, 8.335, 11.33, 7.58, 12.82, 11.15,
                      14.48, 9.761, 11.07, 7.793, 13.77, 9.651, 13.3, 7.171,
                      13.13, 8.161, 14.4, 8.072, 13.27, 8.685, 12.5, 9.816,
                      14.54, 10.4, 12.24, 10.75, 15.17, 11.83, 13.18, 9.485,
                      11.44, 8.348, 13.07, 9.595, 12.85, 10.35, 9.926, 9.565,
                      13.17, 11.48, 11.74, 10.28, 11.32, 7.707, 10.35, 12.28,
                      10.96, 10.9, 11.49, 10.62, 10.72, 10.47, 15.29, 9.581,
                      14.77, 10.9, 11.9, 9.788, 10.87, 9.039, 15.07, 10.63,
                      12.46), nrow = 2, ncol = 37),
  CD4.sigma.sd = 1.469,
  beta = matrix(c(7.59, 6.547, 5.031, 4.346, 9.251, 5.613, 6.525, 2.635, 4.053,
                  10.55, 8.262, 13.92, 3.59, 2.925, 10.02, 3.226, 4.638, 4.588,
                  4.503, 19.19, 0.2846, 0.09968, 0.2863, -0.04011, -0.01887,
                  0.6014, -0.1435, -0.3375, 0.6405, 0.06098, 0.3629, -0.06609,
                  0.09016, 0.07878, 0.04751, 0.2625, 0.05186, 0.02015, 0.6716,
                  -0.1723), nrow = 20, ncol = 2),
  beta.mu = c(6.55, 0.2097),
  beta.pred = matrix(c(8.684, 13.55, 0.1315, 0.06223), nrow = 2, ncol = 2),
  beta.sigma.prec = matrix(c(0.08727, -0.008303, -0.008303, 10.23),
                           nrow = 2, ncol = 2),
  beta3 = -0.1595
))
parameters.to.save <- c("beta.mu", "beta.sigma.cov", "CD4.pred", "CD4.pred.mu")
aids.sim <- bugs(data = aids.data,
                 inits = aids.inits,
                 model = aids.model,
                 n.chains = 1,
                 n.burnin = 10000,
                 n.iter = 25000,
                 n.thin = 20,
                 parameters.to.save = parameters.to.save,
                 DIC = FALSE)
aids.coda <- as.mcmc.list(aids.sim)

# beta.mu[1]            4.64980  6.01075  6.81850  7.54550  8.7580
# beta.mu[2]           -0.34600 -0.16660 -0.06913  0.03044  0.2322
# beta.sigma.cov[1,1]  11.48175 16.34750 20.93000 25.58750 39.1557
# beta.sigma.cov[1,2]  -2.04195 -0.89485 -0.47995 -0.20312  0.3500
# beta.sigma.cov[2,1]  -2.04195 -0.89485 -0.47995 -0.20312  0.3500
# beta.sigma.cov[2,2]   0.06616  0.11315  0.15205  0.20850  0.4407
