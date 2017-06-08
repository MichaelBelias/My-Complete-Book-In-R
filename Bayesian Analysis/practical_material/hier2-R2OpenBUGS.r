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

# Question 3: random slope and intercept model

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
aids.id <- droplevels(subset(aids.id, patient %in% as.character(1:20)))

# clearer/shorter names
colnames(aids)[2] <- "surv.time"
colnames(aids)[5] <-  "t.obs"

aids.CD4.df <- dcast(patient ~ t.obs, value.var = "CD4", data = aids)[, -1]
aids.drug.df <- dcast(patient ~ t.obs, value.var = "drug", data = aids)[, 1:2]
aids.drug.ddI <- ifelse(aids.drug.df[, 2] == "ddI", 1, 0)

aids.npatients <- nrow(aids.CD4.df)
aids.n.t.obs <- length(aids.t.obs)

aids.t.obs <- sort(unique(aids$t.obs))
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

# Question 2
