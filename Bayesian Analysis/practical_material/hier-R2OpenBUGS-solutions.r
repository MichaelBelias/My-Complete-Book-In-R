library(R2OpenBUGS)
library(coda)

##################################################
### QUESTION 1:  META-ANALYSIS

# control data is in the first column, treatment in the second
sepsis.data <- list(
  Ns = 10,
  y = matrix(c(23, 8, 5, 14, 209, 5, 13, 13, 8, 39, 20, 2, 0, 8, 186, 4,
               10, 19, 3, 40),
               nrow = 10,
               ncol = 2),
  n = matrix(c(65, 43, 59, 32, 1212, 50, 34, 41, 40, 381, 61, 43, 56, 34,
               1204, 100, 68, 40, 40, 372),
               nrow = 10,
               ncol = 2)
)

sepsis.re.model <- function(){
  for (s in 1:Ns){ # for each study, Ns = total number of studies
    for (a in 1:2){ # for each arm
      y[s, a] ~ dbin(p[s, a], n[s, a])
      logit(p[s, a]) <- alpha[s] + mu[s, a]
    }
    # treatment effect
    mu[s, 1] <- 0
    mu[s, 2] ~ dnorm(delta, prec.mu)

    ## study-specific baselines, vague priors
    alpha[s] ~ dnorm(0, 0.01)
  }

  delta ~ dnorm(0, 0.01)

  prec.mu <- 1/(sd.mu * sd.mu)
  sd.mu ~ dunif(0, 10)
}

ini <- list(list(delta = 0, alpha = rep(0, 10)))

sepsis.re.sim <- bugs(data = sepsis.data,
                      inits = ini,
                      model = sepsis.re.model,
                      n.chains = 1,
                      n.burnin = 1000,
                      n.iter = 50000,
                      n.thin = 1,
                      parameters.to.save = c("p", "alpha", "mu", "delta"),
                      DIC=FALSE)
sepsis.re.coda <- as.mcmc.list(sepsis.re.sim)

summary(sepsis.re.coda)
summary(sepsis.re.coda[, "delta"])

# The 95% posterior CI for the log odds ratio, delta, is (-1.35, -0.01),
# suggesting that IVIG is reducing the chances of in-hospital infection.

# Question 1(b) Monitoring the odds ratio

sepsis.re.model2 <- function(){
  for (s in 1:Ns){ # for each study, Ns = total number of studies
    for (a in 1:2){ # for each arm
      y[s, a] ~ dbin(p[s, a], n[s, a])
      logit(p[s, a]) <- alpha[s] + mu[s, a]
    }
    # treatment effect
    mu[s, 1] <- 0
    mu[s, 2] ~ dnorm(delta, prec.mu)

    ## study-specific baselines, vague priors
    alpha[s] ~ dnorm(0, 0.01)
  }

  delta ~ dnorm(0, 0.01)
  # Add the following line to get odds ratios rather than log odds ratios
  delta.or <- exp(delta)

  prec.mu <- 1/(sd.mu * sd.mu)
  sd.mu ~ dunif(0, 10)
}
parameters.to.save <- c("p", "alpha", "mu", "delta", "delta.or")

sepsis.re.sim <- bugs(data = sepsis.data,
                      inits = ini,
                      model = sepsis.re.model2,
                      n.chains = 1,
                      n.burnin = 1000,
                      n.iter = 50000,
                      n.thin = 1,
                      parameters.to.save = parameters.to.save,
                      DIC=FALSE)
sepsis.re.coda <- as.mcmc.list(sepsis.re.sim)

summary(sepsis.re.coda[, "delta.or"])

# The 95% posterior CI for the odd ratio, delta.or, is (0.26, 0.99),
# with median 0.58
# Note that rather than re-running the whole model, we could simply
# exponentiate the stored samples for delta. This approach is very useful for
# complex models that take a long time to run!

# Question 2 - Fixed effects model
sepsis.fe.model <- function(){
  for (s in 1:Ns){ # for each study, Ns = total number of studies
    for (a in 1:2){ # for each arm
      y[s, a] ~ dbin(p[s, a], n[s, a])
      # Now only have a single delta
      logit(p[s, a]) <- alpha[s] + delta[a]
    }

    ## study-specific baselines, vague priors
    alpha[s] ~ dnorm(0, 0.01)
  }

  delta[1] <- 0
  delta[2] ~ dnorm(0, 0.01)
  delta.or <- exp(delta[2])

  prec.mu <- 1/(sd.mu * sd.mu)
  sd.mu ~ dunif(0, 10)
}

parameters.to.save <- c("p", "alpha", "delta", "delta.or")
ini <- list(list(delta = c(NA, 0), alpha = rep(0, 10)))

sepsis.fe.sim <- bugs(data = sepsis.data,
                      inits = ini,
                      model = sepsis.fe.model,
                      n.chains = 1,
                      n.burnin = 1000,
                      n.iter = 50000,
                      n.thin = 1,
                      parameters.to.save = parameters.to.save,
                      DIC=FALSE)
sepsis.fe.coda <- as.mcmc.list(sepsis.fe.sim)

summary(sepsis.fe.coda[, "delta.or"])

# The 95% posterior CI for the odds ratio is now (0.68, 0.96), median 0.81
# See also the results shown in the lecture
# Under the fixed effect model, the lower bound of the 95% CI for the
# overall odds ratio is considerably higher.
# This is due to the lack of shrinkage in the fixed effects model
# e.g. the Sandberg (2000) study is not shrunk down towards the mean

# Question 3 - US studies
is.us <- c(1, 0, 1, 0, 1, 0, 0, 0, 0, 1)
sepsis.data$is.us <- is.us

sepsis.re.model <- function(){
  for (s in 1:Ns){ # for each study, Ns = total number of studies
    for (a in 1:2){ # for each arm
      y[s, a] ~ dbin(p[s, a], n[s, a])
      logit(p[s, a]) <- alpha[s] + mu[s, a]
    }
    # treatment effect
    mu[s, 1] <- 0
    mu[s, 2] ~ dnorm(mean.mu[s], prec.mu)
    mean.mu[s] <- delta + beta * is.us[s]

    ## study-specific baselines, vague priors
    alpha[s] ~ dnorm(0, 0.01)
  }

  delta ~ dnorm(0, 0.01)
  beta ~ dnorm(0, 0.01)

  prec.mu <- 1/(sd.mu * sd.mu)
  sd.mu ~ dunif(0, 10)
}

ini <- list(
  list(
    delta = 0,
    alpha = rep(0, 10)
  )# ,
#   list(
#     delta = 0.1,
#     alpha = c(1,-1,-2,0,0,-2,1,0,2,2)
#   )
)
parameters.to.save <- c("p", "alpha", "mu", "delta", "beta")

sepsis.re.sim <- bugs(data = sepsis.data,
                      inits = ini,
                      model = sepsis.re.model,
                      n.chains = 1,
                      n.burnin = 1000,
                      n.iter = 50000,
                      n.thin = 10,
                      parameters.to.save = parameters.to.save,
                      DIC=FALSE)
sepsis.re.coda <- as.mcmc.list(sepsis.re.sim)

summary(sepsis.re.coda)
summary(sepsis.re.coda[, "beta"])

# There is a suggestion that the treatment is less effective in the US
# (the median for beta, the log odds ratio, is 0.52), but the 95% CI
# is wide (-0.95, 1.62) and includes 0, so the evidence is really quite weak
# for this effect.
# Fitting a separate model is equivalent, except that prec.mu is shared between
# both the US and no-US studies in the above model, whereas they would be
# unconstrained if US and non-US studies were considered separately. (Of
# course, you could specify a meta-regression model that allowed prec.mu
# to differ between US and non-US studies)

##################################################
### QUESTION 2:  PRIOR SENSITIVITY IN HIERARCHICAL MODELS

Ns <- 6
y <- matrix(c(4, 5, 3, 1, 14, 2, 1, 0, 9, 7, 13, 18), byrow = T, ncol = 2)
n <- matrix(c(11, 13, 14, 13, 19, 17, 12, 13, 13, 17, 47, 48), byrow = T, ncol = 2)

gran.data <- list(Ns = Ns, y = y, n = n)

# Model with vague Normal(0, 10^2) prior
gran.re.vague.model <- function(){
  for (s in 1:Ns){ # for each study, Ns = total number of studies
    for (a in 1:2){ # for each arm
      y[s, a] ~ dbin(p[s, a], n[s, a])
      logit(p[s, a]) <- alpha[s] + mu[s, a]
    }
    # treatment effect
    mu[s, 1] <- 0
    mu[s, 2] ~ dnorm(delta, prec.mu)

    alpha[s] ~ dnorm(0, 0.01)
  }

  delta ~ dnorm(0, 0.01)
  delta.or <- exp(delta)

  prec.mu <- 1/(sd.mu * sd.mu)
  sd.mu ~ dunif(0, 10)
}

ini <- list(
  list(
    delta = 0,
    alpha = rep(0, 6)
  )
)

parameters.to.save <- c("p", "alpha", "mu", "delta", "delta.or")

gran.re.vague.sim <- bugs(data = gran.data,
                          inits = ini,
                          model = gran.re.vague.model,
                          n.chains = 1,
                          n.burnin = 1000,
                          n.iter = 50000,
                          n.thin = 1,
                          parameters.to.save = parameters.to.save,
                          DIC=FALSE)
gran.re.vague.coda <- as.mcmc.list(gran.re.vague.sim)

summary(gran.re.vague.coda)
summary(gran.re.vague.coda[, "delta"])

# With an informative prior
gran.re.scept.model <- function(){
  for (s in 1:Ns){ # for each study, Ns = total number of studies
    for (a in 1:2){ # for each arm
      y[s, a] ~ dbin(p[s, a], n[s, a])
      logit(p[s, a]) <- alpha[s] + mu[s, a]
    }
    # treatment effect
    mu[s, 1] <- 0
    mu[s, 2] ~ dnorm(delta, prec.mu)

    alpha[s] ~ dnorm(0, 0.01)
  }

  # For a Normal prior on the log odds ratio scale want
  # 95% of the region to be between log(0.75)=-0.29 and 0.29, and
  # mean = 0. So SD = 0.29/1.96 = 0.15. And precision = 1/(0.15^2)
  delta ~ dnorm(0, 44)
  delta.or <- exp(delta)

  prec.mu <- 1/(sd.mu * sd.mu)
  sd.mu ~ dunif(0, 10)
}

ini <- list(
  list(
    delta = 0,
    alpha = rep(0, 6)
  )
)

parameters.to.save <- c("p", "alpha", "mu", "delta", "delta.or")
gran.re.scept.sim <- bugs(data = gran.data,
                          inits = ini,
                          model = gran.re.scept.model,
                          n.chains = 1,
                          n.burnin = 1000,
                          n.iter = 50000,
                          n.thin = 1,
                          parameters.to.save = parameters.to.save,
                          DIC=FALSE)
gran.re.scept.coda <- as.mcmc.list(gran.re.scept.sim)

summary(gran.re.scept.coda)
summary(gran.re.scept.coda[, "delta"])

# Model with vague prior
gran.re.vague2.model <- function(){
  for (s in 1:Ns){ # for each study, Ns = total number of studies
    for (a in 1:2){ # for each arm
      y[s, a] ~ dbin(p[s, a], n[s, a])
      logit(p[s, a]) <- alpha[s] + mu[s, a]
    }
    # treatment effect
    mu[s, 1] <- 0
    mu[s, 2] ~ dnorm(delta, prec.mu)

    alpha[s] ~ dnorm(0, 0.01)
  }

  delta ~ dunif(-100, 100)
  delta.or <- exp(delta)

  prec.mu <- 1/(sd.mu * sd.mu)
  sd.mu ~ dunif(0, 10)
}

ini <- list(
  list(
    delta = 0,
    alpha = rep(0, 6)
  )
)

parameters.to.save <- c("p", "alpha", "mu", "delta", "delta.or")

gran.re.vague2.sim <- bugs(data = gran.data,
                          inits = ini,
                          model = gran.re.vague2.model,
                          n.chains = 1,
                          n.burnin = 1000,
                          n.iter = 50000,
                          n.thin = 1,
                          parameters.to.save = parameters.to.save,
                          DIC=FALSE)
gran.re.vague2.coda <- as.mcmc.list(gran.re.vague2.sim)

summary(gran.re.vague2.coda)
summary(gran.re.vague2.coda[, "delta"])

# Under the vague N(0, 10^2) prior delta's posterior media is
# -1.12 (95% CI -3.67, 0.77)
# Under the informative prior, posterior median -0.04 (95% CI -0.33, 0.25)
# Under the U(-100, 100) prior, delta's posterior median is -1.12
# (95% CI -3.80, 0.79)
# The posterior distribution is much more concentrated around 0 under the
# informative prior, but there is no substantive difference between the
# two vague priors


# Question 3
gran.re.model <- function(){
  for (i in 1:Np){
    for (s in 1:Ns){ # for each study, Ns = total number of studies
      for (a in 1:2){ # for each arm
        y[i, s, a] ~ dbin(p[i, s, a], n[i, s, a])
        logit(p[i, s, a]) <- alpha[i, s] + mu[i, s, a]
      }
      # treatment effect
      mu[i, s, 1] <- 0
      mu[i, s, 2] ~ dnorm(delta[i], prec.mu[i])

      alpha[i, s] ~ dnorm(0, 0.01)
    }

    delta[i] ~ dnorm(0, 0.01)
    delta.or[i] <- exp(delta[i])
  }

  sd.mu[1] ~ dunif(0, 10)
  prec.mu[1] <- 1/(sd.mu[1] * sd.mu[1])

  sd.mu[2] ~ dunif(0, 100)
  prec.mu[2] <- 1/(sd.mu[2] * sd.mu[2])

  prec.mu[3] ~ dgamma(1, 1)
  sd.mu[3] <- 1/sqrt(prec.mu[3])

  prec.mu[4] ~ dgamma(0.001, 0.001)
  sd.mu[4] <- 1/sqrt(prec.mu[4])
}

ini <- list(
  list(
    delta = rep(0, 4),
    alpha = matrix(rep(rep(0, 6), 4), nrow = 4),
    prec.mu = c(NA, NA, 1, 1)
  )
)

gran.data.edited <- gran.data
gran.data.edited$Np <- 4
gran.data.edited$y <- array(NA, dim = c(4, dim(gran.data$y)))
gran.data.edited$y[1,,] <- gran.data$y
gran.data.edited$y[2,,] <- gran.data$y
gran.data.edited$y[3,,] <- gran.data$y
gran.data.edited$y[4,,] <- gran.data$y
gran.data.edited$n <- array(NA, dim = c(4, dim(gran.data$n)))
gran.data.edited$n[1,,] <- gran.data$n
gran.data.edited$n[2,,] <- gran.data$n
gran.data.edited$n[3,,] <- gran.data$n
gran.data.edited$n[4,,] <- gran.data$n

parmeters.to.save <- c("p", "alpha", "mu", "delta", "delta.or", "sd.mu",
                       "prec.mu")

gran.re.sim <- bugs(data = gran.data.edited,
                    inits = ini,
                    model = gran.re.model,
                    n.chains = 1,
                    n.burnin = 1000,
                    n.iter = 50000,
                    n.thin = 1,
                    parameters.to.save = parameters.to.save,
                    DIC=FALSE)
gran.re.coda <- as.mcmc.list(gran.re.sim)

summary(gran.re.coda)
summary(gran.re.coda[, "delta[1]"])
summary(gran.re.coda[, "delta[2]"])
summary(gran.re.coda[, "delta[3]"])
summary(gran.re.coda[, "delta[4]"])

# Posterior quantiles under each prior
#
#                   2.5%       25%       50%      75%   97.5%
# delta[1]    -3.640e+00 -1.736000 -1.115000 -0.581475  0.75972
# delta[2]    -3.693e+00 -1.726000 -1.115000 -0.571275  0.77831
# delta[3]    -2.591e+00 -1.485000 -1.035000 -0.605100  0.28281
# delta[4]    -2.991e+00 -1.515000 -1.009000 -0.570175  0.41630
# delta.or[1]  2.626e-02  0.176200  0.328000  0.559025  2.13805
# delta.or[2]  2.491e-02  0.177900  0.328000  0.564800  2.17800
# delta.or[3]  7.494e-02  0.226400  0.355250  0.546000  1.32700
# delta.or[4]  5.023e-02  0.219900  0.364600  0.565425  1.51600
# p[1,1,1]     1.529e-01  0.300300  0.394300  0.492400  0.66880
# p[2,1,1]     1.530e-01  0.300300  0.393100  0.490900  0.67030
# p[3,1,1]     1.713e-01  0.317100  0.408800  0.503900  0.67480
# p[4,1,1]     1.653e-01  0.314975  0.405500  0.500700  0.67300
# p[1,4,2]     1.503e-05  0.001632  0.007224  0.022182  0.09988
# p[2,4,2]     9.776e-06  0.001583  0.007044  0.022190  0.09889
# p[3,4,2]     1.567e-04  0.003243  0.010240  0.026390  0.10140
# p[4,4,2]     6.866e-05  0.002825  0.009885  0.026430  0.10370

# The results under priors 1 and 2 are quite similar, but are fairly different
# to priors 3 and 4, which are fairly similar to each other, with 4 being more
# extreme.
# Note that the probabilities for each under are quite stable across all the
# priors -- the data provide a lot of information about these quantites, and
# so the prior plays a less strong part in these estimates compared to the
# random effect variance parameter.


####

gran.re.model <- function(){
  for (s in 1:Ns){ # for each study, Ns = total number of studies
    for (a in 1:2){ # for each arm
    y[s, a] ~ dbin(p[s, a], n[s, a])
      logit(p[s, a]) <- alpha[s] + mu[s, a]
    }
    # treatment effect
    mu[s, 1] <- 0
    mu[s, 2] ~ dnorm(delta, prec.mu)

    alpha[s] ~ dnorm(0, 0.01)
  }

  delta ~ dnorm(0, 0.01)
  delta.or <- exp(delta)

  prec.mu <- 1/(sd.mu * sd.mu)
  sd.mu <- sqrt(var.d)

  # Turner et al (2012) informative prior
  # This is log normal with mean -3.93 and variance 1.51 on the log scale
  var.d ~ dlnorm(-3.93, 0.44) # 1/(1.51^2)
}

ini <- list(
  list(
    delta = 0,
    alpha = rep(0, 6)
  )
)

parameters.to.save <- c("p", "alpha", "mu", "delta", "delta.or", "sd.mu",
                        "prec.mu")
gran.re.sim <- bugs(data = gran.data,
                    inits = ini,
                    model = gran.re.model,
                    n.chains = 1,
                    n.burnin = 1000,
                    n.iter = 100000,
                    n.thin = 1,
                    parameters.to.save = parameters.to.save,
                    DIC=FALSE)
gran.re.coda <- as.mcmc.list(gran.re.sim)

summary(gran.re.coda)
summary(gran.re.coda[, "delta"])

#                2.5%       25%      50%      75%    97.5%
# delta    -1.7780000 -1.043000 -0.75175 -0.49050  -0.001991
# delta.or  0.1690000  0.352300  0.47150  0.61230   0.998000
# prec.mu   0.5504000  2.057000  5.42200 22.76000 379.307500
# sd.mu     0.0513495  0.209600  0.42950  0.69722   1.348000

# Under the informative prior, the 95% CI for delta is much narrower, since the informative prior allows much lower levels of heterogeneity than our previous vague priors, which leads to greater precision in the treatment effect. The interval also doesn't include 0, in contrast to the other priors!
