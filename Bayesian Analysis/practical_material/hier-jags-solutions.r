library(rjags)

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

sepsis.re.model <- "
model
{
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
"

ini <- list(list(delta = 0, alpha = rep(0, 10)))

sepsis.re.jag <- jags.model(textConnection(sepsis.re.model),
                            data = sepsis.data,
                            inits = ini,
                            n.chains = 1)
update(sepsis.re.jag, n.iter = 1000)
sepsis.re.coda <- coda.samples(model = sepsis.re.jag,
                               variable.names = c("p", "alpha", "mu", "delta"),
                               n.iter = 50000)

# summary(sepsis.re.coda)
summary(sepsis.re.coda[, "delta"])

# The 95% posterior CI for the log odds ratio, delta, is (-1.35, -0.01),
# suggesting that IVIG is reducing the chances of in-hospital infection.

# Part (b) Monitoring the odds ratio
sepsis.re.model2 <- "
model
{
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
"
variable.names <- c("p", "alpha", "mu", "delta", "delta.or")

sepsis.re.jag <- jags.model(textConnection(sepsis.re.model2),
                            data = sepsis.data,
                            inits = ini,
                            n.chains = 1)
update(sepsis.re.jag, n.iter = 1000)
sepsis.re.coda <- coda.samples(model = sepsis.re.jag,
                               variable.names = variable.names,
                               n.iter = 50000)

summary(sepsis.re.coda[, "delta.or"])

# The 95% posterior CI for the odd ratio, delta.or, is (0.26, 0.99), with median 0.58.
# Note that rather than re-running the whole model, we could simply
# exponentiate the stored samples for delta. This approach is very useful for
# complex models that take a long time to run!

# Part (c) Fixed effects model
sepsis.fe.model <- "
model
{
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
"
variable.names <- c("p", "alpha", "delta", "delta.or")
ini <- list(list(delta = c(NA, 0), alpha = rep(0, 10)))

sepsis.fe.jag <- jags.model(textConnection(sepsis.fe.model),
                            data = sepsis.data,
                            inits = ini,
                            n.chains = 1)
update(sepsis.fe.jag, n.iter = 1000)
sepsis.fe.coda <- coda.samples(model = sepsis.fe.jag,
                               variable.names = variable.names,
                               n.iter = 10000)

summary(sepsis.fe.coda[, "delta.or"])
summary(sepsis.fe.coda[, "delta[2]"])


# The 95% posterior CI for the odds ratio is now (0.68, 0.96)
# See also the results shown for each study in the lecture
# Under the fixed effect model, the lower bound of the 95% CI for the
# overall odds ratio is considerably higher.
# This is due to the lack of shrinkage in the fixed effects model
# e.g. the Sandberg (2000) study is not shrunk down towards the mean

# Question 3 - US studies
is.us <- c(1, 0, 1, 0, 1, 0, 0, 0, 0, 1)
sepsis.data$is.us <- is.us

sepsis.re.model <- "
model
{
  for (s in 1:Ns){ # for each study, Ns = total number of studies
    for (a in 1:2){ # for each arm
      y[s, a] ~ dbin(p[s, a], n[s, a])
      logit(p[s, a]) <- alpha[s] + mu[s, a]
    }
    # treatment effect
    mu[s, 1] <- 0
    mu[s, 2] ~ dnorm(delta + beta * is.us[s], prec.mu)

    ## study-specific baselines, vague priors
    alpha[s] ~ dnorm(0, 0.01)
  }

  delta ~ dnorm(0, 0.01)
  beta ~ dnorm(0, 0.01)

  prec.mu <- 1/(sd.mu * sd.mu)
  sd.mu ~ dunif(0, 10)
}
"

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
variable.names <- c("p", "alpha", "mu", "delta", "beta")

sepsis.re.jag <- jags.model(textConnection(sepsis.re.model),
                            data = sepsis.data,
                            inits = ini,
                            n.chains = 1)
update(sepsis.re.jag, n.iter = 1000)
sepsis.re.coda <- coda.samples(model = sepsis.re.jag,
                               variable.names = variable.names,
                               n.iter = 50000,
                               thin = 10)

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

# Model with vague prior
gran.re.vague.model <- "
model
{
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
"

ini <- list(
  list(
    delta = 0,
    alpha = rep(0, 6)
  )
)

variable.names <- c("p", "alpha", "mu", "delta", "delta.or")
gran.re.vague.jag <- jags.model(textConnection(gran.re.vague.model),
                                data = gran.data,
                                inits = ini,
                                n.chains = 1)
update(gran.re.vague.jag, n.iter = 1000)
gran.re.vague.coda <- coda.samples(model = gran.re.vague.jag,
                                   variable.names = variable.names,
                                   n.iter = 50000)

summary(gran.re.vague.coda)
summary(gran.re.vague.coda[, "delta"])

# With an informative prior
gran.re.scept.model <- "
model
{
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
"

ini <- list(
  list(
    delta = 0,
    alpha = rep(0, 6)
  )
)

variable.names <- c("p", "alpha", "mu", "delta", "delta.or")
gran.re.scept.jag <- jags.model(textConnection(gran.re.scept.model),
                                data = gran.data,
                                inits = ini,
                                n.chains = 1)
update(gran.re.scept.jag, n.iter = 1000)
gran.re.scept.coda <- coda.samples(model = gran.re.scept.jag,
                             variable.names = variable.names,
                             n.iter = 50000)

summary(gran.re.scept.coda)
summary(gran.re.scept.coda[, "delta"])

# Vague delta ~ Unif(-100, 100) prior
gran.re.vague2.model <- "
model
{
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
"

ini <- list(
  list(
    delta = 0,
    alpha = rep(0, 6)
  )
)

variable.names <- c("p", "alpha", "mu", "delta", "delta.or")
gran.re.vague2.jag <- jags.model(textConnection(gran.re.vague2.model),
                                data = gran.data,
                                inits = ini,
                                n.chains = 1)
update(gran.re.vague2.jag, n.iter = 1000)
gran.re.vague2.coda <- coda.samples(model = gran.re.vague2.jag,
                                   variable.names = variable.names,
                                   n.iter = 50000)

summary(gran.re.vague2.coda)
summary(gran.re.vague2.coda[, "delta"])

# Under the vague N(0, 10^2) prior delta's posterior media is
# -1.11 (95% CI -3.6, 0.76)
# Under the informative prior, posterior median -0.04 (95% CI -0.33, 0.25)
# Under the U(-100, 100) prior, delta's posterior median is -1.12 (95% CI -3.74, 0.89)
# The posterior distribution is much more concentrated around 0 under the
# informative prior, but there is no substantive difference between the
# two vague priors, although the tails are a bit longer under the Uniform prior

# Question 3
gran.re.model <- "
model
{
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
"

ini <- list(
  list(
    delta = rep(0, 4),
    alpha = matrix(rep(rep(0, 6), 4), nrow = 4)
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

variable.names <- c("p", "alpha", "mu", "delta", "delta.or", "sd.mu", "prec.mu")
gran.re.jag <- jags.model(textConnection(gran.re.model),
                          data = gran.data.edited,
                          inits = ini,
                          n.chains = 1)
update(gran.re.jag, n.iter = 1000)
gran.re.coda <- coda.samples(model = gran.re.jag,
                             variable.names = variable.names,
                             n.iter = 50000)

summary(gran.re.coda)
summary(gran.re.coda[, "delta[1]"])
summary(gran.re.coda[, "delta[2]"])
summary(gran.re.coda[, "delta[3]"])
summary(gran.re.coda[, "delta[4]"])

summary(gran.re.coda[, "p[1,1,1]"])
summary(gran.re.coda[, "p[2,1,1]"])
summary(gran.re.coda[, "p[3,1,1]"])
summary(gran.re.coda[, "p[4,1,1]"])

summary(gran.re.coda[, "p[1,4,2]"])
summary(gran.re.coda[, "p[2,4,2]"])
summary(gran.re.coda[, "p[3,4,2]"])
summary(gran.re.coda[, "p[4,4,2]"])

# Posterior quantiles under each prior
#
#                   2.5%       25%       50%      75%   97.5%
# delta[1]    -3.628e+00 -1.724845 -1.114454 -0.575547  0.79024
# delta[2]    -3.691e+00 -1.735634 -1.115782 -0.570793  0.83823
# delta[3]    -2.633e+00 -1.500165 -1.043157 -0.612633  0.29363
# delta[4]    -2.942e+00 -1.499531 -1.007866 -0.566172  0.42063
# delta.or[1]  2.656e-02  0.178201  0.328094  0.562397  2.20392
# delta.or[2]  2.494e-02  0.176288  0.327659  0.565077  2.31227
# delta.or[3]  7.189e-02  0.223093  0.352341  0.541922  1.34128
# delta.or[4]  5.277e-02  0.223235  0.364997  0.567694  1.52293
# prec.mu[1]   3.096e-02  0.162511  0.343468  0.682540  2.86032
# prec.mu[2]   2.726e-02  0.157636  0.334984  0.682145  2.86903
# prec.mu[3]   1.362e-01  0.418507  0.707247  1.147868  2.68728
# prec.mu[4]   7.138e-02  0.309834  0.641900  1.351689 16.71698
# sd.mu[1]     5.913e-01  1.210420  1.706306  2.480611  5.68349
# sd.mu[2]     5.904e-01  1.210770  1.727778  2.518677  6.05667
# sd.mu[3]     6.100e-01  0.933370  1.189089  1.545784  2.70963
# sd.mu[4]     2.446e-01  0.860125  1.248148  1.796534  3.74306
# p[1,1,1]     1.515e-01  0.301896  0.396404  0.494869  0.67474
# p[2,1,1]     1.572e-01  0.300605  0.393181  0.492128  0.67229
# p[3,1,1]     1.697e-01  0.316381  0.406242  0.501308  0.67905
# p[4,1,1]     1.668e-01  0.315173  0.406265  0.500348  0.67787
# p[1,4,2]     1.462e-05  0.001579  0.007018  0.021949  0.09899
# p[2,4,2]     9.044e-06  0.001626  0.007161  0.022080  0.09941
# p[3,4,2]     1.666e-04  0.003310  0.010437  0.026833  0.10321
# p[4,4,2]     7.585e-05  0.002838  0.009877  0.026353  0.10178

# The results under priors 1 and 2 are quite similar, but are fairly
# different to priors 3 and 4, which are fairly similar to each other,
# with 4 being more extreme.
# Note that the probabilities for each under are quite stable across all the
# priors -- the data provide a lot of information about these quantites, and
# so the prior plays a less strong part in these estimates compared to the
# random effect variance parameter.

# Question 4
gran.re.model <- "
model
{
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
"

ini <- list(
  list(
    delta = 0,
    alpha = rep(0, 6)
  )
)

variable.names <- c("p", "alpha", "mu", "delta", "delta.or", "sd.mu", "prec.mu")
gran.re.jag <- jags.model(textConnection(gran.re.model),
                          data = gran.data,
                          inits = ini,
                          n.chains = 1)
update(gran.re.jag, n.iter = 1000)
gran.re.coda <- coda.samples(model = gran.re.jag,
                             variable.names = variable.names,
                             n.iter = 100000)

summary(gran.re.coda)
summary(gran.re.coda[, "delta"])

#                2.5%       25%      50%      75%    97.5%
# delta    -1.7966980 -1.05570 -0.75246 -0.48757  -0.002604
# delta.or  0.1658456  0.34795  0.47121  0.61412   0.997400
# prec.mu   0.4958213  1.98749  5.11030 20.11280 347.118539
# sd.mu     0.0536736  0.22298  0.44236  0.70933   1.420160

## Under the informative prior, the 95% CI for delta is much narrower, since the informative prior allows much lower levels of heterogeneity than our previous vague priors, which leads to greater precision in the treatment effect. The interval also doesn't include 0, in contrast to the other priors!
