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

# Question 2 - Fixed effects model

# Question 3 - US studies

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

# Model with vague Uniform prior

# Question 3

# to run all the priors in the same model file, you might find the following
# dataset and initial conditions useful - see ask or solutions if you
# are unsure!
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

# Question 4
