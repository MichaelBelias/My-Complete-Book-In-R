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

summary(sepsis.re.coda)
summary(sepsis.re.coda[, "delta"])

# The 95% posterior CI for the log odds ratio, delta, is (-1.35, -0.01),
# suggesting that IVIG is reducing the chances of in-hospital infection.

# Part (b) Monitoring the odds ratio

summary(sepsis.re.coda[, "delta"])
t=as.numeric(sepsis.re.coda[, "delta"][[1]] )
b=exp(t)

# Hierarchical Models
# Question 2 

sepsis.fe.model <- "
model
{
  for (s in 1:Ns){ # for each study, Ns = total number of studies
  for (a in 1:2){ 
  # for each arm
  y[s, a] ~ dbin(p[s, a], n[s, a])
  logit(p[s, a]) <- alpha[s] + delta[a]
  }

  ## study-specific baselines

  alpha[s] ~ dnorm(0, 0.01)
}

  delta[1] <- 0
  delta[2] ~ dnorm(0, 0.01)
}
"

## When a variable is fixed then we should put NA in the initial values

ini <- list(list(delta = c(NA,0), alpha = rep(0, 10)))

sepsis.fe.jag <- jags.model(textConnection(sepsis.fe.model),
                            data = sepsis.data,
                            inits = ini,
                            n.chains = 1)

update(sepsis.fe.jag, n.iter = 1000)

sepsis.fe.coda <- coda.samples(model = sepsis.fe.jag,
                               variable.names = c("p", "alpha", "mu", "delta"),
                               n.iter = 50000)

summary(sepsis.fe.coda)
summary(sepsis.fe.coda[, "delta[2]"])


# Question 3 
# Now we are trying to check whether there is evidence that the treatment effects 
# are higher in the US.

US <- c(1, 0, 1, 0, 1, 0, 0, 0, 0, 1)
sepsis.data$US<- US

US.sepsis.re.model <- "
model
{
  for (s in 1:Ns){ # for each study, Ns = total number of studies
  for (a in 1:2){ # for each arm
  y[s, a] ~ dbin(p[s, a], n[s, a])
  logit(p[s, a]) <- alpha[s] + mu[s, a]
  }
  # treatment effect
  mu[s, 1] <- 0
  mu[s, 2] ~ dnorm(delta + beta * US[s], prec.mu)
  
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

US.sepsis.re.jag <- jags.model(textConnection(US.sepsis.re.model),
                            data = sepsis.data,
                            inits = ini,
                            n.chains = 1)
update(US.sepsis.re.jag, n.iter = 1000)
US.sepsis.re.coda <- coda.samples(model = US.sepsis.re.jag,
                               variable.names = variable.names,
                               n.iter = 50000,
                               thin = 10)

summary(US.sepsis.re.coda)
summary(US.sepsis.re.coda[, "beta"])



# Question 3 - US studies

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

# Vague delta ~ Unif(-100, 100) prior

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
