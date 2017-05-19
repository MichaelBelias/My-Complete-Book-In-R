### run samples

library(rjags)

source('model.R')
dataf <- read.csv('fakeData.csv') # or replace with own data set

# here we use sum-to-zero effects coding
options(contrasts=c('contr.sum', 'contr.sum'))

contrasts(dataf$J)
contrasts(dataf$K)

# design matrix
X <- model.matrix(~ J*K, data = dataf)[,2:4] # don't include intercept (dealt with via B0 parameter)

# put data in a list
datal <- list(
  y = dataf$y,
  n = nrow(dataf),
  X = X,
  nEff = ncol(X),
  id = dataf$ppt,
  S = length(unique(dataf$ppt))
  )

# settings
nAdapt = 1000
nBurn = 1000
nChains = 4 # as JAGS has 4 random number generators
nSave = 10^4 # may want to increase
nThin = 1
nIter = ceiling((nSave*nThin)/nChains)

# create JAGS model
mod = jags.model(textConnection(modelString), data = datal, n.chains = nChains, n.adapt = nAdapt)

# burn in
update(mod, n.iter = nBurn)

# MCMC samples
samp = coda.samples(mod, variable.names=params, n.iter=nIter, thin=nThin)

# check posterior
gelman.diag(samp)
effectiveSize(samp) # ideally should be > 10000 (see Kruschke, 2015)

# summary of posterior
summary(samp)

# extract to matrix
mcmc <- as.matrix(samp, chains = T)

