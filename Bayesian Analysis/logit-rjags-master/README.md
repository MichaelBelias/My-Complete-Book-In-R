## Logit model in rjags

A simple mixed effects model for Bayesian analysis of binary data (responses are 0 or 1).

* makeData.R creates a simulated data set with two repeated measures independent variables (A and B). The simulated data has two main effects (on logit scale) and no interaction with a random participant effect ~ N(0, 1). The data set is written to fakeData.csv to be read by run.R.

* model.R contains the JAGS model string and the to-be-monitored parameters.

* run.R structures data for rjags and runs mcmc samples.

See stephenrho.github.io/rjags-model.html for more detail.