## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=12, fig.height=8,
                      echo=FALSE, warning=FALSE, message=FALSE,comment = "")

## ----echo=FALSE----------------------------------------------------------
set.seed(51)
library(ggplot2)
library(knitr)
x= sample(1:6, size = 10000, replace = T)
barplot(table(x), col = "darkgreen" , axes = F)

## ------------------------------------------------------------------------
set.seed(34)

n = 100
x = numeric(n)

for (i in 2:n) {
  x[i] = rnorm(1, mean=x[i-1], sd=1.0)
}

plot.ts(x)

## ------------------------------------------------------------------------
Q = matrix(c(0.0, 0.5, 0.0, 0.0, 0.5,
             0.5, 0.0, 0.5, 0.0, 0.0,
             0.0, 0.5, 0.0, 0.5, 0.0,
             0.0, 0.0, 0.5, 0.0, 0.5,
             0.5, 0.0, 0.0, 0.5, 0.0), 
           nrow=5, byrow=TRUE)

Q %*% Q # Matrix multiplication in R. This is Q^2.

## ------------------------------------------------------------------------
Q5 = Q %*% Q %*% Q %*% Q %*% Q # h=5 steps in the future
round(Q5, 3)

## ------------------------------------------------------------------------
Q10 = Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q # h=10 steps in the future
round(Q10, 3)

## ------------------------------------------------------------------------
Q30 = Q
for (i in 2:30) {
  Q30 = Q30 %*% Q
}
round(Q30, 3) # h=30 steps in the future

## ------------------------------------------------------------------------
c(0.2, 0.2, 0.2, 0.2, 0.2) %*% Q

## ------------------------------------------------------------------------
n = 5000
x = numeric(n)
x[1] = 1 # fix the state as 1 for time 1
for (i in 2:n) {
  x[i] = sample.int(5, size=1, prob=Q[x[i-1],]) # draw the next state from the intergers 1 to 5 with probabilities from the transition matrix Q, based on the previous value of X.
}

## ------------------------------------------------------------------------
table(x) / n

## ------------------------------------------------------------------------
set.seed(38)

n = 1500
x = numeric(n)
phi = -0.6

for (i in 2:n) {
  x[i] = rnorm(1, mean=phi*x[i-1], sd=1.0)
}

plot.ts(x)

## ----echo=FALSE----------------------------------------------------------
hist(x, freq=FALSE)
curve(dnorm(x, mean=0.0, sd=sqrt(1.0/(1.0-phi^2))), col="red", add=TRUE)
legend("topright", legend="theoretical stationary\ndistribution", col="red", lty=1, bty="n")

## ------------------------------------------------------------------------
lg = function(mu, n, ybar) {
  mu2 = mu^2
  n * (ybar * mu - mu2 / 2.0) - log(1 + mu2)
}

## ------------------------------------------------------------------------
mh = function(n, ybar, n_iter, mu_init, cand_sd) {
  ## Random-Walk Metropolis-Hastings algorithm
  
  ## step 1, initialize
  mu_out = numeric(n_iter)
  accpt = 0
  mu_now = mu_init
  lg_now = lg(mu=mu_now, n=n, ybar=ybar)
  
  ## step 2, iterate
  for (i in 1:n_iter) {
    ## step 2a
    mu_cand = rnorm(n=1, mean=mu_now, sd=cand_sd) # draw a candidate
    
    ## step 2b
    lg_cand = lg(mu=mu_cand, n=n, ybar=ybar) # evaluate log of g with the candidate
    lalpha = lg_cand - lg_now # log of acceptance ratio
    alpha = exp(lalpha)
    
    ## step 2c
    u = runif(1) # draw a uniform variable which will be less than alpha with probability min(1, alpha)
    if (u < alpha) { # then accept the candidate
      mu_now = mu_cand
      accpt = accpt + 1 # to keep track of acceptance
      lg_now = lg_cand
    }
    
    ## collect results
    mu_out[i] = mu_now # save this iteration's value of mu
  }
  
  ## return a list of output
  list(mu=mu_out, accpt=accpt/n_iter)
}

## ------------------------------------------------------------------------
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
ybar = mean(y)
n = length(y)
hist(y, freq=FALSE, xlim=c(-1.0, 3.0)) # histogram of the data
curve(dt(x=x, df=1), lty=2, add=TRUE) # prior for mu
points(y, rep(0,n), pch=1) # individual data points
points(ybar, 0, pch=19) # sample mean

## ----comment=""----------------------------------------------------------
set.seed(43) # set the random seed for reproducibility
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=3.0)
str(post)

library("coda")
traceplot(as.mcmc(post$mu))

## ------------------------------------------------------------------------
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=0.05)
post$accpt
traceplot(as.mcmc(post$mu))

post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=0.9)
post$accpt

traceplot(as.mcmc(post$mu))

## ------------------------------------------------------------------------
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=30.0, cand_sd=0.9)
post$accpt

## ------------------------------------------------------------------------
traceplot(as.mcmc(post$mu))

## ------------------------------------------------------------------------
post$mu_keep = post$mu[-c(1:100)] # discard the first 200 samples
plot(density(post$mu_keep, adjust=2.0), main="", xlim=c(-1.0, 3.0), xlab=expression(mu)) # plot density estimate of the posterior
curve(dt(x=x, df=1), lty=2, add=TRUE) # prior for mu
points(ybar, 0, pch=19) # sample mean

curve(0.017*exp(lg(mu=x, n=n, ybar=ybar)), from=-1.0, to=3.0, add=TRUE, col="blue") # approximation to the true posterior in blue

## ----echo=FALSE----------------------------------------------------------
update_mu = function(n, ybar, sig2, mu_0, sig2_0) {
  sig2_1 = 1.0 / (n / sig2 + 1.0 / sig2_0)
  mu_1 = sig2_1 * (n * ybar / sig2 + mu_0 / sig2_0)
  rnorm(n=1, mean=mu_1, sd=sqrt(sig2_1))
}

## ----echo=FALSE----------------------------------------------------------
update_sig2 = function(n, y, mu, nu_0, beta_0) {
  nu_1 = nu_0 + n / 2.0
  sumsq = sum( (y - mu)^2 ) # vectorized
  beta_1 = beta_0 + sumsq / 2.0
  out_gamma = rgamma(n=1, shape=nu_1, rate=beta_1) # rate for gamma is shape for inv-gamma
  1.0 / out_gamma # reciprocal of a gamma random variable is distributed inv-gamma
}

## ----echo=FALSE----------------------------------------------------------
gibbs = function(y, n_iter, init, prior) {
  ybar = mean(y)
  n = length(y)
  
  ## initialize
  mu_out = numeric(n_iter)
  sig2_out = numeric(n_iter)
  
  mu_now = init$mu
  
  ## Gibbs sampler
  for (i in 1:n_iter) {
    sig2_now = update_sig2(n=n, y=y, mu=mu_now, nu_0=prior$nu_0, beta_0=prior$beta_0)
    mu_now = update_mu(n=n, ybar=ybar, sig2=sig2_now, mu_0=prior$mu_0, sig2_0=prior$sig2_0)
    
    sig2_out[i] = sig2_now
    mu_out[i] = mu_now
  }
  
  cbind(mu=mu_out, sig2=sig2_out)
}

## ----echo=FALSE----------------------------------------------------------
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
ybar = mean(y)
n = length(y)

## prior
prior = list()
prior$mu_0 = 0.0
prior$sig2_0 = 1.0
prior$n_0 = 2.0 # prior effective sample size for sig2
prior$s2_0 = 1.0 # prior point estimate for sig2
prior$nu_0 = prior$n_0 / 2.0 # prior parameter for inverse-gamma
prior$beta_0 = prior$n_0 * prior$s2_0 / 2.0 # prior parameter for inverse-gamma

hist(y, freq=FALSE, xlim=c(-1.0, 3.0)) # histogram of the data
curve(dnorm(x=x, mean=prior$mu_0, sd=sqrt(prior$sig2_0)), lty=2, add=TRUE) # prior for mu
points(y, rep(0,n), pch=1) # individual data points
points(ybar, 0, pch=19) # sample mean

## ----comment="",echo=FALSE-----------------------------------------------
set.seed(53)

init = list()
init$mu = 0.0

post = gibbs(y=y, n_iter=1e3, init=init, prior=prior)
head(post)

library("coda")
plot(as.mcmc(post))


summary(as.mcmc(post))

## ------------------------------------------------------------------------
data("PlantGrowth")

head(PlantGrowth)


## ------------------------------------------------------------------------
boxplot(weight ~ group, data=PlantGrowth)

## ----comment="",echo=FALSE-----------------------------------------------
lmod = lm(weight ~ group, data=PlantGrowth)
summary(lmod)

par(mfrow=c(2,2))
plot(lmod) # for graphical residual analysis

## ----comment="",echo=FALSE-----------------------------------------------
anova(lmod)

## ----comment="",echo=FALSE-----------------------------------------------

library("rjags")

mod_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(mu[grp[i]], prec)
    }
    
    for (j in 1:3) {
        mu[j] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(5/2.0, 5*1.0/2.0)
    sig = sqrt( 1.0 / prec )
} "

set.seed(82)

data_jags = list(y=PlantGrowth$weight, 
              grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
    inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains

## ----comment="",echo=FALSE,fig.height=9----------------------------------
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
effectiveSize(mod_sim)


(pm_params = colMeans(mod_csim))

## ----comment="",echo=FALSE,fig.height=9----------------------------------
yhat = pm_params[1:3][data_jags$grp]
resid = data_jags$y - yhat
plot(resid)

## ----comment="",echo=FALSE,fig.height=9----------------------------------

plot(yhat, resid)

## ----comment="",echo=FALSE,fig.height=9----------------------------------
summary(mod_sim)

## ----comment="",echo=FALSE,fig.height=9----------------------------------
HPDinterval(mod_csim)

## ----comment="",echo=FALSE,fig.height=9----------------------------------
mean(mod_csim[,3] > mod_csim[,1])

## ----comment="",echo=FALSE,fig.height=9----------------------------------
mean(mod_csim[,3] > 1.1*mod_csim[,1])

## ----echo=FALSE, comment=""----------------------------------------------
data("warpbreaks")

head(warpbreaks)

## ----echo=FALSE, comment=""----------------------------------------------
table(warpbreaks$wool, warpbreaks$tension)

boxplot(breaks ~ wool + tension, data=warpbreaks)


## ----echo=FALSE, comment=""----------------------------------------------
boxplot(log(breaks) ~ wool + tension, data=warpbreaks)

## ---- echo=FALSE,comment="", warning=FALSE,message=FALSE,fig.height=9, echo=FALSE----
library("rjags")
mod1_string = " model {
    for( i in 1:length(y)) {
        y[i] ~ dnorm(mu[tensGrp[i]], prec)
    }
    
    for (j in 1:3) {
        mu[j] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(5/2.0, 5*2.0/2.0)
    sig = sqrt(1.0 / prec)
} "

set.seed(83)
str(warpbreaks)

data1_jags = list(y=log(warpbreaks$breaks), tensGrp=as.numeric(warpbreaks$tension))

params1 = c("mu", "sig")

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5e3)

## convergence diagnostics
plot(mod1_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)

## ---- echo=FALSE,comment="", warning=FALSE,message=FALSE, echo=FALSE-----
summary(mod1_sim)

## ---- echo=FALSE,comment="", warning=FALSE,message=FALSE, echo=FALSE-----
dic1 = dic.samples(mod1, n.iter=1e3)

## ---- echo=FALSE,comment="", warning=FALSE,message=FALSE, echo=FALSE-----
X = model.matrix(lm(breaks ~ wool + tension, data=warpbreaks))
head(X)

## ---- echo=FALSE,comment="", warning=FALSE,message=FALSE, echo=FALSE-----
tail(X)

## ---- echo=FALSE,comment="", warning=FALSE,message=FALSE, echo=FALSE-----
mod2_string = " model {
    for( i in 1:length(y)) {
        y[i] ~ dnorm(mu[i], prec)
        mu[i] = int + alpha*isWoolB[i] + beta[1]*isTensionM[i] + beta[2]*isTensionH[i]
    }
    
    int ~ dnorm(0.0, 1.0/1.0e6)
    alpha ~ dnorm(0.0, 1.0/1.0e6)
    for (j in 1:2) {
        beta[j] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(3/2.0, 3*1.0/2.0)
    sig = sqrt(1.0 / prec)
} "

data2_jags = list(y=log(warpbreaks$breaks), isWoolB=X[,"woolB"], isTensionM=X[,"tensionM"], isTensionH=X[,"tensionH"])

params2 = c("int", "alpha", "beta", "sig")

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)

## convergene diagnostics
plot(mod2_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)

## ---- echo=FALSE,comment="", warning=FALSE,message=FALSE, echo=FALSE-----
summary(mod2_sim)
(dic2 = dic.samples(mod2, n.iter=1e3))
dic1

## ---- echo=FALSE,comment="", warning=FALSE,message=FALSE, echo=FALSE-----
boxplot(log(breaks) ~ wool + tension, data=warpbreaks)

## ---- echo=FALSE,comment="", warning=FALSE,message=FALSE, echo=FALSE-----
lmod2 = lm(log(breaks) ~ .^2, data=warpbreaks)
summary(lmod2)

## ---- echo=FALSE,fig.height=9,comment="", warning=FALSE,message=FALSE, echo=FALSE----

mod3_string = " model {
    for( i in 1:length(y)) {
        y[i] ~ dnorm(mu[woolGrp[i], tensGrp[i]], prec)
    }
    
    for (j in 1:max(woolGrp)) {
        for (k in 1:max(tensGrp)) {
            mu[j,k] ~ dnorm(0.0, 1.0/1.0e6)
        }
    }
    
    prec ~ dgamma(3/2.0, 3*1.0/2.0)
    sig = sqrt(1.0 / prec)
} "

str(warpbreaks)

data3_jags = list(y=log(warpbreaks$breaks), woolGrp=as.numeric(warpbreaks$wool), tensGrp=as.numeric(warpbreaks$tension))

params3 = c("mu", "sig")

mod3 = jags.model(textConnection(mod3_string), data=data3_jags, n.chains=3)
update(mod3, 1e3)

mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=5e3)
mod3_csim = as.mcmc(do.call(rbind, mod3_sim))

plot(mod3_sim, ask=TRUE)

## convergence diagnostics
gelman.diag(mod3_sim)
autocorr.diag(mod3_sim)
effectiveSize(mod3_sim)
raftery.diag(mod3_sim)

## ---- echo=FALSE, warning=FALSE,message=FALSE, echo=FALSE----------------
(dic3 = dic.samples(mod3, n.iter=1e3))
dic2
dic1

## ---- echo=FALSE, warning=FALSE,message=FALSE, echo=FALSE----------------
summary(mod3_sim)
HPDinterval(mod3_csim)
par(mfrow=c(3,2)) # arrange frame for plots
densplot(mod3_csim[,1:6], xlim=c(2.0, 4.5))

## ---- echo=FALSE, warning=FALSE,message=FALSE, echo=FALSE----------------
prop.table( table( apply(mod3_csim[,1:6], 1, which.min) ) )

## ----echo=FALSE, comment=""----------------------------------------------
library("car")
Lein=car::Leinhardt
head(Lein)

## ----echo=FALSE, comment=""----------------------------------------------
str(Lein)

## ----echo=FALSE, comment=""----------------------------------------------
pairs(Lein)

## ----echo=FALSE, comment=""----------------------------------------------
plot(infant ~ income, data=Lein)

## ----echo=FALSE, comment=""----------------------------------------------
hist(Leinhardt$infant, main = "Histogram of infant mortality")
hist(Leinhardt$income, main = "Histogram of income")

## ----echo=FALSE, comment=""----------------------------------------------
Leinhardt$loginfant = log(Leinhardt$infant)
Leinhardt$logincome = log(Leinhardt$income)

plot(loginfant ~ logincome, data=Leinhardt)

## ----echo=FALSE, comment=""----------------------------------------------
lmod = lm(loginfant ~ logincome, data=Leinhardt)
summary(lmod)


## ----echo=FALSE, comment="",warning=FALSE, fig.height=8------------------
par(mfrow=c(2,2))
plot(lmod)
dev.off()

## ---- echo=FALSE, warning=FALSE,message=FALSE, echo=FALSE----------------
dat = na.omit(Leinhardt)



mod1_string = " model {
    for (i in 1:n) {
        y[i] ~ dnorm(mu[i], prec) ### y comes from a normal distribution with mu and prec
        mu[i] = b[1] + b[2]*log_income[i] ## The mu is calculated as a linear prediction
    }
    
    for (i in 1:2) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)  ## Both coefficients have uninformative priors
    }
    
    prec ~ dgamma(5/2.0, 5*10.0/2.0) ## precision is modeled from a gamma distribution
    sig2 = 1.0 / prec                ## So variance is from inverse gamma
    sig = sqrt(sig2)                 ## and the standard deviation too 
} "

set.seed(72)
data1_jags = list(y=dat$loginfant, n=nrow(dat), 
              log_income=dat$logincome)

params1 = c("b", "sig")

inits1 = function() {
    inits = list("b"=rnorm(2,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, inits=inits1, n.chains=3)
update(mod1, 1000) # burn-in

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5000)

mod1_csim = do.call(rbind, mod1_sim) # combine multiple chains

## ---- echo=FALSE, warning=FALSE,message=FALSE, echo=FALSE, comment=""----
plot(mod1_sim)
gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

summary(mod1_sim)

## ---- echo=FALSE, warning=FALSE,message=FALSE, echo=FALSE, comment="", fig.height=7----
lmod0 = lm(infant ~ income, data=Leinhardt)
plot(resid(lmod0)) # to check independence (looks okay)


plot(predict(lmod0), resid(lmod0)) # to check for linearity, constant variance (looks bad)

qqnorm(resid(lmod0)) # to check Normality assumption (we want this to be a straight line)

## ---- echo=FALSE, warning=FALSE,message=FALSE, echo=FALSE, comment=""----
X = cbind(rep(1.0, data1_jags$n), data1_jags$log_income)
head(X)
(pm_params1 = colMeans(mod1_csim)) # posterior mean


yhat1 = drop(X %*% pm_params1[1:2])
resid1 = data1_jags$y - yhat1
plot(resid1) # against data index


plot(yhat1, resid1) # against predicted values

qqnorm(resid1) # checking normality of residuals

plot(predict(lmod), resid(lmod)) # to compare with reference linear model

rownames(dat)[order(resid1, decreasing=TRUE)[1:5]] # which countries have the largest positive residuals?

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
library("rjags")
set.seed(75)
mod2_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(mu[i], prec)
        mu[i] = b[1] + b[2]*log_income[i] + b[3]*is_oil[i]
    }
    
    for (i in 1:3) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(5/2.0, 5*10.0/2.0)
    sig = sqrt( 1.0 / prec )
} "


set.seed(75)
data2_jags = list(y=dat$loginfant, log_income=dat$logincome,
                  is_oil=as.numeric(dat$oil=="yes"))
data2_jags$is_oil

params2 = c("b", "sig")

inits2 = function() {
    inits = list("b"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, inits=inits2, n.chains=3)
update(mod2, 1e3) # burn-in

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)

mod2_csim = as.mcmc(do.call(rbind, mod2_sim)) # combine multiple chains

## ----echo=FALSE,fig.height=9, warning=FALSE, message=FALSE, comment=""----
plot(mod2_sim)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)

autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)
summary(mod2_sim)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
X2 = cbind(rep(1.0, data1_jags$n), data2_jags$log_income, data2_jags$is_oil)
head(X2)


## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
(pm_params2 = colMeans(mod2_csim)) # posterior mean

yhat2 = drop(X2 %*% pm_params2[1:3])
resid2 = data2_jags$y - yhat2
plot(resid2) # against data index

## ----echo=FALSE, warning=FALSE, message=FALSE, comment="", fig.height=9----
par(mfrow=c(2,1))
plot(yhat1, resid1, main = "Residuals from the first model") # residuals from the first model
plot(yhat2, resid2, main= "Residuals from the second model") # against predicted values

sd(resid2) # standard deviation of residuals

## ----echo=FALSE, warning=FALSE, message=FALSE, comment="", fig.height=9----
mod3_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dt( mu[i], tau, df )
        mu[i] = b[1] + b[2]*log_income[i] + b[3]*is_oil
    }
    
    for (i in 1:3) {
        b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    df = nu + 2.0 # we want degrees of freedom > 2 to guarantee existence of mean and variance
    nu ~ dexp(1.0)
    
    tau ~ dgamma(5/2.0, 5*10.0/2.0) # tau is close to, but not equal to the precision
    sig = sqrt( 1.0 / tau * df / (df - 2.0) ) # standard deviation of errors
} "



## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
dic.samples(mod1, n.iter=1e3)
dic.samples(mod2, n.iter=1e3)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
library("boot")
data("urine")
dat = na.omit(urine)
head(urine)


## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
pairs(dat)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
library("corrplot")
Cor = cor(dat)
corrplot(Cor, type="upper", method="ellipse", tl.pos="d")
corrplot(Cor, type="lower", method="number", col="black", 
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
X = scale(dat[,-1], center=TRUE, scale=TRUE)
head(X[,"gravity"])


## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
colMeans(X)
apply(X, 2, sd)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
ddexp = function(x, mu, tau) {
  0.5*tau*exp(-tau*abs(x-mu)) 
}
curve(ddexp(x, mu=0.0, tau=1.0), from=-5.0, to=5.0, ylab="density", main="Double exponential\ndistribution") # double exponential distribution
curve(dnorm(x, mean=0.0, sd=1.0), from=-5.0, to=5.0, lty=2, add=TRUE) # normal distribution
legend("topright", legend=c("double exponential", "normal"), lty=c(1,2), bty="n")

## ----echo=FALSE,fig.height=9, warning=FALSE, message=FALSE, comment=""----
library("rjags")
mod1_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*ph[i] + b[3]*osmo[i] + b[4]*cond[i] + b[5]*urea[i] + b[6]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:6) {
        b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
    }
} "

set.seed(92)
head(X)

data_jags = list(y=dat$r, gravity=X[,"gravity"], ph=X[,"ph"], osmo=X[,"osmo"], cond=X[,"cond"], urea=X[,"urea"], calc=X[,"calc"])

params = c("int", "b")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

## convergence diagnostics
plot(mod1_sim, ask=TRUE)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
summary(mod1_sim)
par(mfrow=c(3,2))
densplot(mod1_csim[,1:6], xlim=c(-3.0, 3.0))

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
colnames(X) # variable names

## ----echo=FALSE, warning=FALSE,fig.height=9, message=FALSE, comment=""----
mod2_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*cond[i] + b[3]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:3) {
        b[j] ~ dnorm(0.0, 1.0/25.0) # noninformative for logistic regression
    }
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim, ask=TRUE)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

dic2 = dic.samples(mod2, n.iter=1e3)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
dic1
dic2
summary(mod2_sim)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
HPDinterval(mod2_csim)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
par(mfrow=c(3,1))
densplot(mod2_csim[,1:3], xlim=c(-3.0, 3.0))
colnames(X)[c(1,4,6)] # variable names

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
(pm_coef = colMeans(mod2_csim))

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
pm_Xb = pm_coef["int"] + X[,c(1,4,6)] %*% pm_coef[1:3]
phat = 1.0 / (1.0 + exp(-pm_Xb))
head(phat)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
plot(phat, jitter(dat$r))

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
(tab0.5 = table(phat > 0.5, data_jags$y))

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
sum(diag(tab0.5)) / sum(tab0.5)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
(tab0.3 = table(phat > 0.3, data_jags$y))
sum(diag(tab0.3)) / sum(tab0.3)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
library("COUNT")
data("badhealth")
head(badhealth)
any(is.na(badhealth))

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
hist(badhealth$numvisit, breaks=20)

## ------------------------------------------------------------------------
plot(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==0, xlab="age", ylab="log(visits)")
points(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==1, col="red")

## ------------------------------------------------------------------------
library("rjags")
mod_string = " model {
    for (i in 1:length(numvisit)) {
        numvisit[i] ~ dpois(lam[i])
        log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
    }
    
    int ~ dnorm(0.0, 1.0/1e6)
    b_badh ~ dnorm(0.0, 1.0/1e4)
    b_age ~ dnorm(0.0, 1.0/1e4)
    b_intx ~ dnorm(0.0, 1.0/1e4)
} "

set.seed(102)

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)

## ------------------------------------------------------------------------
X = as.matrix(badhealth[,-1])
X = cbind(X, with(badhealth, badh*age))
head(X)

## ------------------------------------------------------------------------
(pmed_coef = apply(mod_csim, 2, median))

## ------------------------------------------------------------------------
llam_hat = pmed_coef["int"] + X %*% pmed_coef[c("b_badh", "b_age", "b_intx")]
lam_hat = exp(llam_hat)

hist(lam_hat)

## ------------------------------------------------------------------------
resid = badhealth$numvisit - lam_hat
plot(resid) # the data were ordered

## ------------------------------------------------------------------------
plot(lam_hat, badhealth$numvisit)
abline(0.0, 1.0)

## ------------------------------------------------------------------------
plot(lam_hat[which(badhealth$badh==0)], resid[which(badhealth$badh==0)], xlim=c(0, 8), ylab="residuals", xlab=expression(hat(lambda)), ylim=range(resid))
points(lam_hat[which(badhealth$badh==1)], resid[which(badhealth$badh==1)], col="red")

## ------------------------------------------------------------------------
var(resid[which(badhealth$badh==0)])

## ------------------------------------------------------------------------
var(resid[which(badhealth$badh==1)])

## ------------------------------------------------------------------------
summary(mod_sim)

## ------------------------------------------------------------------------
x1 = c(0, 35, 0) # good health
x2 = c(1, 35, 35) # bad health

## ------------------------------------------------------------------------
head(mod_csim)

## ------------------------------------------------------------------------
loglam1 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x1
loglam2 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x2

## ------------------------------------------------------------------------
lam1 = exp(loglam1)
lam2 = exp(loglam2)

## ------------------------------------------------------------------------
(n_sim = length(lam1))

## ------------------------------------------------------------------------
y1 = rpois(n=n_sim, lambda=lam1)
y2 = rpois(n=n_sim, lambda=lam2)

plot(table(factor(y1, levels=0:18))/n_sim, pch=2, ylab="posterior prob.", xlab="visits")
points(table(y2+0.1)/n_sim, col="red")

## ------------------------------------------------------------------------
mean(y2 > y1)

## ----echo=F--------------------------------------------------------------
dat = read.table(file="Data/cookies.dat", header=TRUE)
kable(head(dat,10), caption = "First 10 values" ) 

## ----echo=F--------------------------------------------------------------
h<-   as.matrix(table(dat$location))
kable(t(h),caption =  "number of cookies per location")

## ----echo=F--------------------------------------------------------------
hist(dat$chips, main = "Histogram of  chocolate chips in total")

## ----echo=F--------------------------------------------------------------
boxplot(chips ~ location, data=dat,main="Boxplot of Cookie production")

## ----echo=F--------------------------------------------------------------
set.seed(112)
n_sim = 500
alpha_pri = rexp(n_sim, rate=1.0/2.0)
beta_pri = rexp(n_sim, rate=5.0)
mu_pri = alpha_pri/beta_pri
sig_pri = sqrt(alpha_pri/beta_pri^2)

summary(mu_pri)

## ----echo=F--------------------------------------------------------------
summary(sig_pri)

## ----echo=F--------------------------------------------------------------
lam_pri = rgamma(n=n_sim, shape=alpha_pri, rate=beta_pri)
summary(lam_pri)

## ----echo=F, comment=""--------------------------------------------------
(lam_pri = rgamma(n=5, shape=alpha_pri[1:5], rate=beta_pri[1:5]))

## ----echo=F, comment=""--------------------------------------------------
(y_pri = rpois(n=150, lambda=rep(lam_pri, each=30)))

## ----echo=F, message=FALSE, warning=FALSE--------------------------------
library("rjags")

mod_string = " model {
for (i in 1:length(chips)) {
  chips[i] ~ dpois(lam[location[i]])
}

for (j in 1:max(location)) {
  lam[j] ~ dgamma(alpha, beta)
}

mu ~ dgamma(2.0, 1.0/5.0) ###Gamma mean = a/b
sig ~ dexp(1.0)           ### Gamma sig = a/b^2 


## we can solve the above to get 

alpha = mu^2 / sig^2 
beta = mu / sig^2   


} "

## ----echo=F, message=FALSE, comment=""-----------------------------------
set.seed(113)

data_jags = as.list(dat)

params = c("lam", "mu", "sig")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))



par(mar = rep(3, 4))
## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)

## ---- comment=""---------------------------------------------------------
## observation level residuals
(pm_params = colMeans(mod_csim))

## ---- comment=""---------------------------------------------------------
yhat = rep(pm_params[1:5], each=30)
resid = dat$chips - yhat
plot(resid)

## ---- comment=""---------------------------------------------------------
plot(jitter(yhat), resid)

## ---- comment=""---------------------------------------------------------
var(resid[yhat<7])

## ---- comment=""---------------------------------------------------------
var(resid[yhat>11])

## ---- comment=""---------------------------------------------------------
## location level residuals
lam_resid = pm_params[1:5] - pm_params["mu"]
plot(lam_resid)
abline(h=0, lty=2)

## ----comment=""----------------------------------------------------------
summary(mod_sim)

## ----echo=FALSE,comment=""-----------------------------------------------
library("car")
data("Leinhardt")

str(Leinhardt)

pairs(Leinhardt)

head(Leinhardt)

## ----echo=FALSE,comment=""-----------------------------------------------
dat = na.omit(Leinhardt)
dat$logincome = log(dat$income)
dat$loginfant = log(dat$infant)
str(dat)

## ----echo=FALSE, warning=FALSE, message=FALSE, comment=""----------------
library("rjags")

mod_string = " model {
  for (i in 1:length(y)) {
    y[i] ~ dnorm(mu[i], prec)
    mu[i] = a[region[i]] + b[1]*log_income[i] + b[2]*is_oil[i]
  }
  
  for (j in 1:max(region)) {
    a[j] ~ dnorm(a0, prec_a)
  }
  
  a0 ~ dnorm(0.0, 1.0/1.0e6)
  prec_a ~ dgamma(1/2.0, 1*10.0/2.0)
  tau = sqrt( 1.0 / prec_a )
  
  for (j in 1:2) {
    b[j] ~ dnorm(0.0, 1.0/1.0e6)
  }
  
  prec ~ dgamma(5/2.0, 5*10.0/2.0)
  sig = sqrt( 1.0 / prec )
} "

set.seed(116)
data_jags = list(y=dat$loginfant, log_income=dat$logincome,
                  is_oil=as.numeric(dat$oil=="yes"), region=as.numeric(dat$region))
data_jags$is_oil
table(data_jags$is_oil, data_jags$region)

params = c("a0", "a", "b", "sig", "tau")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3) # burn-in

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)

mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combine multiple chains

par(mar = rep(3, 4))
## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## ------------------------------------------------------------------------
library("COUNT")
library("rjags")

## ------------------------------------------------------------------------
data("badhealth")

mod_string = " model {
    for (i in 1:length(numvisit)) {
        numvisit[i] ~ dpois(lam[i])
        log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
    }
    
    int ~ dnorm(0.0, 1.0/1e6)
    b_badh ~ dnorm(0.0, 1.0/1e4)
    b_age ~ dnorm(0.0, 1.0/1e4)
    b_intx ~ dnorm(0.0, 1.0/1e4)
} "

set.seed(102)

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))
plot(density(mod_csim[,"b_badh"]), main = "Density plot for bad health coefficient" )

## ------------------------------------------------------------------------
mod2_string = " model {
    for (i in 1:length(numvisit)) {
        numvisit[i] ~ dpois(lam[i])
        log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
    }
    
    int ~ dnorm(0.0, 1.0/1e6)
    b_badh ~ dnorm(0.0, 1.0/0.2^2)
    b_age ~ dnorm(0.0, 1.0/1e4)
    b_intx ~ dnorm(0.0, 1.0/0.01^2)
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

## ------------------------------------------------------------------------
curve(dnorm(x, mean=0.0, sd=sqrt(1e4)), from=-3.0, to=3.0, ylim=c(0.0, 3.0), lty=2,
      main="b_badh", ylab="density", xlab="b_badh")
curve(dnorm(x, mean=0.0, sd=0.2), from=-3.0, to=3.0, col="red", lty=2, add=TRUE)
lines(density(mod_csim[,"b_badh"]))
lines(density(mod2_csim[,"b_badh"]), col="red")
legend("topleft", legend=c("noninformative prior", "posterior", "skeptical prior", "posterior"),
       lty=c(2,1,2,1), col=rep(c("black", "red"), each=2), bty="n")

## ------------------------------------------------------------------------
curve(dnorm(x, mean=0.0, sd=sqrt(1e4)), from=-0.05, to=0.05, ylim=c(0.0, 140.0), lty=2,
      main="b_intx", ylab="density", xlab="b_intx")
curve(dnorm(x, mean=0.0, sd=0.01), from=-0.05, to=0.05, col="red", lty=2, add=TRUE)
lines(density(mod_csim[,"b_intx"]))
lines(density(mod2_csim[,"b_intx"]), col="red")
legend("topleft", legend=c("noninformative prior", "posterior", "skeptical prior", "posterior"),
       lty=c(2,1,2,1), col=rep(c("black", "red"), each=2), bty="n")

## ------------------------------------------------------------------------
mean(mod2_csim[,"b_intx"] > 0) # posterior probability that b_intx is positive

