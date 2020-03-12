## Logistic regression with additive measurement error
# Logistic regression with measurement error in the independent variables (assumed to be normally distributed)
#
# Set the parameters 
number_of_sims <- 1000
beta0 <- 1 ; beta1 <- 1 ; beta2 <- 0  # When beta2=0, true model is just xi1
.Random.seed <- c(0, 324,48,6085)     # Seed for random number generation


corr <- .90                           # Correlation to be used in the correlation matrix of X 
Phi <- rbind(c(1,corr), c(corr,1))    # Covariance matrix of latent vars
reli1 <- .50                          # Reliability of xj1 (measuring xi1)
reli2 <- .95                          # Reliability of xj2 (measuring xi2)
n    <- 250                           #
sdverr1 <- sqrt((1-reli1)/reli1)      #
sdverr2 <- sqrt((1-reli2)/reli2)      #
#######################################



rmvn <- function(nn,mu,sigma)
  # Returns an nn by kk matrix, rows are independent MVN(mu,sigma)
{
  kk <- length(mu)
  dsig <- dim(sigma)
  if(dsig[1] != dsig[2]) stop("Sigma must be square.")
  if(dsig[1] != kk) stop("Sizes of sigma and mu are inconsistent.")
  ev <- eigen(sigma,symmetric=T)
  sqrl <- diag(sqrt(ev$values))
  PP <- ev$vectors
  ZZ <- rnorm(nn*kk) ; dim(ZZ) <- c(kk,nn)
  rmvn <- t(PP%*%sqrl%*%ZZ+mu)
  rmvn
}# End of function rmvn

#######################################
nsig <- 0
cat("  \n")
# Simulate some data (in loop)

for(i in 1:number_of_sims)
{
  Xi <- rmvn(n,c(0,0),Phi)                      # Generate the 
  logodds <- beta0 + beta1*Xi[,1] + beta2*Xi[,2]
  pi <- exp(logodds)/(1+exp(logodds)) 
  y <- rbinom(n,1,pi) # Binary data from true logistic regression model.
  delta1 <- rnorm(n,mean=0,sd=sdverr1)
  delta2 <- rnorm(n,mean=0,sd=sdverr2)
  x1 <- Xi[,1] + delta1 # Xi1 measured with error 
  x2 <- Xi[,2] + delta2 # Xi2 measured with error
  fullmod <- glm(y ~ x1+x2, family=binomial)
  redmod <- glm(y ~ x1, family=binomial)
  sig <- deviance(redmod) - deviance(fullmod) > 3.841
  cat("          ",sig,"\n")
  if(sig) nsig <- nsig+1
}# End nsim simulations     

nsig/number_of_sims
