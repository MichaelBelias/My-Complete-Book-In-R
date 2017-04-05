#
# Read the data from txt file
 simex1<-read.table('ex1.dat')
 names(simex1) <- c('y', paste('x',1:15,sep=''))
 simex1[1:5,]

#
#
full.model.mle <- lm( y~., data=simex1 )
summary(full.model.mle)
#

#
# function for the posterior estimates/summaries 
lm.bayes <- function( y, x, prior.mean=NULL, V=NULL, a=0.001, b=0.001, q=0.05, add.constant=TRUE, MLE=TRUE, digits=3 ){ 
	n<-length(y)
	if (add.constant) x<-cbind(1,x)
	p <- ncol(x)
	if (is.null(prior.mean)) prior.mean <- rep(0,p)
	if (is.null(V)) V <- 1000*diag(p)
	mu <- prior.mean
	inv.V <- solve(V)		
	inv.tilde.Sigma <- t(x) %*% x + inv.V
	tilde.Sigma <- solve( inv.tilde.Sigma )
	M <- t(x) %*% y + inv.V %*% mu
	tilde.beta <- as.vector(tilde.Sigma %*% M) 
	SS <- as.vector( t(y) %*% y - t(tilde.beta)%*%inv.tilde.Sigma%*%tilde.beta + t(mu)%*%inv.V%*%mu )
	tilde.a <- 0.5*n + a
	tilde.b <- 0.5*SS + b

	post.params <- list( beta=tilde.beta, S=tilde.Sigma, a=tilde.a, b=tilde.b )

	post.s2 <-list()

	post.s2$mean <- tilde.b/(tilde.a-1)
	post.s2$mode <- tilde.b/(tilde.a+1)
	post.s2$sd   <- post.s2$mean * sqrt( 1/(tilde.a-2) )
	
	post.beta<-matrix( nrow=ncol(x), ncol=4 )
	post.beta[,1]<- tilde.beta
	tilde.Sigmajj <- diag(tilde.Sigma)	
	post.beta[,2]<- sqrt( post.s2$mean*tilde.Sigmajj )
	se <- sqrt( (tilde.b/tilde.a)*tilde.Sigmajj )
	post.beta[,3]<- tilde.beta+qt(q/2,2*tilde.a)*se        
	post.beta[,4]<- tilde.beta+ qt((1-q)/2, 2*tilde.a )*se
	

	rownames(post.beta) <- colnames(x)
	colnames(post.beta) <- c( 'Mean', 'S.D.', paste(round(100*q/2,1),'%Q',sep=''), paste(round(100*(1-q)/2,1),'%Q',sep='') )

	if (MLE){ 
		mles <- lm( y~x-1 )
		post.beta <- cbind(post.beta, as.matrix(summary(full.model.mle)$coef[,1:2]))
		colnames(post.beta)[5:6]<-c('MLE', 'S.E.')
	}
	
	print (round(post.beta, digits))
	return( list( beta=post.beta, s2=post.s2, post.params=post.params) )
}

res<-lm.bayes( simex1$y, as.matrix(simex1[,-1] ) )
#------------------------------------------------
j<-5
npoints <- 100
b.min<- res$beta[j,1]-4*res$beta[j,2]
b.max<- res$beta[j,1]+4*res$beta[j,2]
bj <- seq( b.min, b.max, length.out=npoints ) 
a<-res$post.params$a
b<-res$post.params$b
beta<-res$post.params$beta[j]
se <- sqrt( (b/a)* diag(res$post.params$S)[j])
density.bj<- dt( bj/se, 2*a, ncp=beta/se)/se

plot( bj, density.bj, type='l')

mlej <- res$beta[j,5]
sej  <- res$beta[j,6]
lines( bj, dnorm(bj, mlej, sej), col=2, lty=2 )

#------------------------------------------------

dinvgamma <-function( x, a, b, log=FALSE ){ 
	res <- a*log(b)-lgamma(a) -(a+1)*log(x) - b/x
	if (!log) res<-exp(res)
	return(res)
}

q <- 0.001 
a<-res$post.params$a
b<-res$post.params$b
s2.min <- 1/qgamma( (1-q/2), a, b )
s2.max <- 1/qgamma( (q/2), a, b )
npoints<-100
s2val <- seq( s2.min, s2.max, length.out=npoints ) 
s2dens<- dinvgamma( s2val, a, b )
plot(s2val,s2dens,type='l')

#------------------------------------------------

library(arm)
bayesarm.model<-bayesglm(y~., data=simex1, prior.scale=1000, prior.df=Inf, n.iter=3000)
bayesarm.model
summary(bayesarm.model)
names(bayesarm.model)

mle <- lm(y~., data=simex1)
mcmc.b <- coef( sim(mle, n.sims=3000) )
dim(mcmc.b)
head(mcmc.b)

par(mfrow=c(1,2))
hist(mcmc.b[,6])
plot(density(mcmc.b[,6]))

mcmc.s <- sigma.hat( sim(mle, n.sims=3000) )
dim(mcmc.s)
head(mcmc.s)

par(mfrow=c(1,2))
hist(mcmc.s)
plot(density(mcmc.s))

 par(mfrow=c(1,2))
 z<-mcmc.b[,6]
 n<-length(z)
 plot(z, type='l', main='Coefficient of X_5')
 z<-mcmc.s
 n<-length(z)
 plot(z, type='l', main='Sigma')

 par(mfrow=c(1,2))
 z<-mcmc.b[,6]
 n<-length(z)
 plot(z, type='l', main='Coefficient of X_5')
 lines(1:n,cumsum(z)/1:n, col=2)
 z<-mcmc.s
 n<-length(z)
 plot(z, type='l', main='Sigma')
 lines(1:n,cumsum(z)/1:n, col=2)


par(mfrow=c(1,2))
 z<-mcmc.b[,6]
 n<-length(z)
 plot(1:n,cumsum(z)/1:n, , type='l', main='Coefficient of X_5')
 z<-mcmc.s
 n<-length(z)
 plot(1:n,cumsum(z)/1:n, type='l', main='Sigma')




#------------------------------------------------
library(MCMCpack)
mcmcpack.model<-MCMCregress(y~., data=simex1, b0=0, B0=1/1000, c0=0.001, d0=0.001)
head(mcmcpack.model)
dim(mcmcpack.model)

summary(mcmcpack.model)
plot(mcmcpack.model)


par(mfrow=c(1,2))
 z<-as.vector(mcmcpack.model[,6])
 n<-length(z)
 plot(z, type='l', main='Coefficient of X_5')
 lines(1:n,cumsum(z)/1:n, col=2)
 z<- as.vector(mcmcpack.model[,17])
 n<-length(z)
 plot(z, type='l', main='Sigma')
 lines(1:n,cumsum(z)/1:n, col=2)


par(mfrow=c(1,2))
 z<-as.vector(mcmcpack.model[,6])
 n<-length(z)
 plot(1:n,cumsum(z)/1:n, , type='l', main='Coefficient of X_5')
 z<-as.vector(mcmcpack.model[,17])
 n<-length(z)
 plot(1:n,cumsum(z)/1:n, type='l', main='Sigma')
#------------------------------------------------
library(bayesm)
simex1.list <- list(y=simex1[,1], X=as.matrix(cbind(1,simex1[,-1])) )
p<-ncol(simex1.list$X)
prior1 <- list( betabar=rep(0,p), A=0.01*diag(p), nu=0.02, ssq=1) 
mcmc1  <- list( R=1100, keep=1 )
bayesm.model <- runireg(Data=simex1.list, Prior=prior1, Mcmc=mcmc1)

bayesm.model <- runireg(Data=simex1.list, Mcmc=list(R=1100))

names(bayesm.model)
head(bayesm.model$betadraw)
dim(bayesm.model$betadraw)
head(bayesm.model$sigmasq)
dim(bayesm.model$sigmasq)

bayesm.model <- runireg(Data=simex1.list, Prior=prior1, Mcmc=mcmc1)
summary(bayesm.model$betadraw, burnin=100)
summary(bayesm.model$sigmasq, burnin=100)

z1<- apply(bayesm.model$betadraw, 2, mean)
z2<- apply(bayesm.model$betadraw, 2, sd)
z3<- c( mean(bayesm.model$sigmasqdraw), sd(bayesm.model$sigmasqdraw) )
z <- rbind( cbind(z1,z2), z3)
rownames(z) <- c('Const', paste('X', 1:15, sep=''), 'sigmasq') 
colnames(z) <- c('Mean', 'SD') 
round(z,3)

 

summary(bayesm.model$betadraw, burnin=100)
summary(bayesm.model$sigmasq, burnin=100)
 
 
plot(bayesm.model$betadraw)

plot.bayesm.mat(bayesm.model$betadraw[,1:3])
plot(bayesm.model$sigmasq) 



#------------------------------------------------


library(bayesm)
simex1.list <- list(y=simex1[,1], X=as.matrix(cbind(1,simex1[,-1])) )
p<-ncol(simex1.list$X)
prior1 <- list( betabar=rep(0,p), A=0.01*diag(p), nu=0.02, ssq=1) 
mcmc1  <- list( R=5100, keep=1 )
bayesm.model2 <- runiregGibbs(Data=simex1.list, Prior=prior1, Mcmc=mcmc1)
              
  
names(bayesm.model2)
head(bayesm.model2$betadraw)
dim(bayesm.model2$betadraw)
head(bayesm.model2$sigmasq)
dim(bayesm.model2$sigmasq)
         
plot(bayesm.model$betadraw)
 plot.bayesm.mat(bayesm.model$betadraw[,1:3])
 plot(bayesm.model$sigmasq) 
 
             
#------------------------------------------------
#G-prior
#------------------------------------------------
y<-simex1$y
x<-cbind( 1, as.matrix(simex1[,-1] ) )
n<-length(y)
p<-ncol(x)

g<-n
S <- g*solve( t(x) %*% x )
res<-lm.bayes( y, x, prior.mean=rep(0,p), V=S, add.constant=FALSE )

#----------------------------------------------------
library(bayesm)
y<-simex1$y
x<-cbind( 1, as.matrix(simex1[,-1] ) )
n<-length(y)
p<-ncol(x)
g<-n
S <- g*solve( t(x) %*% x )
simex1.list <- list(y=y, X=x )
prior1 <- list( betabar=rep(0,p), A=S, nu=0.02, ssq=1) 
mcmc1  <- list( R=5000, keep=1 )
bayesm.model <- runireg(Data=simex1.list, Prior=prior1, Mcmc=mcmc1)
summary(bayesm.model$betadraw)

round(cbind(res$beta, apply( bayesm.model$betadraw,2,mean ) 
              , apply( bayesm.model$betadraw,2,sd )), 3)


summary(bayesm.model$sigmasqdraw)
res$s2

#------------------------------------------------
# Jeffreys prior
#------------------------------------------------


y<-simex1$y
x<-cbind( 1, as.matrix(simex1[,-1] ) )
n<-length(y)
p<-ncol(x)

g<- 1000*n
S <- g*solve( t(x) %*% x )
res.j<-lm.bayes( y, x, prior.mean=rep(0,p), V=S, add.constant=FALSE, a=(-16/2+1), b=0.000001)

res$s2
res.j$s2
MLE <- lm(y~., data=simex1)
summary(MLE)$sigma^2

#------------------------------------------------
# Jeffreys prior
#------------------------------------------------

library(bayesm)
y<-simex1$y
x<-cbind( 1, as.matrix(simex1[,-1] ) )
n<-length(y)
p<-ncol(x)
g<-n
S <- solve( t(x) %*% x )/g
simex1.list <- list(y=y, X=x )
prior1 <- list( betabar=rep(0,p), A=S, nu=-14, ssq=0) 
mcmc1  <- list( R=5000, keep=1 )
bayesm.model <- runireg(Data=simex1.list, Prior=prior1, Mcmc=mcmc1)
summary(bayesm.model$betadraw)

round(cbind(res.j$beta, apply( bayesm.model$betadraw,2,mean ) 
              , apply( bayesm.model$betadraw,2,sd )), 3)


summary(bayesm.model$sigmasqdraw)
res.j$s2

