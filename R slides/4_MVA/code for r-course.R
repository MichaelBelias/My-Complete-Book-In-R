# mentioned in this tutorial
RatAUEB <- c("HDclassif", "corrgram", "cluster", "mclust","FactMixtAnalysis","nnet","class","tree")
# Install the packages
install.packages(RatAUEB, repos = "http://cran.rstudio.com/")

#packages installed
# we need now to call them

library(corrgram)
library(HDclassif)
library(cluster)
library(mclust)
library(FactMixtAnalysis)
library(nnet)
library(class)
library(tree)

###### START



#### distances in R
dist1<-dist(winedata[,2:4])
dist2<-daisy(winedata[,2:4])
x0<-winedata[,-1]
dM2 <- as.dist(apply(x0, 1, function(i) mahalanobis(x0, i, cov = cov(x0))))


corrgram(winedata[,-1])
pairs(winedata[,-1])
pairs(winedata[,-1], col=winedata[,1])
hclust(dist(winedata[,-1]),method="complete")
hc1<-hclust(dist(winedata[,-1]),method="complete")
summary(hc1)
plot(hc1$height)



### run with other linkage
hc2<-hclust(dist(winedata[,-1]),method="ward.D")
hc3<-hclust(dist(winedata[,-1]),method="single")
hc4<-hclust(dist(winedata[,-1]),method="average")



#### create clasifications
clas1<-cutree(hc1, k=2:5)
clas2<-cutree(hc2, k=2:5)
clas3<-cutree(hc3, k=2:5)
clas4<-cutree(hc4, k=2:5)


##### random noise data  ########

fakedata<-mvrnorm(100,mu=rep(0,5), Sigma=diag(1,5))
fake<- hclust(dist(fakedata),method="ward.D")
plot(fake)
plot(fake$height)



####replicate fake data
fakedata<-mvrnorm(100,mu=rep(0,5), Sigma=diag(1,5))
fake<- hclust(dist(fakedata),method="ward.D")
points(fake$height)

##### data with two clusters
fakedata<-mvrnorm(100,mu=rep(0,5), Sigma=diag(1,5))
fakedata<-rbind(fakedata,mvrnorm(100,mu=rep(5,5), Sigma=diag(1,5)))
fake<- hclust(dist(fakedata),method="ward.D")
plot(fake)
plot(fake$height)


#####################################################

plot(silhouette(clas1[,1], dist(winedata[,-1])))


wine<-as.matrix(winedata[,-1])
m <- manova(wine~clas1[,2])
summary(m,test="Wilks")

mywilks<- summary(m,test="Wilks")$stats[1,2]


######################################################

adjustedRandIndex(clas1[,2],clas2[,2])



indices<-NULL
for  ( i in 3:28) {
usevar<- winedata[,2:i]
myclust<-hclust(dist(usevar),method="complete")
myclass<- cutree(myclust,3)
indices<-c(indices,adjustedRandIndex(myclass,winedata[,1]))
}


#############################################################
### K-MEANS 
#############################################################

km1<- kmeans(wine, 3)
km1$cluster	
km1$centers	
km1$totss	
km1$withinss	
km1$tot.withinss	
km1$betweenss	
km1$size	
km1$iter	
km1$ifault	

table(km1$cluster,clas1[,2])


km2<- kmeans(wine, 3, algorithm="Lloyd")
 table(km1$cluster, km2$cluster)


#################################################################
###  mclust
#################################################################

wine<-winedata[,2:7]
mc1<-Mclust(wine, G=2:5, modelNames=c("EII", "VII", "EEI", "EVI", "VEI", "VVI"))


mc1$G	#The optimal number of mixture components.
mc1$BIC	#All BIC values.
mc1$bic	#Optimal BIC value.
mc1$loglik	#The loglikelihood corresponding to the optimal BIC.
mc1$df	 #The number of estimated parameters.
mc1$parameters	#A list with the following components:
mc1$pro #A vector whose kth component is the mixing proportion for the kth component 
mc1$mean #The mean for each component. 
mc1$variance  #A list of variance parameters for the model. 
mc1$z	 #posterior probabilities.
mc1$classification	 #map(z): The classification corresponding to z.
mc1$uncertainty	#The uncertainty associated with the classification.

plot(mc1)
summary(mc1)


#### FActore analyzers
mymfa<-fma(wine, 3, 2)

mymfa$H   #The estimated factor loading matrix.
mymfa$lik  #The log-likelihood computed at each iteration of the EM algorithm.
mymfa$w   #A matrix with the estimated weights of the mixture.
mymfa$Beta  # estimated component means of the mixture.
mymfa$phi #  coefficients of the covariates
mymfa$sigma  #An array of dimension k x r x r which contains the estimated component covariance
mymfa$psi  #The noise diagonal variance matrix.
mymfa$index  #The allocation vector.
mymfa$bic   #The BIC value.
mymfa$aic  #The AIC value.
mymfa$elapsed #Computational time in seconds.

#####  how many parameters we need to estimate?
mymfaall<-fma(winedata[,2:14], 4, 2)





###############################################################


ncol <- 304
nrow <- 268

##pict <- matrix(scan("c:\\paulia.rgb"), ncol=3, byrow=T)
row<-rep(1:nrow,each=ncol)
column<-rep(1:ncol,nrow)
pixels <- which(apply(pict, 1, function(x) any(x)>0 & all(x<255)))
npixels <- length(pixels)

metr<-0
plot(c(0,nrow),c(0,ncol),type="n",xlab="",ylab="")
for (i in 1:nrow) {
  for (j in 1:ncol){
metr<-metr+1
points(i,j,col=rgb(pict[metr,1]/255,pict[metr,2]/255,pict[metr,3]/255))
}}

plot(c(0,nrow),c(0,ncol),type="n",xlab="",ylab="",axes=FALSE)
points(row[pixels],column[pixels],col=rgb(pict[pixels,1]/255,
       pict[pixels,2]/255,pict[pixels,3]/255))

mydata<-cbind(pict[pixels,],row[pixels],column[pixels])

### Change the next lines as appropriate
setsize <- 2000
trset <- sample(npixels, setsize, replace=FALSE)
minclus <- 3
maxclus <- 15
Modelnames<- c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV")


### Go!
X <- pict[pixels[trset], ]
bicvals <- mclustBIC(X, G=minclus:maxclus, modelNames = Modelnames)
sumry <- summary(bicvals, X)
znew <- do.call("estep", c(list(data=pict[pixels,]), sumry))
classif <- map(znew$z)

newcolors<-cbind(tapply(pict[pixels,1],classif,mean),
tapply(pict[pixels,2],classif,mean),
tapply(pict[pixels,3],classif,mean) )

plot(c(0,nrow),c(0,ncol),type="n",xlab="",ylab="",axes=FALSE)
points(row[pixels],column[pixels],col=rgb(newcolors[classif,1]/255,
       newcolors[classif,2]/255,newcolors[classif,3]/255))


##################################################################################



wines<-as.data.frame(winedata[,c(1,2,15,28)])




########### STEP 2: LDA
m1<-lda(Type~., data=wines)
m2<-predict(m1)
table(wines[,1],m2$class)
m3<-lda(Type~., data=wines, CV=TRUE)
table(wines[,1],m3$class)

par(mfrow=c(3,1))
plot(m3$posterior[,1], col=winedata[,1])
plot(m3$posterior[,2], col=winedata[,1])
plot(m3$posterior[,3], col=winedata[,1])

plot(m2$x,col=winedata[,1])
 plot(m1)

########### QDA
mq1<-qda(Type~., data=wines)
mq2<-predict(mq1)
table(wines$Type,mq2$class)


###########################  multinomial logistic
mult <- multinom(Type ~ ., wines)
summary(mult)
plot(mult$fitted,m2$posterior)
mult.class<- apply(mult$fitted,1,which.max)
table(mult.class,m2$class)



############################  knn
library(class)
km1<-knn(winedata[,2:7],winedata[,2:7], cl=winedata[,1],k=4)
table(winedata[,1],km1)

km2<-knn(winedata[,2:7],winedata[,2:7], cl=winedata[,1],k=7)
table(winedata[,1],km2)

knn.cv(winedata[,2:7], cl=winedata[,1],k=7)



########################  tree
wines<-as.data.frame(winedata[,1:15])



fit1<-tree(as.factor(Type)~.,data=wines)
fit1$frame
fit1$where
table(wines$Type, predict(fit1,type='class'))
plot(fit1)
summary(fit1)
 text(fit1)

fit2<-tree(as.factor(Type)~.,data=wines, split="gini")
plot(fit2)
text(fit2)


predict(fit1,type='vector')
 predict(fit1,type='class')
 predict(fit1,type='tree')



#####################################################################
###  variable selection using the diabetes data
###  the variables are
###
### Number of times pregnant 
### Plasma glucose concentration a 2 hours in an oral glucose tolerance test 
### Diastolic blood pressure (mm Hg) 
### Triceps skin fold thickness (mm) 
### 2-Hour serum insulin (mu U/ml) 
### Body mass index (weight in kg/(height in m)^2) 
### Diabetes pedigree function 
### Age (years) 
### Class variable (1 or 2) 




re<-NULL


for (k in 1:20) {
print(k)
t<-NULL

for (i in 1:100) {

te<- sample(1:768,k)


train <- diabetes[-te,2:4]
test <-   diabetes[te,2:4]
cl <- factor(diabetes[-te,9])
z <-  lda(train, cl)
pr<-  predict(z, test)$class
t<- c(t,   sum(diabetes[te,9] == pr)    /dim(test)[1])
}


re<-c(re,mean(t))
}


###########################################################
###
###   k-fold cross validation
###
############################################################


deiktes<-sample(1:768)
variab<-2:8


re<-NULL


for (k in c(1,2,3,4,6,8,12,16,24)) {
print(k)

omades<- 768/k

t<-NULL

for (i in 1:omades) {
te<- deiktes[ ((i-1)*k+1):(i*k)]
train <- diabetes[-te,variab]
test <-   diabetes[te,variab]
cl <- factor(diabetes[-te,9])
z <-  qda(train, cl)
pr<-  predict(z, test)$class
t<- c(t,   sum(diabetes[te,9] == pr)    /dim(test)[1])
}


re<-c(re,mean(t))
}


###########################################################
###
###  variable selction using  WILKS LAMBDA  
###  (we can use minor changes to use accuracy 
###
############################################################


allcomb<-combn(2:8,7)
tw<-NULL
for (i in 1:dim(allcomb)[2]) {
model<-manova(as.matrix(diabetes[,allcomb[,i]])~diabetes[,9])
tw<-c(tw,summary(model,test="Wilks")$stats[1,2])
}
cbind(t(allcomb),tw)
cbind(t(allcomb),tw)[which.min(tw),]

