library(KMsurv)
library(nlme)
angina <- read.csv("h:/teaching/Athens/BIO223 (Survival-Yiannoutsos)/data/angina.csv")

died<- gsummary(angina[angina$status==1,], sum, groups=angina[angina$status==1,]$years)
censored<- gsummary(angina[angina$status==0,], sum, groups=angina[angina$status==0,]$years)
# Time length must be one plus the length of everyone else's
lt<-length(died[,1])
years<-died[,1]
years[lt+1]<-NA

LT.angina<-data.frame(cbind(died[,1], died[,3], censored[,3]))
names(LT.angina)<-c("years", "died", "censored")

#Life table for the angina example
lifetab(years, sum(LT.angina$died+LT.angina$censored),
        LT.angina$censored, LT.angina$died)

############################# Life table when the data are not grouped ###########################
nurshome <- read.csv("h:/teaching/Athens/BIO223 (Survival-Yiannoutsos)/data/nurshome.csv")
nurshome<-nurshome[nurshome$rx==1,]
nurshome[1:5,]

los100<-floor(nurshome$los/100)
died<- gsummary(as.data.frame(nurshome$fail), sum, groups=los100)
total<- gsummary(as.data.frame(nurshome$fail), length, groups=los100)
censored<-total-died

LT.treated<-data.frame(cbind(unique(los100), died, censored))
names(LT.treated)<-c("LOS100", "died", "censored")

#los100 must have one more length than everyone else
los100<-sort(LT.treated$LOS100)
lt<-length(los100)
los100[lt+1]<-NA

#Life table for the angina example
nursLT<-lifetab(los100, sum(total), LT.treated$censored, LT.treated$died)

summary(nursLT)

#Plot of the survival
plot(los100[0:lt],nursLT[,5], type="s", 
     ylab="Survival", xlab="LOS (100-day intervals)")

#Plot of the hazard
plot(los100[0:lt],nursLT[,7], type="s", 
     ylab="Hazard", xlab="LOS (100-day intervals)")
