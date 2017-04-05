##########################
### Lecture 3 - R code ###
##########################

# Leukemia Data
leukem <- read.csv("C:/Users/christofer/Dropbox/Applied_Survival_Analysis_Jan2016/Data/leukem.csv")
head(leukem)

# KM estimtes for both arms
library(survival)
leukem.fit = survfit( Surv(weeks,remiss) ~ trt, data = leukem)
summary(leukem.fit)

# Plot of KM estimates
plot(leukem.fit, mark.time = F,main = "Comparison of Treatments for Leukemia",
     ylab = "Survival Probability",
     xlab = "Time from Remission to Relapse (weeks)",lty = 1:2,col = c("black","red"))
legend("topright",lty = 1:2,col = c("black","red"),
       legend = c("Control (N=21)","6-MP (N=21)"),bty = "n")

# Logrank - the default
survdiff(Surv(weeks,remiss) ~ trt, data = leukem)

# Peto & Peto modification of the Gehan-Wilcoxon test (not the same with STATA)
survdiff(Surv(weeks,remiss) ~ trt, data = leukem, rho = 1)


####################
# P sample logrank #
####################
# Time taken to finish a test with 3 different noise distractions. All
# tests were stopped after 12 minutes.
data = data.frame(time = c(9.0,10.0,12.0,
                           9.5,12.0,12,
                           9.0,12,12,
                           8.5,11.0,12,
                           10.0,12.0,12,
                           10.5,10.5,12),delta = c(rep(1,3),c(1,1,0),c(1,0,0),rep(c(1,1,0),3)),
                  group = rep(1:3,6))

# Logrank
survdiff(Surv(time,delta)~group,data)

# Wilcxon
survdiff(Surv(time,delta)~group,data,rho = 1)


######################
# Stratified logrank #
######################
nurshome <- read.csv("C:/Users/christofer/Dropbox/Applied_Survival_Analysis_Jan2016/Data/nurshome.csv")
nurshome$age1 = 1*(nurshome$age>85)

survdiff( Surv(los,fail) ~ age1 + strata(gender),data = nurshome)


# Several tables required for comptution logrank test in leukemia data
log.rank = function(surv,fail,group)
{
  # Sort data with respect to time
  data = data.frame(time = surv,fail = fail,group = group)
  data = data[order(data$time),]
  
  time = data$time
  delta = data$fail
  trt = data$group
  
  # Distinct event times
  times.uni = unique(time[delta==1])
  K = length(times.uni)
  
  # Number of failures for each event time
  dj = table(time[delta == 1])
  
  dj.group0 = table(time[delta == 1 & trt == 0])
  dj.group1 = table(time[delta == 1 & trt == 1])
  
  dj.grp0 = rep(0,K)
  dj.grp0[match(as.numeric(names(dj.group0)),times.uni)] = dj.group0
  
  dj.grp1 = dj - dj.grp0
  
  # Number at risk for each event time
  rj = rep(NA,K)
  rj.grp0 = rep(NA,K)
  rj.grp1 = rep(NA,K)
  
  for (i in 1:K)
  {
    rj[i] = sum(time >= times.uni[i])
    rj.grp0[i] = sum(time[trt==0] >= times.uni[i])
  }
  
  rj.grp1 = rj - rj.grp0
  
  out = list()
  
  for (i in 1:K)
  {
    mat = matrix(nrow = 2,ncol = 2)
    colnames(mat) = c("Event","No event")
    rownames(mat) = c("Group0","Group1")
    
    mat[1,1] = dj.grp0[i]
    mat[2,1] = dj.grp1[i]
    mat[1,2] = rj.grp0[i] - dj.grp0[i]
    mat[2,2] = rj.grp1[i] - dj.grp1[i]
    
    out[[i]] = mat
  }
  names(out) = paste("RiskSet",times.uni)
  
  out$chi2 = sum(dj.grp0 - rj.grp0*dj/rj)^2/sum( rj.grp0*rj.grp1*dj*(rj-dj)/(rj^2*(rj-1)))
  
  return(out)
}

log.rank(leukem$weeks,leukem$remiss,leukem$trt)
