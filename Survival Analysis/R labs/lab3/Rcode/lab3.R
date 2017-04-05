#######################################################
### Lab 3: Comparing survival curves between groups ###
#######################################################

# (a): Artificial data from 2 groups
surv = c(15,18,19,19,20,16,18,20,23,24)
fail = c(rep(1,5),c(0,1,0,1,0))
group = c(rep(0,5),rep(1,5))

data = data.frame(time = surv, fail = fail, group = group)

# Logrank test
library(survival)
survdiff( Surv(time,fail) ~ group, data = data)

# Compare the median survival between the two groups
quantile( survfit(Surv(time,fail) ~ group,data = data), probs = 0.5,conf.int = F )

# OPTIONAL: function to compute the logrank test
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

log.rank(data$time,data$fail,data$group)

# (b): Leukemia Data
leukem <- read.csv("C:/Applied_Survival_Analysis_Jan2016/lab3/data/leukem.csv")
leukem[1:4,]

# Encode trt as a factor
leukem$trt[leukem$trt == 1] = "6-MP"
leukem$trt[leukem$trt == 0] = "Control"

# Specify the order of the levels
# so that the control group will be the reference category in regression models
leukem$trt = factor(leukem$trt,levels = c("Control","6-MP"))

# (b)-(i)
# KM estimates in both groups
leukem.fit = survfit( Surv(weeks,remiss) ~ trt, data = leukem)
summary(leukem.fit)

# (b)-(ii)
# Plotting the KM estimates
setwd("C:/Applied_Survival_Analysis_Jan2016/lab3/graphs")

pdf("leukemKM.pdf",height = 6,width = 6)
plot(leukem.fit, mark.time = F,main = "Comparison of Treatments for Leukemia",
     ylab = "Survival Probability",
     xlab = "Time from Remission to Relapse (weeks)",lty = 1:2,col = c("black","red"))
legend("topright",lty = 1:2,col = c("black","red"),
       legend = c("Control (N=21)","6-MP (N=21)"),bty = "n")
dev.off()

# (b)-(iii)
# Logrank and Wilcoxon tests

# Logrank - the default
survdiff(Surv(weeks,remiss) ~ trt, data = leukem)

# Peto & Peto modification of the Gehan-Wilcoxon test
survdiff(Surv(weeks,remiss) ~ trt, data = leukem, rho = 1)

# (b)-(iv)
# Plot of estimated log cumulative hazard
pdf("leukemlogHaz.pdf",height = 6,width = 6)
plot(leukem.fit, mark.time = F,fun="cloglog",
     ylab = "Log cumulative hazard",
     xlab = "Time from Remission to Relapse (weeks)",lty = 1:2,col = c("black","red"))
legend("topleft",lty = 1:2,col = c("black","red"),
       legend = c("Control (N=21)","6-MP (N=21)"),bty = "n")
dev.off()


