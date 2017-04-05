########################################
### LAB 9: Time-dependent covariates ###
########################################

# (a)-(i): Read the famous Stanford heart transplant data set in R
stanford = read.csv("C:/Applied_Survival_Analysis_Jan2016/lab9/data/stanford.csv")
head(stanford)

# Fit a naive Cox model
library(survival)
naiveCox = coxph(Surv(time,fail) ~ transplant,data = stanford)
summary(naiveCox)

# (a)-(ii): Create time-updated transplantation status
# We need to trasnform the dataset

# First, let's have a look at the original data
# Please, note that NA means missing in R
stanford[stanford$patid %in% c(12,16,38,80),
         c("patid","transplant","waitime","time","fail")]

# We want to replace each observation in the dataset with 2 copies of the observation
# only for the patients who got a new heart
n = nrow(stanford)
stanford2 = stanford[rep(1:n,each = 2),]

# Indexing observations by patient
stanford2$ord = rep(1:2,n)

# Keep 1 observation for patients without a transplantation
stanford2 = stanford2[-which(stanford2$tr == 0 & stanford2$ord == 2),]

# see
stanford2[stanford2$patid %in% c(12,16,38,80),
         c("patid","transplant","waitime","time","fail","ord")]

# First row of patients with transplant == 1
sub = stanford2$ord == 1 & stanford2$tr == 1

# Second row of patients with trasplant == 1
sub2 = stanford2$ord == 2 & stanford2$tr == 1

##############################################
# Now we need to create entry and exit times #
##############################################

# Create exit times and modify the transplantation and failure history
stanford2$tstop = stanford2$time

stanford2$tstop[sub] = stanford2$waitime[sub]
stanford2$fail[sub] = 0
stanford2$transplant[sub] = 0

stanford2$tstart = 0
stanford2$tstart[sub2] = stanford2$tstop[sub]

# see
stanford2[stanford2$patid %in% c(12,16,38,80),]
Surv(stanford2$tstart,stanford2$tstop,stanford2$fail)[stanford2$patid 
                                                      %in% c(12,16,38,80)]

# We fit a model involving the time-updated transplantation status (see ?Surv)
fitCox = coxph( Surv(tstart,tstop,fail) ~ transplant,data = stanford2)
summary(fitCox)

##############################
### Time dependent effects ###
##############################

# Leukemia Data
leukem <- read.csv("C:/Applied_Survival_Analysis_Jan2016/lab3/data/leukem.csv")
leukem[1:4,]

# Encode trt as a factor
leukem$trt[leukem$trt == 1] = "6-MP"
leukem$trt[leukem$trt == 0] = "Control"

# Specify the order of the levels
# so that the control group will be the reference category in regression models
leukem$trt = factor(leukem$trt,levels = c("Control","6-MP"))

# KM estimates in both groups
leukem.fit = survfit( Surv(weeks,remiss) ~ trt, data = leukem)

setwd("C:/Applied_Survival_Analysis_Jan2016/lab9/graphs")

# Plot of estimated log cumulative hazard
pdf("leukemlogHaz.pdf",height = 6,width = 6)
plot(leukem.fit, mark.time = F,fun="cloglog",
     ylab = "Log cumulative hazard",
     xlab = "Time from Remission to Relapse (weeks)",lty = 1:2,col = c("black","red"))
legend("topleft",lty = 1:2,col = c("black","red"),
       legend = c("Control (N=21)","6-MP (N=21)"),bty = "n")
dev.off()

# (b)-(ii): Incorrect model !!!
fit = coxph( Surv(weeks,remiss) ~ trt + I(weeks*(trt=="6-MP")),data = leukem)
summary(fit)

# You have to use the tt argumnent, 
# Decode trt as a numeric variable
leukem$trt = 1*(leukem$trt == "6-MP")
fit = coxph( Surv(weeks,remiss) ~ trt + tt(trt),data = leukem, 
             tt = function(x,t,...){x*t})
summary(fit)
