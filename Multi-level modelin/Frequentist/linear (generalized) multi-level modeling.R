# Or to install the dev version
#library(devtools)
#install_github("lme4", user = "lme4")


library(lme4)  # load library
library(arm)  # convenience functions for regression in R
lmm.data <- read.table("http://www.unt.edu/rss/class/Jon/R_SC/Module9/lmm.data.txt", 
                       header = TRUE, sep = ",", na.strings = "NA", dec = ".", strip.white = TRUE)
# summary(lmm.data)
head(lmm.data)

