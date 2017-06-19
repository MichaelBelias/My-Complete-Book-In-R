###Function for checking if pkg is installed 
is_installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 

####Packages needed
if(!is_installed("devtools"))  
{  
  install.packages("devtools", dependencies = T)  
}  
library("devtools",character.only=TRUE,quietly=TRUE,verbose=FALSE) 


package_names= c("hglm","foreign","knitr","caret","rstan","rstanarm","ggplot2","ggthemes","ggThemeAssist","GGally","reshape2",
                 "lme4","nlme","compiler","metafor","sas7bdat","memisc","parallel",
                 "boot","nlme", "MASS","lme4","Rglpk","mvmeta","netmeta",
                 "gemtc", "msm","Formula","gems", "Epi", "simMSM","pander",
                 "gemtc","coda", "ipdmeta","R2jags","ape","knitr","ggmcmc",
                 "RCurl","XML","ips","RMySQL","svDialogs","Rphylip",
                 "xlsx","devtools","httpuv","png","grid","phytools","phangorn",
                 "car","survival","corrplot","COUNT",
                 "knitr","diagram","stringi","metafor","rdrop2")





###Package installer and installation
for(package_name in package_names)  
{  
  if(!is_installed(package_name))  
  {  
    install.packages(package_name, dependencies = T)  
  }  
  library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE) 
  cat("Package :",package_name," Acquired","\n")
  
}  

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()

rm(package_name,package_names,is_installed,detachAllPackages)
search()
