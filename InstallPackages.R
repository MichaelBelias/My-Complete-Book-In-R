###Function for checking if pkg is installed 
is_installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 

####Packages needed
if(!is_installed("devtools"))  
{  
  install.packages("devtools", dependencies = T)  
}  
library("devtools",character.only=TRUE,quietly=TRUE,verbose=FALSE) 


package_names= c("hglm","foreign","knitr","ggplot2","ggthemes","ggThemeAssist","GGally","reshape2",
                 "lme4","nlme","compiler","metafor","sas7bdat","memisc","parallel",
                 "boot","nlme", "MASS","lme4","Rglpk","mvmeta","netmeta",
                 "gemtc", "msm","Formula","gems", "Epi", "simMSM","pander",
                 "gemtc","coda", "ipdmeta","R2jags","ape","knitr","ggmcmc",
                 "RCurl","XML","ips","RMySQL","svDialogs","Rphylip","rdrop2",
                 "xlsx","devtools","httpuv","png","grid","phytools","phangorn",
                 "rstan","rstanarm","car","survival","corrplot","COUNT",
                 "knitr","diagram","stringi","metafor","caret")





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

rm(package_name,package_names,is_installed)
search()
