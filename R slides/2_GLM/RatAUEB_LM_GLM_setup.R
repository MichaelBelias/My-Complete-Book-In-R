## This script installs the packages that will be used or will be
## mentioned in the module "Linear and Generalized Linear Models" of
## the R Summer School at AUEB.

## Author: Ioannis Kosmidis
## Date: 23/06/2014


# Please run the chunk PACKAGES below before attending the course

## ----PACKAGES--------------------------------------------
# Update the installed packages first
update.packages(ask = FALSE, repos = 'http://cran.rstudio.com/')
# Install the packages for this tutorial
RatAUEBpackages <- c('alr3', 'brglm', 'car', 'caret', 'elasticnet', 'fortunes', 'ggplot2', 'gridExtra', 'GGally', 'leaps', 'lmtest', 'mfp', 'relimp')
install.packages(RatAUEBpackages, repos = 'http://cran.rstudio.com/')
