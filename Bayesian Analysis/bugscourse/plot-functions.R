# A few convenient plotting functions for CODA objects.
#
# These functions are not designed to be entirely
# genera\l, so they might not work on examples more
# complex than those in the course - but they be a
# useful starting point.


# CODA's densplot() function is misleading for discrete parameters
# densplot.discrete() is a simple alternative
densplot.discrete <- function(x,  ...){
  xx <- as.matrix(x)
  for (i in 1:nvar(x)){
    y <- xx[,  i,  drop = TRUE]
    barplot(prop.table(table(y)),  yaxs = "r",  ...)
    }
  return(invisible(x))
}

# The following functions depend on the R package "denstrip":
# install.packages("denstrip") # if necessary
library(denstrip)

## Get CODA samples of a vector-valued variable "varname".
## Rows are MCMC iterations, with multiple chains stacked
## on top of each other
## Columns are indices of the vector
## The functions below depend on this
get.coda.matrix <- function(coda, varname){
    var.cols <- grep(paste("^",varname,"\\[",sep=""), varnames(coda))
    if (length(var.cols) > 0){
        m <- do.call("rbind", coda[,var.cols])
        indnames <- gsub(paste0(varname, "\\[([0-9]+)\\]"), "\\1", colnames(m))
        m[, order(as.numeric(indnames))]
    } else {
        stop("No variable '", varname, "' found in coda object")
    }
}

## One possible R implementation of Inference->Compare->box plots
## in Win/OpenBUGS
bugs.boxplot <- function(coda, varname,ordered=FALSE, ...){
    res <- get.coda.matrix(coda, varname)
    qs <- apply(res, 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
    if (ordered) qs <- qs[,order(qs[3,])]
    plot(0, type="n", ylim=range(qs) * c(0.96, 1.04), xlim=c(0,ncol(res)), ...)
    cistrip(t(qs[c(3,2,4),]), horiz=FALSE, at=1:ncol(res),
            col="blue", d=0, cex=0, lwd=5)
    cistrip(t(qs[c(3,1,5),]), horiz=FALSE, at=1:ncol(res))
}

## R implementation of density strips
## (Inference->Compare->density strips in OpenBUGS --
## not in WinBUGS)
## A smooth alternative to boxplots that illustrates
## posterior density by shading darkness.
## Have a look at help(denstrip) for options to
## tweak the appearance
density.strips <- function(coda, varname, ...){
    res <- get.coda.matrix(coda, varname)
    qs <- apply(res, 2, quantile, c(0.025, 0.5, 0.975))
    plot(0, type="n", ylim=c(-4, 4), xlim=c(0,ncol(res)), ...)
    for (i in 1:ncol(res)){
        denstrip(res[,i], at=i, horiz=FALSE, width=0.5, colmax="blue",
        ticks=qs[c(1,3),i], mticks=qs[2,i], tcol="black", mcol="black")
    }
}


## R implementation of Inference->Compare->caterpillar plot
## in Win/OpenBUGS
catplot <- function(coda, varname, ordered=FALSE){
    res <- get.coda.matrix(coda, varname)
    qs <- apply(res, 2, quantile, c(0.025, 0.5, 0.975))
    if (ordered) qs <- qs[,order(qs[2,])]
    plot(0, type="n", yaxt="n", ylab="", xlab=varname, xlim=c(0, 1),
         ylim=c(0,ncol(res)))
    abline(v=mean(qs[2,]), col="red")
    cistrip(t(qs[c(2,1,3),]), horiz=TRUE, at=ncol(res):1, d=0)
    indnames <- gsub(varname, "", colnames(qs))
    text(x=qs[3,], y=ncol(res):1, labels=indnames, col="blue",
         cex=0.5, pos=4)
}

## R implementation of the Win/OpenBUGS Inference->Compare->Model Fit
## function for plotting posterior fitted values from a regression model
model.fit <- function(y, # Data for outcome in regression (numeric)
                      x, # Data for predictor in regression (numeric)
                      coda, # CODA R object from BUGS/JAGS model fit
                      mu, # Variable name (character) for the
                          # fitted values (E(Y)) in the
                          # regression model
                      quantiles=c(0.025, 0.5, 0.975), # Quantiles of
                                                      # the fitted
                                                      # values to plot
                      ... # Other arguments to control the
                          # appearance of the plot, passed
                          # to the R functions plot()
                          # (plot.default()) and lines().
                          # See example below
                      ){
    res <- get.coda.matrix(coda, mu)
    qs <- apply(res, 2, quantile, quantiles)
    plot(x, y, ...)
    lines(x, qs[1,], col="lightblue", ...)
    lines(x, qs[2,], col="red", ...)
    lines(x, qs[3,], col="lightblue", ...)
}
