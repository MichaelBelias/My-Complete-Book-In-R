#############################################################
## FUNCTIONS FOR PROCESSING POSTERIORS AND MODEL CRITICISM ##
#############################################################


## Extract particular node from mcmc.list
getNode <- function(mcList, node)
{
  tmp <- lapply(mcList, function(x, node) { x[, grep(node, colnames(x))] }, node = node)
  if(is.null(dim(tmp[[1]])))
    do.call(c, tmp)
  else
    do.call(rbind, tmp)
}

## Check traces of posterior samples for convergence
## node is a string giving a regular expression that can be supplied to grep
chkTraces <- function(mcList, node)
{
  ## get posterior samples for specified parameter
  tmp <- as.mcmc.list(lapply(mcList, function(x, node) { mcmc(x[, grep(node, colnames(x))]) }, node = node))
  plot(tmp, density = T, smooth = F, ask = T)
}


## Get summary statistics when posterior samples stored in a matrix (rows = samples, columns = parameters)
sumStats <- function(post)
{
    mean <- apply(post, 2, mean)
    sd <- apply(post, 2, sd)
    medcrI <- apply(post, 2, quantile, probs = c(0.5,0.025,0.975))

    ret <- t(data.frame(mean, sd, t(medcrI)))
    rownames(ret) <- c("Mean","SD","Median","2.5%ile","97.5%ile")
    ret
}




####################################
## (Saturated) deviance functions ##
####################################

## Poisson
## y is observed, theta is mean
dev.pois <- function(theta, y){
    ret <- 2 * ((theta - y) - y*(log(theta)-log(y)))
    if(0 %in% y)
    {
        idx <- which(y == 0, arr.ind = T)
        ret[idx] <- 2 * (theta[idx] - y[idx])
    }
    ret
}

## Binomial
## y observed out of n, n*theta is mean
dev.bin <- function(theta, y, n){
    ret <- 2 * ((y * (log(y) - log(n*theta))) + ((n-y) * (log(n-y) - log(n - (n*theta)))))
    if(0 %in% y)
    {
        idx <- which(y == 0, arr.ind = T)
        ret[idx] <- 2 * (n[idx]-y)[idx] * (log(n[idx]-y[idx]) - log(n[idx] - (n[idx]*theta[idx])))
    }
    if(any(y == n))
    {
        idx <- which(y == n, arr.ind = T)
        ret[idx] <- 2 * y[idx] * (log(y[idx]) - log(n[idx]*theta[idx]))
    }
    ret
}

## Negative Binomial
## y failures observed before r successes observed, r*(1-theta)/theta is mean
dev.negbin <- function(theta, r, y){
    2 * (r*(log(r) - log(theta * (r+y))) + y*(log(y) - log((1 - theta)*(r+y))))
}

## Binomial when given eta = logit(p)
dev.bin.lgt <- function(eta, y, n){
    -2 * ((y * (eta - log(y/n) + log((n-y)/n))) - (n * (log(1 + exp(eta)) + log((n-y)/n))))
}


## Negative Binomial when given s = log(r), eta = logit(theta)
## y failures observed before r successes observed, r*(1-theta)/theta is mean
dev.negbin.lgt <- function(eta, s, y){
    2 * ((exp(s)*(s - eta)) - ((exp(s) + y) * (log(exp(s) + y) - log(1 + exp(eta)))) + (y * log(y)))
}


## Inverse logit function
expit <- function(x) { exp(x) / (1 + exp(x)) }



## Function to calculate deviance summaries, if likelihoods among Poisson, Binomial, NegBin
##
## post = posterior samples in mcmc.list format
## daty = list of numerators if proportion data, or alternatively counts
## datn = list of denominators if proportion data
## pNames = vector of parameter names needed for calculating deviances
## likTypes = vector of strings giving likelihoods assumed, possible values:
## "poisson", "binomial", "negbin", "poisson.log", "binomial.lgt", "negbin.lgt"
## dependent on whether want to calculate deviances on original or log/logit scale of parameters
## aux = vector of auxiliary parameter names if needed, i.e. overdispersion parameter for negative binomial
## devstem = string used in grep command to extract monitored deviances from mcmc.list
## yhatstem = string used in grep command to extract monitored predictions from mcmc.list
## plugin = string, either "mean" or "median", to choose whether posterior mean or median is
## used as the plug-in estimate in calculating Dhat, the deviance calculated at the plug-in estimate
##
## Outputs summary table of observed data, corresponding estimates, posterior mean deviance Dbar, plug-in deviance Dhat, effective number of parameters pD, and Deviance Information Criterion DIC
getDevs <- function(post, daty, datn, pNames, likTypes, aux = NULL, devstem, yhatstem, plugin)
    {
        require("gtools")
        
        ## Data
        yall <- unlist(lapply(daty, t))
        nall <- unlist(lapply(datn, t))
        
        ## MODEL CRITICISM
        dev <- getNode(post, devstem)
        yhat <- getNode(post, yhatstem)

        p <- list()
        for(m in 1:length(pNames))
        {
            p[[m]] <- getNode(post, pNames[m])
        }

        if(!is.null(aux))
        {
            auxParams <- list()
            for(a in 1:length(aux))
            {
                auxParams[[a]] <- getNode(post, aux[a])
            }
            auxParams <- do.call(cbind, auxParams)
            auxParams <- auxParams[,mixedorder(colnames(auxParams))]
        }

        ## sort columns
        dev <- dev[,mixedorder(colnames(dev))]
        yhat <- yhat[,mixedorder(colnames(yhat))]
        p <- do.call(cbind, p)
        ms <- mixedorder(colnames(p))
        p <- p[,ms]
        if(length(pNames) > 1)
        {
            yall <- yall[ms]
            nall <- nall[ms]
        }

        ## posterior mean deviances
        Dbar <- apply(dev, 2, mean)

        ## deviances calculated at posterior mean or median of parameters p
        if(plugin == "median")
        {
            yhatbar <- apply(yhat, 2, median)
            pbar <- apply(p, 2, median)
            if(!is.null(aux))
                auxbar <- apply(auxParams, 2, median)
        }
        else if(plugin == "mean")
        {
            yhatbar <- apply(yhat, 2, mean)
            pbar <- apply(p, 2, mean)
            if(!is.null(aux))
                auxbar <- apply(auxParams, 2, mean)
        }

        ## Dhat
        Dhat <- Dbar
        for(i in 1:length(likTypes))
        {
            idx <- grep(pNames[i], colnames(p))
            if(!is.null(aux))
                idxA <- grep(aux[i], colnames(auxParams))
            if(likTypes[i] == "binomial")
            {
                Dhat[idx] <- dev.bin(theta = pbar[idx], y = yall[idx], n = nall[idx])
            }
            else if(likTypes[i] == "poisson")
            {
                Dhat[idx] <- dev.pois(theta = pbar[idx], y = yall[idx])
            }
            else if(likTypes[i] == "negbin")
            {
                Dhat[idx] <- dev.negbin(theta = pbar[idx], r = auxbar[idxA], y = yall[idx])
            }
            else if(likTypes[i] == "binomial.lgt")
            {
                Dhat[idx] <- dev.bin.lgt(eta = pbar[idx], y = yall[idx], n = nall[idx])
                pbar[idx] <- apply(expit(p[,idx]), 2, mean)
            }
            else if(likTypes[i] == "poisson.log")
            {
                Dhat[idx] <- dev.pois.log(eta = pbar[idx], y = yall[idx])
                pbar[idx] <- apply(exp(p[,idx]), 2, mean)
            }
            else if(likTypes[i] == "negbin.lgt")
            {
                Dhat[idx] <- dev.negbin.lgt(eta = pbar[idx], s = auxbar[idxA], y = yall[idx])
                pbar[idx] <- apply(expit(p[,idx]), 2, mean)
                auxbar[idxA] <- apply(exp(auxParams[,idxA]), 2, mean)
            }
        }
        
        ## effective number of parameters (for use in calculating DIC)
        pD <- Dbar - Dhat

        ## Deviance Information Criterion (DIC) = Dhat + 2pD = Dhat + (Dbar - Dhat) + pD = Dbar + pD
        DIC <- Dbar + pD

        ## summary table
        sumtab <- data.frame(y = yall, n = nall, pobs = yall / nall, yexpbar = yhatbar, thetabar = pbar, Dbar = Dbar, Dhat = Dhat, pD = pD, DIC = DIC, row.names = colnames(dev))
        sumtab <- rbind(sumtab, c(NA,NA,NA,NA,NA, sum(Dbar), sum(Dhat), sum(pD), sum(DIC)))
        rownames(sumtab)[dim(sumtab)[1]] <- "TOTAL"
        sumtab
    }

