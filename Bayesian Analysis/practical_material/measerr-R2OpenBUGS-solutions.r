library(R2OpenBUGS)


##################################################
### QUESTION 1:  BINARY MISCLASSIFICATION

dat <- data.frame(
d = c(1,1,1,1,0,0,0,0,1,1,0,0),
z = c(0,0,1,1,0,0,1,1,NA,NA,NA,NA),
x = c(0,1,0,1,0,1,0,1,0,1,0,1)
)
count = c(13,3,5,18,33,11,16,16,318,375,701,535)

## full data
cervix_data <- c(as.list(dat[rep(1:12,count),]), n = sum(count))
## gold standard subset only
cervix_data <- c(as.list(dat[rep(1:8,count[1:8]),]), n = sum(count[1:8]))

cervix_inits <-
  list(list(psi = 0.5, beta0 = 0, beta = 0,
       phi = structure(.Data = c(0.5,0.5,0.5,0.5), .Dim = c(2,2))))

cervix.sim <- bugs(data=cervix_data, inits=cervix_inits,
                   model = "cervix_mod.txt",
                   working.directory = "models",
                   param = c("beta","or.beta","phi"),
                   n.burnin=1000, n.iter=20000, n.chains=1)
cervix.sim

##           mean  sd  2.5%   25%   50%   75% 97.5%
## beta       0.6 0.4  -0.1   0.4   0.6   0.9   1.4
## or.beta    2.0 0.8   0.9   1.4   1.9   2.4   4.1
## phi[1,1]   0.3 0.1   0.1   0.2   0.3   0.3   0.4
## phi[1,2]   0.2 0.1   0.1   0.2   0.2   0.3   0.4
## phi[2,1]   0.5 0.1   0.3   0.4   0.5   0.6   0.7
## phi[2,2]   0.8 0.1   0.6   0.7   0.8   0.8   0.9
## deviance 443.7 3.6 438.6 441.1 443.1 445.7 452.4

## With complete data only: the posterior precision of the OR is not much worse: 95% CI (0.9, 4.1).  Though the measurement error is big, and the precision improvements would have been greater if the measurement error were smaller.

### While we do get slightly more precise estimates from the full model, the potential disadvantages are the assumptions that it relies on.  Specifically, the assumption that the complete cases have the same parameter values as the incomplete cases.  Particularly the misclassification probabilities and odds ratios.   These assumptions are uncheckable from data.


##################################################
### QUESTION 2:  CLASSICAL MEASUREMENT ERROR (dugongs)

dat <- list(x = c(1.0,  1.5,  1.5,  1.5, 2.5,   4.0,  5.0,  5.0,  7.0,
                  8.0,  8.5,  9.0,  9.5, 9.5,  10.0, 12.0, 12.0, 13.0,
                  13.0, 14.5, 15.5, 15.5, 16.5, 17.0, 22.5, 29.0, 31.5),
            y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47,
                  2.19, 2.26, 2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43,
                  2.47, 2.56, 2.65, 2.47, 2.64, 2.56, 2.70, 2.72, 2.57),
            n = 27)

inits <- list(list(alpha = 3, beta = 2, gamma = 0.9, log.sigma = -5))

dugongs.sim <- bugs(data=dat, inits=inits, model="dugongs_mod.txt",
                    working.directory="models",
                    param = c("mu","z","alpha","beta","gamma","sigma2"),
                    n.burnin=1001, n.iter=10000, n.chains=1)


## Extract posterior median of true age

zmedian <- dugongs.sim$median$z
plot(zmedian, dat$y, xlab="Observed age (years)",
     ylab="Length (m)", pch=19, col="blue", xlim=c(0,35), ylim=c(1.5, 2.75))
# points(zmedian, dat$y, col="red") # after increasing the error...

## Solution:
## You should find that after increasing the measurement error, the dugong length is a smooth nonlinear function of the fitted true ages z, with less noise than in the original model.   The reason is that the estimates of z are now influenced much more by their likelihood contribution from the length data y through the nonlinear model, compared to the likelihood contribution from the observed ages x.   Increasing the measurement error makes the information provided by x about z weaker.   Also notice that the uncertainty around the fitted true ages increases when the measurement error is increased (compare the posterior standard deviations of z).



##################################################
### QUESTION 3:  BERKSON MEASUREMENT ERROR (pollution)

### As presented in the lectures, without calibration data

poll_data <- list(y = c(21, 20, 15), n = c(48, 34, 21),
                  x = c(10, 30, 50), alpha = 4.48, beta = 0.76)
poll_inits <- list(list(theta = c(0.0, 0.0), z = c(0.0, 0.0, 0.0)))
poll.sim <- bugs(data=poll_data, inits=poll_inits,
                 model = "poll_mod.txt",
                 working.directory="models",
                 param = c("theta","or10"),
                 n.burnin=1000, n.iter=100000, n.chains=1)
poll.sim

## Current: 1 chains, each with 1e+05 iterations (first 1000 discarded)
## Cumulative: n.sims = 99000 iterations saved
##          mean  sd 2.5%  25%  50%  75% 97.5%
## theta[1] -0.6 0.8 -2.2 -0.9 -0.5 -0.2   0.4
## theta[2]  0.0 0.0  0.0  0.0  0.0  0.0   0.1
## or10      1.5 0.5  1.0  1.2  1.4  1.6   2.6
## deviance 14.1 2.1 11.8 12.6 13.6 15.1  19.7

## Including calibration data relating observed exposures xcal to true exposures zcal
poll_data <- list(y = c(21, 20, 15), n = c(48, 34, 21),
                  x = c(10, 30, 50),
                  xcal = c(3, 6, 9, 10, 11, 13, 15, 20, 20.5, 21,
                           22, 23, 24, 24.5, 28, 30, 33, 41, 47, 53, 60, 70, 90),
                  zcal = c(1.1, 10.7, 3.8, 26.4, 15.8, 7, 20.3, 26.3, 25.2, 17.7, 34.8, 25.5,
                           17.1, 3.2, 35.9, 26.9, 29.4, 44.1, 47.6, 50.1, 58.4, 64.7, 73.6),
                  ncal = 23
                  )

## Solution: model including calibration data

poll_inits <- list(list(theta = c(0.0, 0.0), alpha=0.2, beta=1, sig=1))
poll.cal.sim <- bugs(data=poll_data, inits=poll_inits,
                     model = "poll_cal_mod_solution.txt",
                     working.directory="models",
                     param = c("theta","or10","alpha","beta","sig"),
                     n.burnin=1000, n.iter=100000, n.chains=1)
poll.cal.sim

## Current: 1 chains, each with 1e+05 iterations (first 1000 discarded)
## Cumulative: n.sims = 99000 iterations saved
##           mean  sd  2.5%   25%   50%   75% 97.5%
## theta[1]  -0.6 0.5  -1.8  -0.9  -0.6  -0.2   0.3
## theta[2]   0.0 0.0   0.0   0.0   0.0   0.0   0.1
## or10       1.4 0.3   1.0   1.2   1.4   1.5   2.2
## alpha      4.2 2.9  -1.6   2.3   4.3   6.1   9.9
## beta       0.8 0.1   0.7   0.8   0.8   0.9   1.0
## sig        8.2 1.3   6.1   7.3   8.0   9.0  11.3
## deviance 175.0 3.4 170.4 172.5 174.3 176.7 183.3

### Including the calibration data does not make much difference to the posterior distribution of the odds ratio, compared to using known error coefficients alpha,beta,sigma.  (which are similar to the posterior estimates of alpha, beta and sigma from the data).

## If anything, the estimates appear more precise with the calibration data included - it is unclear why - perhaps because the calibration model permits the chance of smaller (as well as larger) measurement errors, and smaller errors would lead to more precise estimates.

## Though the point of this exercise is just to get used to coding the evidence synthesis in BUGS!
