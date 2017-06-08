library(R2OpenBUGS)


##################################################
### QUESTION 1:  BINARY MISCLASSIFICATION

## [ Code in models/cervix_mod.txt ]

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

### Just run on both datasets, compare the results and discuss.



##################################################
### QUESTION 2:  CLASSICAL MEASUREMENT ERROR (dugongs)

## [ Code in models/dugongs_mod.txt ]

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

zmedian <- dugongs.sim$median$z ## posterior median of true age
plot(zmedian, dat$y, xlab="Observed age (years)",
     ylab="Length (m)", pch=19, col="blue", xlim=c(0,35), ylim=c(1.5, 2.75))

## Exercise:
## Increase the prior measurement error SD (i.e. lower precision),
## run the model again, extracting the posterior median true ages,
## and overlay on the previous plot, something like this:
# points(zmedian, dat$y, col="red")



##################################################
### QUESTION 3:  BERKSON MEASUREMENT ERROR (pollution)

### Code in models/poll_mod.txt, as in the lectures, without calibration data

poll_data <- list(y = c(21, 20, 15), n = c(48, 34, 21),
                  x = c(10, 30, 50), alpha = 4.48, beta = 0.76)
poll_inits <- list(list(theta = c(0.0, 0.0), z = c(0.0, 0.0, 0.0)))
poll.sim <- bugs(data=poll_data, inits=poll_inits,
                 model = "poll_mod.txt",
                 working.directory="models",
                 param = c("theta","or10"),
                 n.burnin=1000, n.iter=100000, n.chains=1)
poll.sim

## Including calibration data relating observed exposures xcal to true exposures zcal
poll_data <- list(y = c(21, 20, 15), n = c(48, 34, 21),
                  x = c(10, 30, 50),
                  xcal = c(3, 6, 9, 10, 11, 13, 15, 20, 20.5, 21,
                           22, 23, 24, 24.5, 28, 30, 33, 41, 47, 53, 60, 70, 90),
                  zcal = c(1.1, 10.7, 3.8, 26.4, 15.8, 7, 20.3, 26.3, 25.2, 17.7, 34.8, 25.5,
                           17.1, 3.2, 35.9, 26.9, 29.4, 44.1, 47.6, 50.1, 58.4, 64.7, 73.6),
                  ncal = 23
                  )

## Exercise: adapt the model to include this data
