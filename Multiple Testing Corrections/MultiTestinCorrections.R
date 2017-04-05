set.seed(6970222)


## Pick 2 mix of two distributions
x <- c(rnorm(900), rnorm(100, mean = 3))
## Take a random
p <- pnorm(x, lower.tail = F)


test <- p > 0.05
summary(test[1:900])
summary(test[901:1000])
plot(p)
hist(p)
bonftest <- p > 0.00005
summary(bonftest[1:900])
summary(bonftest[901:1000])

library(cherry)
data(NAEP)
