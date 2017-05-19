### make data

nTrials = 12
nSubs = 30

fake <- expand.grid(trial = 1:nTrials, ppt = 1:nSubs, J = c('J1', 'J2'), K = c('K1', 'K2'), stringsAsFactors = T)

B0 <- 1.5 # grand mean
B1 <- .5 # main effect of J
B2 <- -.3 # main effect of K
s <- rnorm(nSubs, 0, 1)

# there are two true main effects and no interaction
logitSuccess <- with(fake, B0 + c(B1, -B1)[J] + c(B2, -B2)[K] + s[ppt])

logistic <- function(x){
  1/(1 + exp(-x))
}

fake$y <- rbinom(n = length(logitSuccess), size = 1, prob = logistic(logitSuccess))

write.csv(fake, 'fakeData.csv')

with(aggregate(y ~ J + K, data = fake, FUN = mean), barplot(y, names.arg = paste(J,K, sep = '-')))

