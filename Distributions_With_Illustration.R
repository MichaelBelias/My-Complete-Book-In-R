#### Distributions Illustrated

# Display the Student's t distributions with various
# degrees of freedom and compare to the normal distribution

x <- seq(-4, 4, length=100)
hx <- dnorm(x)

Degrees_Of_Freedom <- c(1, 3, 8, 30)
colors <- c("red", "blue", "darkgreen", "gold", "black")
labels <- c("df=1", "df=3", "df=8", "df=30", "normal")

plot(x, hx, type="l", lty=2, xlab="x value",
     ylab="Density", main="Comparison of t Distributions")

for (i in 1:4){
  lines(x, dt(x,Degrees_Of_Freedom[i]), lwd=2, col=colors[i])
}

legend("topright", inset=.05, title="Distributions",cex = 0.75,
       labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)


# normal distribution with a mean of 100 and a standard deviation of 15. 
# What proportion is between 80 and 120?

mean=100; sd=15 ; lb=80; ub=120

x <- seq(-4,4,length=100)*sd + mean ## Vector of x values
hx <- dnorm(x,mean,sd)

plot(x, hx, type="n", xlab="Values", ylab="",
     main="Normal Distribution", axes=FALSE)

i <- x >= lb & x <= ub
lines(x, hx)
polygon(c(lb,x[i],ub), c(0,hx[i],0), col="red") 

area <- pnorm(ub, mean, sd) - pnorm(lb, mean, sd)
result <- paste("P(",lb,"< Values <",ub,") =",
                signif(area, digits=3))
mtext(result,3)
axis(1, at=seq(40, 160, 20), pos=0)

