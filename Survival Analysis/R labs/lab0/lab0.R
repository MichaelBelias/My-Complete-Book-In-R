################################
### LAB 0: Introduction to R ###
################################

# To see the help file of some function
?sqrt

# Basic algebra in R
6 + 2
9 - 3
3*5
6/4
sqrt(9)
log(1)
exp(log(1))
exp(log(3))

# Generating number sequences
1:10
seq(1,9,by = 2)
rep(3,5)
x = 1:3
x
rep(x,each = 5)

# Vectors
x = c(8,11,23,2,9,6,21,17,31)
x
x[5]
x[1:5]
x[c(2,3,8)]

# Some important functions defined for vectors
min(x)
max(x)
c(min(x),max(x))
range(x)
length(x)

# To sort a vector
sort(x)
x[order(x)]

# Character vector - string
vec.c = c("Athens","London","NY","Paris")
vec.c

# Concatenate strings
paste(vec.c,collapse = ";")
y = paste(vec.c,1:4)
y

# logical operators in R
x = 1:10
x

x > 5
x > 5 & x < 8

#see
x[x > 5]
x[x > 5 & x < 8]

x == 5
x != 5

# Basic matrix algebra
z = 1:20
matrix(z,nrow = 4,ncol = 5)
matrix(z,nrow = 4,ncol = 5,byrow = T)

A = rbind(c(1,2,3),c(4,5,6))
A
A[,1]
A[,2]
A[1,]

B = cbind(c(1,3),c(2,4))
B
B[1,]
B[2,]

B%*%A      # matrix multiplication
B^2        # elements to the power
B%*%B      # matrix multiplication
t(A)       # transpose

ncol(A)
nrow(A)
c(nrow(A),ncol(A))
dim(A)

# Extract diagonal elements
x = 1:9
X = matrix(x,3,3)
X
diag(X)

diag(1:5)

C = matrix(c(16,7,3,21),2,2)
C
solve(C)   # inverse of C
solve(C) %*% C
C %*% solve(C)

# Data frames: Extremely important
# A data frame is essentially a dataset

a = c("Athens","London","NY","Paris")
a
b = c("Greece","UK","USA","France")
b
c = c(5,8,19,10)
c
geo = data.frame(country = b,city = a,population = c)
geo

# See how we get the variables as vectors
geo$country
geo$city
geo$population

# Load data from R
data(esoph)
?esoph
esoph[1:10,]

# Create one-way and two-way frequency tables
table(esoph$agegp)
table(esoph$agegp,esoph$alcgp)

# Import library
library(survival)
?aml
aml

# Mean time by group
tapply(aml$time,aml$x,mean)

mean(aml$time[aml$x == "Maintained"])
mean(aml$time[aml$x == "Nonmaintained"])

# Median times
tapply(aml$time,aml$x,median)

# Just for fun, fit a linear model
lmfit = lm(time ~ x,data = aml)
summary(lmfit)

coef(lmfit)

# Plotting ...
?pbc
pbc[1:5,]

plot(pbc$age,pbc$chol,xlab = "Age in years",ylab = "Serum cholesterol (mg/dl)")

# or
plot(chol ~ age,xlab = "Age in years",ylab = "Serum cholesterol (mg/dl)",data = pbc)

