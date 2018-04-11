## Install packages ##
library(metafor)

## Load the Data 
library(haven)
IPDMA <- read_sas("/home/mike/Documents/My-Complete-Book-In-R/Bayesian Analysis/Multilevel/Maroeska/IPDMA.sas7bdat")
head(IPDMA)


IPDMA$OUTCOME =(IPDMA$POUTCOME -1)^2
IPDMA$AGED = (IPDMA$AGED -1 )^2



summary(IPDMA)


Baseline=matrix(nrow=33
                ,ncol = 3
                ,dimnames = list(names(IPDMA), c("Control","Treatment","Total")))



Baseline[22,]= c(table(IPDMA[which(IPDMA$AGE<2 ),]$TREAT),sum(table(IPDMA[which(IPDMA$AGE<2 ),]$TREAT))); row.names(Baseline)[22] = "Age <2"


for(i in c(1:3,5:21,23:33)){
  
  Baseline[i,]= c(table(IPDMA[which(IPDMA[,i]==1 ),]$TREAT),
                  sum(table(IPDMA[which(IPDMA[,i]==1 ),]$TREAT)))
}
Baseline=Baseline[-4,];rm(i)


write.csv(Baseline,file = "Tested Data/Baseline.csv")


### 2x2 Matrix Age<2

a=matrix(c(table(IPDMA[which(IPDMA$AGE<2 & IPDMA$POUTCOME==0),]$TREAT),table(IPDMA[which(IPDMA$AGE<2 & IPDMA$POUTCOME==1),]$TREAT)),nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"
a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

fit= glm(formula = OUTCOME~TREAT, data= IPDMA[which(IPDMA$AGE<2),],
         family = poisson)

exp(-coef(fit))

### 2x2 Matrix Age>=2

a=matrix(c(table(IPDMA[which(IPDMA$AGE>=2 & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$AGE>=2 & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

fit= glm(formula = OUTCOME~TREAT, data= IPDMA[which(IPDMA$AGE>=2),],
         family = binomial("log"))

exp(-coef(fit))
summary(fit)


### 2x2 Matrix Age Interaction

fit= glm(formula = OUTCOME~TREAT*AGED , data= IPDMA,
         family = binomial)

summary(fit)
exp(-coef(fit))


### 2x2 Matrix Bilateral AOM No

a=matrix(c(table(IPDMA[which(IPDMA$BILAT_0==0 & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$BILAT_0==0 & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

fit= glm(formula = OUTCOME~TREAT, data= IPDMA[which(IPDMA$BILAT_0 == 0),], family = poisson)

exp(-coef(fit))

### 2x2 Matrix Bilateral AOM Yes

a=matrix(c(table(IPDMA[which(IPDMA$BILAT_0==1 & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$BILAT_0==1 & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

fit= glm(formula = OUTCOME~TREAT, data= IPDMA[which(IPDMA$BILAT_0 == 1),], family = poisson)

exp(-coef(fit))


### 2x2 Matrix Bilateral AOM Yes and Age<2

a=matrix(c(table(IPDMA[which(IPDMA$BILAT_0==1 &IPDMA$AGE<2 & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$BILAT_0==1 &IPDMA$AGE<2 & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

### 2x2 Matrix Bilateral AOM No and Age<2

a=matrix(c(table(IPDMA[which(IPDMA$BILAT_0==0 &IPDMA$AGE<2 & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$BILAT_0==0 &IPDMA$AGE<2 & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

### 2x2 Matrix Bilateral AOM Yes and Age>=2

a=matrix(c(table(IPDMA[which(IPDMA$BILAT_0==1 &IPDMA$AGE>=2 & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$BILAT_0==1 &IPDMA$AGE>=2 & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

### 2x2 Matrix Bilateral AOM Yes and Age>=2

a=matrix(c(table(IPDMA[which(IPDMA$BILAT_0==0 &IPDMA$AGE>=2 & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$BILAT_0==0 &IPDMA$AGE>=2 & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))



#### Interaction
fit= glm(formula = OUTCOME~TREAT *AGED* BILAT_0 , data= IPDMA, family = binomial)
### ??
summary(fit)
###

### 2x2 Matrix Bilateral Otorrhea Yes

a=matrix(c(table(IPDMA[which(IPDMA$OTO_0==1  & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$OTO_0==1  & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))


fit= glm(formula = OUTCOME~TREAT, 
    data= IPDMA[which(IPDMA$OTO_0 == 1),], 
    family = poisson)

exp(-coef(fit))


### 2x2 Matrix Bilateral Otorrhea No

a=matrix(c(table(IPDMA[which(IPDMA$OTO_0==0  & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$OTO_0==0  & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

fit= glm(formula = OUTCOME~TREAT, 
         data= IPDMA[which(IPDMA$OTO_0 == 0),], 
         family = poisson)
summary()
exp(-coef(fit))


### Pain at 3–7 days
### Pain at 3–7 days 2x2 Matrix Age<2

a=matrix(c(table(IPDMA[which(IPDMA$AGE<2 & IPDMA$PAIN_2==0),]$TREAT),table(IPDMA[which(IPDMA$AGE<2 & IPDMA$PAIN_2==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"
a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

fit= glm(formula = OUTCOME~TREAT, data= IPDMA[which(IPDMA$AGE<2),],
         family = poisson)

exp(-coef(fit))

### Pain at 3–7 days 2x2 Matrix Age>=2

a=matrix(c(table(IPDMA[which(IPDMA$AGE>=2 & IPDMA$POUTCOME==0),]$TREAT),
           table(IPDMA[which(IPDMA$AGE>=2 & IPDMA$POUTCOME==1),]$TREAT)),
         nrow=2,ncol=2, dimnames = list(c("Control","Treatment"),
                                        c("No-Event","Event")))
a=cbind(a,apply(a,1,sum));colnames(a)[3]="Total"

a=rbind(a,apply(a,2,sum));row.names(a)[3]="Total"


dat=escalc("RR",  ai=a[1,1], bi= a[1,2], ci= a[2,1], di= a[2,2])

exp(summary(dat))

fit= glm(formula = OUTCOME~TREAT, data= IPDMA[which(IPDMA$AGE>=2),],
         family = binomial("log"))

exp(-coef(fit))
summary(fit)


### Pain at 3–7 days 2x2 Matrix Age Interaction

fit= glm(formula = PAIN_2~TREAT*AGED , data= IPDMA,
         family = binomial)

summary(fit)
exp(-coef(fit))


### Pain at 3–7 days 2x2 Matrix Age Interaction

fit= glm(formula = PAIN_2~TREAT*BILAT_0 , data= IPDMA,
         family = binomial)

summary(fit)
exp(-coef(fit))


### Pain at 3–7 days 2x2 Matrix Age Interaction

fit= glm(formula = PAIN_2~TREAT* BILAT_0 + TREAT* AGED+ I(BILAT_0*TREAT* AGED)   , data= IPDMA,
         family = binomial)

summary(fit)
exp(-coef(fit))

