## Install packages ##
library(metafor)
library(ggmcmc)
library(devtools)
library(rstan)
library(lme4)
library(devtools)
library(shinystan)
library(glmer2stan)
library(rstan)
library(lmerTest)
library(haven)

## Load the Data 

IPDMA <- read_sas("~/Desktop/Statistics/My-Complete-Book-In-R/Bayesian Analysis/Multilevel/Maroeska/IPDMA.sas7bdat")
head(IPDMA)


IPDMA$OUTCOME =(IPDMA$POUTCOME -1)^2
IPDMA$AGED = (IPDMA$AGED -1 )^2

#### Bayesian Approach 


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


ageGLM = glm2stan(formula = OUTCOME~ TREAT*AGED , data= IPDMA,
                     family = "binomial", varpriors = "weak",
                     iter=1000,warmup = 100, chains=4,
                     calcWAIC=T)

summary(ageGLM)

fit= glm(formula = OUTCOME~TREAT*AGED , data= IPDMA,
         family = binomial)

summary(fit)


sum(ageGLMM@sim$samples[[1]]$beta_TREAT_X_AGED[1:1000]>0)/1000




ageGLMM = glmer2stan(formula = OUTCOME~ TREAT*AGED + (1|STUDY) , data= IPDMA,
           family = "binomial",
           iter=1000,warmup = 100, chains=4,
           calcWAIC=T)



summary(ageGLMM)[1]$summary[1:4,1:3]
summary(fit2)

fit2= glmer(formula = OUTCOME~TREAT*AGED + (1|STUDY), data= IPDMA,
         family = binomial)

summary(fit2)




a="data{
   int N;
   int OUTCOME[N];
   real TREAT[N];
   real AGED[N];
   int STUDY[N];
   real TREAT_X_AGED[N];
   int bin_total[N];
   int N_STUDY;
   }
   
   parameters{
   real Intercept;
   real beta_TREAT;
   real beta_AGED;
   real beta_TREAT_X_AGED;
   real vary_STUDY[N_STUDY];
   real<lower=0> sigma_STUDY;
   }
   
   model{
   real vary[N];
   real glm[N];
   // Priors
   Intercept ~ double_exponential_log( 0.9 , 1 );
   beta_TREAT ~ double_exponential_log( 0.71 , 1);
   beta_AGED ~ double_exponential_log( 0.64 , 1 );
   beta_TREAT_X_AGED ~ double_exponential_log( 0 , 1 );
   sigma_STUDY ~ gamma( 2 , 1e-4 );
   // Varying effects
   for ( j in 1:N_STUDY ) vary_STUDY[j] ~ normal( 0 , sigma_STUDY );
   // Fixed effects
   for ( i in 1:N ) {
   vary[i] <- vary_STUDY[STUDY[i]];
   glm[i] <- vary[i] + Intercept
   + beta_TREAT * TREAT[i]
   + beta_AGED * AGED[i]
   + beta_TREAT_X_AGED * TREAT_X_AGED[i];
   glm[i] <- inv_logit( glm[i] );
   }
   OUTCOME ~ binomial( bin_total , glm );
   }
   
   generated quantities{
   real dev;
   real vary[N];
   real glm[N];
   dev <- 0;
   for ( i in 1:N ) {
   vary[i] <- vary_STUDY[STUDY[i]];
   glm[i] <- vary[i] + Intercept
   + beta_TREAT * TREAT[i]
   + beta_AGED * AGED[i]
   + beta_TREAT_X_AGED * TREAT_X_AGED[i];
   dev <- dev + (-2) * binomial_log( OUTCOME[i] , bin_total[i] , inv_logit(glm[i]) );
   }
   }
   
   "


myageGLMM =  myglmer2stan(formula = OUTCOME~ TREAT*AGED + (1|STUDY) , data= IPDMA,
                          family = "binomial",
                          iter=1000,warmup = 100, chains=4,
                          calcWAIC=T, mymodel = a)


summary(myageGLMM)








shinystan::launch_shinystan(myageGLMM)

(pm_coef = colMeans(ageGLMM))



a=summary(ageGLMM)
a$summary[1:4,1:3]
stan_plot(ageGLMM)
stan_trace(ageGLMM)
stan_hist(ageGLMM)
stan_dens(ageGLMM)
stan_diag(ageGLMM)
stan_rhat(ageGLMM)
stan_ess(ageGLMM)
stan_mcse(ageGLMM)
stan_ac(ageGLMM)

rstan_ggtheme_options(panel.background = element_rect(fill = "white"),
                      legend.position = "center")
rstan_gg_options(fill = "black", color = "black", pt_color = "red")




traceplot(ageGLMM2)


ageGLMM2$model

summary(ageGLMM2, waic = TRUE)



### Frequentist approach

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



### 2x2 Matrix Age Interaction

fit= glm(formula = OUTCOME~TREAT*BILAT_0, data= IPDMA,
         family = binomial)

summary(fit)
exp(-coef(fit))

#### Bayesian Approach 

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


BILAT_GLM = glm2stan(formula = OUTCOME~ TREAT*BILAT_0 , data= IPDMA,
                  family = "binomial", varpriors = "weak",
                  iter=1000,warmup = 100, chains=4,
                  calcWAIC=T)

summary(BILAT_GLM)



sum(ageGLMM@sim$samples[[1]]$beta_TREAT_X_AGED[1:1000]>0)/1000


### 2x2 Matrix Age Interaction with random intercept

fit2= glmer(formula = OUTCOME~TREAT*BILAT_0 + (1|STUDY), data= IPDMA,
         family = binomial)

summary(fit2)
exp(-coef(fit2))


#### Bayesian Approach 
BILAT_GLMM = glmer2stan(formula = OUTCOME~ TREAT*BILAT_0 + (1|STUDY) , data= IPDMA,
                     family = "binomial",
                     iter=1000,warmup = 100, chains=4,
                     calcWAIC=T)



summary(BILAT_GLMM)[1]$summary[1:4,1:3]



### 2x2 Matrix Age Interaction with random slope

fit3= glmer(formula = OUTCOME~TREAT*BILAT_0+ (BILAT_0|STUDY), data= IPDMA,
            family = binomial)

summary(fit3)
exp(coef(fit3)$STUDY)



#### Bayesian Approach 
BILAT_GLMM2 = glmer2stan(formula = OUTCOME~ TREAT*BILAT_0 + (BILAT_0|STUDY) , data= IPDMA,
                        family = "binomial",
                        iter=1000,warmup = 100, chains=4,
                        calcWAIC=T)



summary(BILAT_GLMM2)[1]$summary[1:4,1:3]

