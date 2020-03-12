#### Simulation study

y = c(rep(0,40),rep(1,60))

fit.logit = glm(y~1,family = binomial("logit"))
table(predict(fit.logit,type= "response"))

fit.log = glm(y~1,family = binomial("log"))
table(predict(fit.log,type= "response"))

fit.identity = glm(y~1,family = binomial("logit"))
table(predict(fit.identity,type= "response"))




df =  data.frame(X =  rep(0:1,each=5000), Treat =  rep(0:1,5000))

df$Y = NA
df[df$X == 0 & df$Treat==0,]$Y = c(rep(0,1250),rep(1,1250))
df[df$X == 0 & df$Treat==1,]$Y =  c(rep(0,2000),rep(1,500))

df[df$X == 1 & df$Treat==0,]$Y = c(rep(0,500),rep(1,2000))
df[df$X == 1 & df$Treat==1,]$Y =  c(rep(0,2000),rep(1,500))


fit.logit = glm(Y~X*Treat,family = binomial("logit"),data = df)
table(round(predict(fit.logit,type= "response"),2))

fit.log = glm(Y~X*Treat,family = binomial("log"),data = df)
table(round(predict(fit.log,type= "response"),2))

fit.identity = glm(Y~X*Treat,family = binomial("identity"),data = df)
table(round(predict(fit.identity,type= "response"),2))

library(dplyr)
df%>%
  group_by(X,Treat)%>%
  summarise(Prop =  mean(Y))




fit.log = glm(y~1,family = binomial("log"))
table(predict(fit.log,type= "response"))




library(ggpubr)



df =  data.frame(X =  seq(0,1,length.out =10000), Treat =  rbinom(10000, 1, 0.5))

df$Y =  with(df, 0.3 +  X*0.15 + 0.3*Treat )




df$Treat =  as.factor(df$Treat)
df$X =  as.numeric(df$X)


ggscatter(df, "X","Y",color = "Treat")


max(df$Y); min(df$Y)




df$Y_binary =  rbinom(10000,1,prob = df$Y)

ggscatter(df, "X","Y_binary",color = "Treat", add = "loess")




fit.logit = glm(Y_binary~ X+Treat, family = binomial("logit"), data = df)


pred_fit.logit =  predict(fit.logit, type = "response")

library(mgcv)
fit.log = glm(Y_binary~ X*Treat, family = binomial("log"),  data = df)

pred_fit.log =  predict(fit.log, type = "response")

fit.identity= glm(Y_binary~  X+Treat, family = binomial("identity"), data = df)

pred_fit.identity =  predict(fit.identity, type = "response")


real.values =  df$Y



df = do.call("rbind", replicate(4, df[,c(1:2,4)], simplify = FALSE))
df$Prediction =  c(real.values,pred_fit.logit,pred_fit.log, pred_fit.identity)
df$Type =  rep(c("real.values","pred_fit.logit","pred_fit.log", "pred_fit.identity"),each=10000)
ggscatter(df, "X","Prediction", color = "Type", size = .5)






require(MASS)
summary(fit.NB <- glm.nb(Y_binary~  X+Treat,  data = df))

pred_fit.NB =  predict(fit.NB, type = "response")



table(df[df$Treat == 0,]$Y_binary)




