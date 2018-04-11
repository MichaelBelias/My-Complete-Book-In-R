## ------------------------------------------------------------------------
library(pander)

df<- read.table("F:/My Complete Book In R/K-means/k-means.txt")

pander(head(df,10), caption = "10 first observations in our Data-Set" , 
       split.table = 60)

## ------------------------------------------------------------------------
df1 <- scale(df[-8]) 
wssplot <- function(data, nc=15, seed=697){
               wss <- (nrow(data)-1)*sum(apply(data,2,var))
               for (i in 2:nc){
                    set.seed(seed)
                    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
                plot(1:nc, wss, type="b", xlab="Number of Clusters",
                     ylab="Within groups sum of squares")}


wssplot(df1,7)


## ----comment="", message=FALSE-------------------------------------------
library(NbClust)
set.seed(6970222)

nc <- NbClust(df1, min.nc=2, max.nc=7, method="kmeans")
pander(table(nc$Best.n[1,]), caption = "")


## ---- message=FALSE------------------------------------------------------

library(ggplot2)

ggplot(data=as.data.frame(table(nc$Best.n[1,])),aes(x=Var1, y=Freq))+ 
    geom_bar(stat= "identity",fill="black",colour="orange",alpha=0.75)+
    xlab("Number of Clusters") + ylab("Number of Criteria") +
  ggtitle("Number of Clusters Chosen \n by 26 Criteria")



## ----comment="", message=FALSE-------------------------------------------
set.seed(6970222)
df1 <- scale(df[-8]) 
kmeans.result <- kmeans(x = df[,-8],centers = 3)
table(df$species, kmeans.result$cluster)

## ----comment="", message=FALSE-------------------------------------------
library(scatterplot3d)

scatterplot3d(x = df$kernel.length,
  df$kernel.width,
  df$area,
  main = " The Original Separation \n Vs The K-means",
  xlab = "Area",
  ylab = "Kernel Length",
  zlab = "Perimeter",
  color= as.integer(df$species),
  pch = kmeans.result$cluster,angle = 60
  
)


library(flexclust)

randIndex(table(df$species, kmeans.result$cluster))

## ------------------------------------------------------------------------
set.seed(6970222) 
library(fpc)

pamk.result <- pamk(df1)
pander(table(pamk.result$pamobject$clustering, df$species) )

## ------------------------------------------------------------------------


layout(matrix(c(1, 2), 1, 2)) # 2 graphs per page

plot(pamk.result$pamobject)



## ---- echo=FALSE---------------------------------------------------------
library(ggplot2)

iris2= iris3

ggplot(iris2,aes(Sepal.Length, Sepal.Width)) + geom_point(color=kmeans.result$cluster)

plot(iris2[c("Sepal.Length", "Sepal.Width")], col = kmeans.result$cluster)

points(kmeans.result$centers[, c("Sepal.Length", "Sepal.Width")],

col = 1:3, pch = 8, cex = 2) # plot cluster centers

## ------------------------------------------------------------------------
library(cluster)

# group into 3 clusters

pam.result <- pam(df[,-8], 3)

pander(table(pam.result$clustering, df$species))


## ------------------------------------------------------------------------
layout(matrix(c(1, 2), 1, 2)) # 2 graphs per page

plot(pam.result)

