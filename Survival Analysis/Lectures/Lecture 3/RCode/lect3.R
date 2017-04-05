par(mfrow=c(1,1),ask=T)

subject<-unique(timedat$id)
for(i in subject){
  plot(timedat$fuptime[timedat$id==i],
       (fitted(cd4.lme)[timedat$id==i])^2,
       type="l", col="red", ylab=expression(paste("CD4 count (cells/",mu,"L)")), 
       xlab="Years on ART", 
       xlim=c(0,max(timedat$fuptime, na.rm=TRUE)), 
       ylim=c(0, (max(fitted(cd4.lme)^2, na.rm=TRUE))))
  points(timedat$fuptime[timedat$id==i], (timedat$cd4.sqrt[timedat$id==i])^2)
  text(4,30, pos=4, sprintf('Subject ID=%s',i))
}

basdat[]