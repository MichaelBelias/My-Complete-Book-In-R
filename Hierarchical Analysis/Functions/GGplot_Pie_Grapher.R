function(x){
  
  table(Myopia$MYOPIC)
  
  pie <- ggplot(x) + geom_bar(position = "stack")
  
  
  at <- as.numeric(cumsum(table(x))-0.5*table(x))
  labels <- paste(round(table(x)/ sum(table(x)),2 )*100 , "%" , sep = "")
  
  pie + coord_polar(theta = "y") + ylab("Frequencies") + labs(title = "Smoking percentage")  + 
    annotate(geom = "text", y = at , x = 1, label = labels ) + xlab(" ") + 
    scale_fill_manual(values= c("Green" , "Brown"),name  ="Smoking",labels=c("No", "Yes") ) 
  
}


x=Myopia$MYOPIC
