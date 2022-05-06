ancova.hsd <- function(TG1,TG2)
{
  stor.hsd = list()
  
  g1 = TG1; g2 = TG2; g = rbind(g1,g2); #g$Age = as.factor(g$Age); g$Weight = as.factor(g$Weight); g$Height = as.factor(g$Height)
  Y = c(rep('1',1,length(g1[,1])), rep('2',1,length(g2[,1])) ) 
  DF  = as.data.frame(cbind(as.factor(Y),g)); #DF[is.na(DF)] <- 0;
  
  gg=1;
  for(tt in 2:(ncol(DF))){
    # print(tt)
    #res.aov <- aov(DF[,tt]  ~ Y + Age + Weight + Height, data = DF)
    #res.aov <- aov(DF[,tt]  ~ Y + Age, data = DF)
    res.aov <- aov(DF[,tt]  ~ Y, data = DF)
    res.hsd = TukeyHSD(res.aov, which = 'Y', conf.level=0.95)
    
    stor.hsd[[gg]] =  res.hsd$`Y`[4];
    gg=gg+1;
  }
  
  fc = colMeans(data.matrix(g1)) - colMeans(data.matrix(g2)); #fc = fc[-c(which(names(fc) %in% c('Age','Weight','Height','VO2max')))]
  return(list(fc=fc, p.adj = unlist(stor.hsd)));
  
}
