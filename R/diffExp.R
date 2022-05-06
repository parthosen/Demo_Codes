diffExp <- function(T0d,T0d1)
{
  
  pval = rep(1,dim(T0d1)[1]);
  FC = rep(0,dim(T0d1)[1]);
  
  
  for(i in 1:dim(T0d1)[1]){
    
    ##print(i)
    
    A = T0d[i,];
    B = T0d1[i,];
    
    LL = which(A %in% B);
    
    if(sum(A) != sum(B) &  length(LL)==0 & !all(A == A[1]) & !all(B == B[1]) ){
    
    #pval[i] = t.test(A,B[1:length(A)],paired = TRUE)$p.value; ## paired
    pval[i] = t.test(A,B,paired = FALSE, var.equal = FALSE)$p.value; ## unpaired
    
    ml = min(c(length(A),length(B)))
    FC[i] = mean(A[1:ml] - B[1:ml]);	
    
    }
    
  }
  
  pval.adj = p.adjust(pval, method = "fdr", n = length(pval))
  
  return(list(FC=FC, pval = pval, pval.adj=pval.adj));
  

}
  