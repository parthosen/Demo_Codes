PLSDA.Geno.Pairs <- function(DF1,Grps,grp.name1,grp.name2,grp.col,vip.cutoff)
{
  
  pdf(paste0("../Results/Simulations/figure/PLSDA_Geno/PLSDA_Genotype_Personalized_",grp.name1,'_',grp.name2,".pdf"), width = 30, height = 11, paper = 'special')
  oripar=par()$mar
  par(xpd=F, mar=oripar+c(2,18,0,6))
  par(mfrow=c(1,2))
  
  
    B = t( cbind(DF1[,which(colnames(DF1) %in% grp.name1)], DF1[,which(colnames(DF1) %in% grp.name2)]) ); rownames(B) <- NULL;

    Grps = Grps[which(Grps %in% c(grp.name1,grp.name2))]
    Grps[Grps==grp.name1] = 0;
    Grps[Grps==grp.name2] = 1;
    
    Grps = as.numeric(Grps);

    plsda <- plsda(B,as.factor(Grps), ncomp = 2);
    obj <- auroc(plsda, roc.comp = 1, plot = F,print = F) 
    
  
    Scores = plsda$variates$X;  rownames(Scores) = Grps;
    loads =  plsda$loadings$X
 
  ############################################################################
  ## % of expained variance ##
  if(1){
    print('Estimate the % of explained variance for each PLS component')
    Rd.YvsU = cor(as.numeric(as.factor(Grps)), plsda$variates$X[, 1:2])
    Rd.YvsU = apply(Rd.YvsU^2, 2, sum)
    Rd.Y = cbind(Rd.YvsU, cumsum(Rd.YvsU))
    colnames(Rd.Y) = c("Proportion", "Cumulative")
  }
  
  
  ellipseplot(Scores[,1],                              # data for x-axis
              Scores[,2],                              # data for y-axis
              as.factor(Grps),                        # factor with classes
              pcol= c(grp.col,'#85C1E9'), # '#D2B4DE','#82E0AA' colors for plotting (must match # of factors) ; '#F01F02','#04F1F9''#999999','#FFE4C4)
              pbgcol=TRUE,                           # point borders black?
              cexsize=3.5,                              # size of points
              ppch=rep(21,2),                         # shape of points (must match # of factors)
              legpos="topright",                      # position of legend
              legcexsize=2,                          # legend text size
              legptsize=3,                           # legend point size
              axissize=2,                            # Set axis text size
              linewidth=2,
              Samp.nom = NULL, #c(samp1,samp2,samp3,samp4)
              elev = 0.98,
              legs= as.factor(c(grp.name1, grp.name2)))
  
  title(xlab=paste0('PLS1 (', round(Rd.Y[1,1]*100,1)[1], '%)'),    # % variance explained on PLS1
        ylab=paste0('PLS2 (', round(Rd.Y[2,1]*100,1)[1], '%)'),    # % variance explained on PLS2
        main= paste0('PLSDA ','[AUC = ',round(obj$Comp1[1],2), ']'),                  # Title
        cex.lab=2,                    # size of label text
        cex.main=2                    # size of title text
  )
  
  ################### Next: Variable Selection by VIP's ################################
  
  print('Getting the VIPs and Regression coefficients from model objects')
  RCV = plsda$mat.c
  VIP = vip(plsda);
  
  ### getting all variables with more than 1 VIPs ####
  vip1 <- sort(VIP[,1]); vip1  <- vip1[vip1>=vip.cutoff];
  rcv <- sort(RCV[,1]); rcv  <- rcv[vip1>=vip.cutoff];
  
  cols = rev(mako(n=50,alpha = 0.8));
  bh = barplot(vip1, horiz = T, col = cols[1:length(vip1)] , las=2, xlim=c(0,max(vip1)+0.20), width = 0.1, space = 0.3, xlab = 'VIPs', main = 'PLSDA:VIPs', las=1, cex.names = 1.2, cex.axis = 2, cex.main = 2, cex.lab=2, border = 'black' );
  box(lty = 1, lwd=2, col = 'black')
  
  dev.off()
  
  return(vip1=vip1)
  
  
}