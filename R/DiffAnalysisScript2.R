##########################################################################################
# EPoS project
# Partho Sen: Jan 12 2019
# RNAseq: This Script reads MetaData files and gets the Counts file to make a common Table
# 64253 transcripts standarized #
# T1=cbind(T1,0); T2=t(apply(T1, 1, function(x) tapply(x, colnames(T1), mean)))
# DEG both by DSEQ and ttest #
# Normalized gene expression data: write to a File #
###########################################################################################
library('nlme') ## gives lme & nlme
library('lme4') ## gives lmer
library('gplots')
library('ggplot2')
library('ade4')
library('corrplot')
library('reshape')
library('Hmisc') ## Also using for imputed fumction
#library(ChemometricsWithR)
#library(ChemometricsWithRData)
library(plyr)
library(car)
library(maptools)
library(rgeos)
library('plotrix')
library(ade4)
library(factoextra)
library(magrittr)
library('DESeq2')
#library('DESeq')


rm(list=ls())
setwd('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Results/')

DFF  = read.table('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Data/RNAseq/cohort_rnaseq_pheno.csv',sep='\t', header = TRUE,check.names=FALSE) 
DFF  = DFF [1:nrow(DFF),]; ##Check rows

Samp = DFF[,1];
FN = DFF[,2];
SPICD = DFF[,3];
Grp = DFF[,5];

load(file="RNAseqTable.rda") ## Load genes count table  log2 transformed ### processed in Script1 ## log2 transformed ##
#RNAseqTable <- RNAseqTable[rowSums(RNAseqTable==0)<5, ]; ## If atleast five sample is zero remove the transcript
#RNAseqTable <- RNAseqTable[rowSums(RNAseqTable)>25000, ]; ## ## If rowSum is less than hundred remove the transcript ##
#RNAseqTable[RNAseqTable == 0] <- NA
RNAseqTable[is.infinite(RNAseqTable)] <- 0.0002; ## Log transformed!
#RNAseqTable[is.na(RNAseqTable)] <- mean(data.matrix(RNAseqTable),1, na.rm = T)
#RNAseqTable = impute(RNAseqTable,fun=minimum);

if(TRUE){
 write.table(RNAseqTable,'../Data/RNAseqTable_Alldata_log2.csv',sep = ',')
}
stop()


id.cirr = which(Grp %in% 'cirrhosis')
id.ctrl  = which(Grp %in% 'control')
id.nas  = which(Grp %in% 'NASH')
id.stea = which(Grp %in% 'Steatosis')

cirr.genes = RNAseqTable[,which(colnames(RNAseqTable) %in% FN[id.cirr])]; colnames(cirr.genes) <- NULL 
ctrl.genes = RNAseqTable[,which(colnames(RNAseqTable) %in% FN[id.ctrl])]; colnames(ctrl.genes) <- NULL 
nas.genes   = RNAseqTable[,which(colnames(RNAseqTable) %in% FN[id.nas])]; colnames(nas.genes) <- NULL 
stea.genes  = RNAseqTable[,which(colnames(RNAseqTable) %in% FN[id.stea])]; colnames(stea.genes) <- NULL 

### ******** QC plots *************** ###

pdf("./Figure/QCplot.pdf", width = 15, height = 8)
oripar=par()$mar
par(xpd=T, mar=oripar+c(2,2,0,6)) ## c(down,left,up,right) take some area to right # XPD = T plots outside
par(mfrow=c(2,2))

boxplot(ctrl.genes,las=2,notch=T,log = "",
        pars = list(boxwex = 0.8, staplewex = 1.5, outwex = 0.5, cex.axis = 0.5, whisklty = 2, staplelwd = 1, outcol="grey",border=c("black")),medcol="black",border = "black",horizontal = FALSE,add = FALSE,at = NULL,col = 'gray80',cex=0.3,ylab="log2(Intensity)",main="Control samples",xpd=T,ylim = c(0, 25));
boxplot(nas.genes,las=2,notch=T,log = "",
        pars = list(boxwex = 0.8, staplewex = 1.5, outwex = 0.5, cex.axis = 0.5, whisklty = 2, staplelwd = 1, outcol="grey",border=c("black")),medcol="black",border = "black",horizontal = FALSE,add = FALSE,at = NULL,col = 'gray80',cex=0.3,ylab="log2(Intensity)",main="NASH samples",xpd=T,ylim = c(0, 25));
boxplot(stea.genes,las=2,notch=T,log = "",
        pars = list(boxwex = 0.8, staplewex = 1.5, outwex = 0.5, cex.axis = 0.5, whisklty = 2, staplelwd = 1, outcol="grey",border=c("black")),medcol="black",border = "black",horizontal = FALSE,add = FALSE,at = NULL,col = 'gray80',cex=0.3,ylab="log2(Intensity)",main="Steatosis samples",xpd=T,ylim = c(0, 25));
boxplot(cirr.genes,las=2,notch=T,log = "",
        pars = list(boxwex = 0.8, staplewex = 1.5, outwex = 0.5, cex.axis = 0.5, whisklty = 2, staplelwd = 1, outcol="grey",border=c("black")),medcol="black",border = "black",horizontal = FALSE,add = FALSE,at = NULL,col = 'gray80',cex=0.3,ylab="log2(Intensity)",main="Cirrhosis samples",xpd=T,ylim = c(0, 25));

dev.off()


### ******** PCA *************** ###

colnames(nas.genes) <- rep('NASH',dim(nas.genes)[2]);
colnames(stea.genes) <- rep('Steatosis',dim(stea.genes)[2]);
colnames(cirr.genes) <- rep('Cirrhosis',dim(cirr.genes)[2]);
colnames(ctrl.genes) <- rep('Control',dim(ctrl.genes)[2]);

A = cbind(nas.genes,stea.genes,cirr.genes,ctrl.genes);
A2 = t(A);

source("/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_Celiac3_new/Analysis/range.R")
source("/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_Celiac3_new/Analysis/ellipseplot.R")
source("/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_Celiac3_new/Analysis/PCAvarAxis.R")

PCs <- prcomp(A2, scale = TRUE, center = TRUE); ## eigs <- pca$sdev^2

explainPCAvar <- PCAvarAxis(PCs)

pdf("./Figure/PCA.pdf", width = 12, height = 9,paper = 'special')
oripar=par()$mar
#par(xpd=F,mar=oripar+c(2,2,0,12)) ## c(down,left,up,right) take some area to right # XPD = T plots outside
par(xpd=T, mar=oripar+c(2,2,0,6))

ellipseplot(PCs$x[,1],                              # data for x-axis
            PCs$x[,2],                              # data for y-axis
            as.factor(rownames(PCs$x)),                        # factor with classes
            pcol=c("#FF008F","#00FF04","#0068FF","#E0EA13"),            # colors for plotting (must match # of factors) ; '#F01F02','#04F1F9''#999999','#FFE4C4)
            pbgcol=FALSE,                           # point borders black?
            cexsize=2,                              # size of points 
            ppch=rep(21,4),                         # shape of points (must match # of factors)
            legpos="topright",                      # position of legend           
            legcexsize=1,                          # legend text size
            legptsize=2,                           # legend point size 
            axissize=1,                            # Set axis text size
            linewidth=2                            # Set axis line size
)  
#draw.ellipse(PCs$x[,1],PCs$x[,2])
grid(lwd=1)

#text(PCs$x[,1],                              # data for x-axis
#PCs$x[,2],
#labels = rownames(A),
#cex = 0.6)

title(xlab=explainPCAvar[["PC1"]],    # % variance explained on PC1
      ylab=explainPCAvar[["PC2"]],    # % variance explained on PC2 
      main="PCA",                  # Title
      cex.lab=2,                    # size of label text
      cex.main=2                    # size of title text
)
#summary(PCs)# it has variance info
dev.off()

################################## PLS-DA SECTION #################################

################################# 'ropls package' #######################
########   ****** PLS-DA plots without averaging  *****  #########

if(FALSE){ ## It takes a while ##

library('ropls')

grps = rownames(A);

plsda <- opls(A,grps,
              predI=4,
              algoC = c("nipals"),
              scaleC = c("standard"),
              plotL=FALSE,
              crossvalI=7,
              permI = 20
)

Scores = plsda@scoreMN;
loads = plsda@loadingMN;

pdf("./Figure/PLSDA_Scores.pdf", width = 12, height = 9,paper = 'special')
oripar=par()$mar
par(xpd=F,mar=oripar+c(2,2,0,12)) ## c(down,left,up,right) take some area to right # XPD = T plots outside

ellipseplot(Scores[,1],                              # data for x-axis
            Scores[,2],                              # data for y-axis
            as.factor(rownames(Scores)),                        # factor with classes
            pcol=c("#FF008F","#00FF04","#0068FF","#E0EA13"),            # colors for plotting (must match # of factors) ; '#F01F02','#04F1F9''#999999','#FFE4C4)
            pbgcol=FALSE,                           # point borders black?
            cexsize=2,                              # size of points 
            ppch=rep(21,4),                         # shape of points (must match # of factors)
            legpos="topright",                      # position of legend           
            legcexsize=1,                          # legend text size
            legptsize=2,                           # legend point size 
            axissize=1,                            # Set axis text size
            linewidth=2                           # Set axis line size
)  
#draw.ellipse(PCs$x[,1],PCs$x[,2])
grid(lwd=1)


title(xlab='t1',    # % variance explained on PC1
      ylab='t2',    # % variance explained on PC2 
      main="PLS-DA",                  # Title
      cex.lab=2,                    # size of label text
      cex.main=2                    # size of title text
)
#summary(PCs)# it has variance info
dev.off()

}


########################################## **** Differential expression analysis **** ##########################################


####  ********** by ttest ************************ ##################

source('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Scripts/diffExp.R')
NC <- diffExp(nas.genes,ctrl.genes); DF.NC = data.frame(log2FoldChange=NC$FC,  pvalue=NC$pval,padj=NC$pval.adj); rownames(DF.NC) <-  rownames(nas.genes);
SC <- diffExp(stea.genes,ctrl.genes); DF.SC = data.frame(log2FoldChange=SC$FC, pvalue=SC$pval,padj=SC$pval.adj); rownames(DF.SC) <- rownames(stea.genes);
CC <- diffExp(cirr.genes,ctrl.genes); DF.CC = data.frame(log2FoldChange=CC$FC, pvalue=CC$pval,padj=CC$pval.adj); rownames(DF.CC) <- rownames(cirr.genes);

source('/Users/pasase/OneDrive/Projects_BTK_SysMed_May2017/Partho_NFALD/Scripts/volcano2.R')

de_threshold=0.0001;


volcano2(DF.NC,"NASH Vs Control",0,de_threshold,"./Figure/NASH_Vs_Control_ttest.pdf",c(-10,10))
volcano2(DF.SC,"Steatosis Vs Control",0,de_threshold,"./Figure/Steatosis_Vs_Control_ttest.pdf",c(-10,10))
volcano2(DF.CC,"Cirrhosis Vs Control",0,de_threshold,"./Figure/Cirrhosis_Vs_Control_ttest.pdf",c(-10,10))


#NC1 <- subset(DF.NC, padj < 0.05 & abs(log2FoldChange) > 0)
#SC1 <- subset(DF.SC, padj < 0.05 & abs(log2FoldChange) > 0)
#CC1 <- subset(DF.CC, padj < 0.05 & abs(log2FoldChange) > 0)

#write.table(data.matrix(NC1),"NASH_DEG_ttest.csv",sep=',', row.names = TRUE, col.names = NA)
#write.table(data.matrix(SC1),"Steatosis_DEG_ttest.csv",sep=',', row.names = TRUE, col.names = NA)
#write.table(data.matrix(CC1),"Cirrhosis_DEG_ttest.csv",sep=',', row.names = TRUE, col.names = NA)



####  ********** by DESEQ not DSEQ2  [no replicates]************************ ##################


Dseq.NC  = read.table('NASH_Ctrl_DEG.csv',sep=',', header = TRUE,check.names=FALSE); rownames(Dseq.NC) <- Dseq.NC[,1]; Dseq.NC = Dseq.NC[,-1];
Dseq.SC  = read.table('Steatosis_Ctrl_DEG.csv',sep=',', header = TRUE,check.names=FALSE); rownames(Dseq.SC) <- Dseq.SC[,1]; Dseq.SC = Dseq.SC[,-1];
Dseq.CC  = read.table('Cirrhosis_Ctrl_DEG.csv',sep=',', header = TRUE,check.names=FALSE); rownames(Dseq.CC) <- Dseq.CC[,1]; Dseq.CC = Dseq.CC[,-1]

de_threshold=0.0001;

volcano2(Dseq.NC,"NASH Vs Control",0,de_threshold,"./Figure/NASH_Vs_Control_DSEQ.pdf",c(-8,8))
volcano2(Dseq.SC,"Steatosis Vs Control",0,de_threshold,"./Figure/Steatosis_Vs_Control_DSEQ.pdf",c(-8,8))
volcano2(Dseq.CC,"Cirrhosis Vs Control",0,de_threshold,"./Figure/Cirrhosis_Vs_Control_DSEQ.pdf",c(-8,8))

NC2 <- subset(Dseq.NC, padj < de_threshold & abs(log2FoldChange) > 0)
SC2 <- subset(Dseq.SC, padj < de_threshold & abs(log2FoldChange) > 0)
CC2 <- subset(Dseq.CC, padj < de_threshold & abs(log2FoldChange) > 0)
#l1 = length(which(Dseq.NC[,6] < 0.05));

write.table(data.matrix(NC2),"NASH_DEG_Dseq0.05.csv",sep=',', row.names = TRUE, col.names = NA)
write.table(data.matrix(SC2),"Steatosis_DEG_Dseq0.05.csv",sep=',', row.names = TRUE, col.names = NA)
write.table(data.matrix(CC2),"Cirrhosis_DEG_Dseq0.05.csv",sep=',', row.names = TRUE, col.names = NA)

################ *Compare DEG among three different conditions* #################################

pdf("./Figure/Venn_Conditions.pdf", width = 9, height = 10)
  Comm <- venn(list('NASH vs Ctrl'= rownames(NC2), 'Steatosis vs Ctrl'=  rownames(SC2), 'Cirrhosis vs Ctrl'=  rownames(CC2)), small=0.7) ## upto 5 leafs
dev.off()

### ******* Here decided that we will use DESeq for the differential Expression analysis ******** ####################

#################################################################### *************************** The END *********************** #############################################################################
